#ifndef RESERVOIRS_HH
#define RESERVOIRS_HH

#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct host {

  host(std::vector<std::string>& data, std::vector<std::string>& header,
       bool tbg = true, bool tbng = false) :
    M(0)
  {
    for (size_t i = 0; i < header.size(); ++i) {
      if (header[i] == "name") {
        name  = data[i];
      } else if (header[i] == "N") {
        std::istringstream s(data[i]);
        s >> N;
      } else if (header[i] == "pos_tbng") {
        if (tbng) {
          size_t posNG;
          std::istringstream s(data[i]);
          s >> posNG;
          M += posNG;
        }
      } else if (header[i] == "pos_tbg") {
        if (tbg) {
          size_t posG;
          std::istringstream s(data[i]);
          s >> posG;
          M += posG;
        }
      } else if (header[i] == "mortality") {
        std::istringstream s(data[i]);
        s >> mu;
      } else if (header[i] == "rec_rate") {
        std::istringstream s(data[i]);
        s >> gamma;
      } else if (header[i] == "abundance") {
        std::istringstream s(data[i]);
        s >> abundance;
      } else if (header[i] == "theta") {
        std::istringstream s(data[i]);
        s >> theta;
      }
    }
  }

  std::string name;
  size_t M, N;
  double mu, gamma, theta, abundance;
};

struct vector {

  vector(std::vector<std::string>& data, std::vector<std::string>& header,
         bool tbg = true, bool tbng = false) :
    M(0)
  {
    for (size_t i = 0; i < header.size(); ++i) {
       if (header[i] == "name") {
         name  = data[i];
       } else if (header[i] == "N") {
         std::istringstream s(data[i]);
         s >> N;
       } else if (header[i] == "pos_tbng") {
         if (tbng) {
           size_t posNG;
           std::istringstream s(data[i]);
           s >> posNG;
           M += posNG;
         }
       } else if (header[i] == "pos_tbg") {
         if (tbg) {
           size_t posG;
           std::istringstream s(data[i]);
           s >> posG;
           M += posG;
         }
       } else if (header[i] == "mortality") {
         std::istringstream s(data[i]);
         s >> mu;
       } else if (header[i] == "rec_rate") {
         std::istringstream s(data[i]);
         s >> gamma;
       } else if (header[i] == "density") {
         std::istringstream s(data[i]);
         s >> density;
       } else if (header[i] == "biting_rate") {
        std::istringstream s(data[i]);
        s >> bitingRate;
       }
    }
  }

  std::string name;
  size_t M, N;
  double prevalence, mu, gamma, density, bitingRate;
};

struct param {
  param(): areaConvert(1.) {}
  param(std::vector<std::string>& data, std::vector<std::string>& header)
  {
    for (size_t i = 0; i < header.size(); ++i) {
      if (header[i] == "area_convert") {
        std::istringstream s(data[i]);
        s >> areaConvert;
      }
    }
  }
  double areaConvert;
};

struct betafunc_params
{

  betafunc_params(std::vector<host> const &hosts,
                  std::vector<vector> const &vectors,
                  param const &params,
                  bool uvp = false,
                  bool bp = false) :
    hosts(hosts),
    vectors(vectors),
    params(params),
    hLambda(std::vector<double>(hosts.size())),
    hPrevalence(std::vector<double>(hosts.size())),
    vPrevalence(std::vector<double>(vectors.size())),
    useVectorPrevalence(uvp)
  {}

  std::vector<host> const &hosts;
  std::vector<vector> const &vectors;
  param const& params;
  
  std::vector<double> hLambda;
  std::vector<double> hPrevalence;
  std::vector<double> vPrevalence;

  bool useVectorPrevalence;
};

// void print_state (size_t iter, gsl_multiroot_fdfsolver * s, size_t n)
void print_state (size_t iter, gsl_multiroot_fsolver * s, size_t n)
{
  printf ("iter = %3u x =", iter);
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3f", gsl_vector_get(s->x, i));
  }
  printf("  f(x) =");
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3e", gsl_vector_get(s->f, i));
  }
  printf("\n");
}

int betafunc_f(const gsl_vector * x, void * p, gsl_vector * f)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  if (!params->useVectorPrevalence) {
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      double sum = 0;
      for (size_t j = 0; j < params->hosts.size(); ++j) {
        sum += gsl_vector_get(x, j) * params->hosts[j].theta *
          params->hPrevalence[j];
      }
      params->vPrevalence[i] = params->vectors[i].bitingRate * sum /
        (params->vectors[i].mu + params->vectors[i].bitingRate * sum);
    }
  }

  double vectorSum = 0;
  for (size_t i = 0; i < params->vectors.size(); ++i) {
    vectorSum += params->vectors[i].bitingRate * params->vectors[i].density *
      params->vPrevalence[i]; 
  }

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    double beta = gsl_vector_get(x, i);
    double y = params->hLambda[i] - beta * params->hosts[i].theta *
      params->params.areaConvert * vectorSum / params->hosts[i].abundance;
    gsl_vector_set(f, i, y);
  }

  return GSL_SUCCESS;
}

int betafunc_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  if (params->useVectorPrevalence) {
    double vectorSum = 0;
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      vectorSum += params->vectors[i].density * params->vPrevalence[i];
    }
    
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      for (size_t j = 0; j < params->hosts.size(); ++j) {
        if (i == j) {
          double y = - params->params.areaConvert * vectorSum /
            params->hosts[i].abundance;
          gsl_matrix_set(J, i, j, y);
        } else {
          gsl_matrix_set(J, i, j, 0);
        }
      }
    }
  } else {
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      for (size_t j = 0; j < params->hosts.size(); ++j) {

        double sum = 0;

        for (size_t l = 0; l < params->hosts.size(); ++l) {
          sum += gsl_vector_get(x, l) *
            params->hosts[l].theta * params->hPrevalence[l];
        }
        
        double enumerator;
        double denominator;
        
        double y = 0;
        
        for (size_t k = 0; k < params->vectors.size(); ++k) {
          
          enumerator = - params->params.areaConvert * params->vectors[k].bitingRate *
            params->hosts[i].abundance * gsl_vector_get(x, i) * params->hosts[i].theta;
          if (i == j) {
            enumerator -= params->hosts[i].abundance * params->params.areaConvert *
              (params->vectors[k].mu + sum) * sum;
          }
          denominator = 1 / ((params->vectors[k].mu + sum) * (params->vectors[k].mu + sum));
          
          y += enumerator /denominator;
        }
        gsl_matrix_set(J, i, j, y);
      }
    }
  }

  return GSL_SUCCESS;
}

int betafunc_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J)
{
  betafunc_f(x, p, f);
  betafunc_df(x, p, J);

  return GSL_SUCCESS;
}

void betaffoiv(void *p, std::vector<double> &beta)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    params->hPrevalence[i] = params->hosts[i].M /
      static_cast<double>(params->hosts[i].N);
    params->hLambda[i] = params->hPrevalence[i] /
      (1-params->hPrevalence[i]) *
      (params->hosts[i].gamma + params->hosts[i].mu);
  }
  
  for (size_t i = 0; i < params->vectors.size(); ++i) {
    params->vPrevalence[i] = params->vectors[i].M /
      static_cast<double>(params->vectors[i].N);
  }
  // std::cout << "rlambda:";
  // for (size_t i = 0; i < params->hosts.size(); ++i) {
  //   std::cout << " " << params->hLambda[i];
  // }
  // std::cout << std::endl;

  // std::cout << "vdensity:";
  // for (size_t i = 0; i < params->vectors.size(); ++i) {
  //   std::cout << " " << params->vectors[i].density;
  // }
  // std::cout << std::endl;

  // std::cout << "rabundance:";
  // for (size_t i = 0; i < params->hosts.size(); ++i) {
  //   std::cout << " " << params->hosts[i].abundance;
  // }
  // std::cout << std::endl;
  
  // std::cout << "theta:";
  // for (size_t i = 0; i < params->hosts.size(); ++i) {
  //   std::cout << " " << params->hosts[i].theta;
  // }
  // std::cout << std::endl;
  
  // std::cout << "biting_rate:";
  // for (size_t i = 0; i < params->vectors.size(); ++i) {
  //   std::cout << " " << params->vectors[i].bitingRate;
  // }
  // std::cout << std::endl;

  // std::cout << "area_convert:";
  // std::cout << " " << params->params.areaConvert;
  // std::cout << std::endl;

  // std::cout << "vprev:";
  // if (params->useVectorPrevalence) {
  //   for (size_t i = 0; i < params->vectors.size(); ++i) {
  //     std::cout << " " << params->vPrevalence[i];
  //   }
  // }
  // std::cout << std::endl;
  
  // std::cout << "rprev:";
  // for (size_t i = 0; i < params->hosts.size(); ++i) {
  //   std::cout << " " << params->hPrevalence[i];
  // }
  // std::cout << std::endl;
  
  // std::cout << "vmu:";
  // if (!params->useVectorPrevalence) {
  //   for (size_t i = 0; i < params->vectors.size(); ++i) {
  //     std::cout << " " << params->vectors[i].mu;
  //   }
  // }
  // std::cout << std::endl;
  
  // const gsl_multiroot_fdfsolver_type * T =
  //   gsl_multiroot_fdfsolver_hybridsj;
  // gsl_multiroot_fdfsolver * s =
  //   gsl_multiroot_fdfsolver_alloc (T, params->hosts.size());
  // gsl_multiroot_function_fdf f = {&betafunc_f, &betafunc_df, &betafunc_fdf,
  //                                 params->hosts.size(), p};

  const gsl_multiroot_fsolver_type * T =
    gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver * s =
    gsl_multiroot_fsolver_alloc (T, params->hosts.size());
  gsl_multiroot_function f = {&betafunc_f, params->hosts.size(), p};

  gsl_vector* x_init = gsl_vector_alloc(params->hosts.size());
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    gsl_vector_set(x_init, i, 1);
  }
  
  // gsl_multiroot_fdfsolver_set(s, &f, x_init);
  gsl_multiroot_fsolver_set(s, &f, x_init);

  size_t iter = 0;
  // print_state (iter, s, params->hosts.size());

  int status;
  do {
    iter++;
    // status = gsl_multiroot_fdfsolver_iterate (s);
    status = gsl_multiroot_fsolver_iterate (s);
    // print_state (iter, s, params->hosts.size());

    if (status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);

  // printf ("status = %s\n", gsl_strerror (status));

  beta.resize(params->hosts.size());

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    beta[i] = gsl_vector_get(s->x, i);
  }
  
  // gsl_multiroot_fdfsolver_free (s);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x_init);
}
#endif
