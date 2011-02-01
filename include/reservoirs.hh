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
                  bool tbg = true,
                  bool tbng = false) :
    nSpecies(hosts.size()),
    nVectorSpecies(vectors.size()),
    area_convert(params.areaConvert),
    useVectorPrevalence(uvp)
  {
    hPrevalence = new double [nSpecies];
    for (size_t i = 0; i < hosts.size(); ++i) {
      hPrevalence[i] = hosts[i].M / static_cast<double>(hosts[i].N);
    }
    hLambda = new double [nSpecies];
    for (size_t i = 0; i < hosts.size(); ++i) {
      hLambda[i] = hPrevalence[i] / (1-hPrevalence[i]) *
        (hosts[i].gamma + hosts[i].mu);
    }

    hAbundance = new double [nSpecies];
    for (size_t i = 0; i < hosts.size(); ++i) {
      hAbundance[i] = hosts[i].abundance;
    }
    
    theta = new double [nSpecies];
    for (size_t i = 0; i < hosts.size(); ++i) {
      theta[i] = hosts[i].theta;
    }

    vDensity = new double [nVectorSpecies];
    for (size_t i = 0; i < vectors.size(); ++i) {
      vDensity[i] = vectors[i].density;
    }
    bitingRate = new double [nVectorSpecies];
    for (size_t i = 0; i < vectors.size(); ++i) {
      bitingRate[i] = vectors[i].bitingRate;
    }

    if (useVectorPrevalence) {
      vPrevalence = new double [nVectorSpecies];
      for (size_t i = 0; i < vectors.size(); ++i) {
        vPrevalence[i] = vectors[i].M / static_cast<double>(vectors[i].N);
      }
    } else {
      vMu = new double [nVectorSpecies];
      for (size_t i = 0; i < vectors.size(); ++i) {
        vMu[i] = vectors[i].mu;
      }
    }
  }

  ~betafunc_params() {
    delete [] hLambda;
    delete [] vDensity;
    delete [] hAbundance;
    delete [] theta;
    delete [] hPrevalence;
    delete [] bitingRate;
    if (useVectorPrevalence) {
      delete [] vPrevalence;
    } else {
      delete [] vMu;
    }
  }
  
  size_t nSpecies;
  size_t nVectorSpecies;

  double* hLambda;
  double* vDensity;
  double* hAbundance;
  double* theta;
  double* hPrevalence;
  double* bitingRate;
  double area_convert;

  bool useVectorPrevalence;
  double* vPrevalence;
  double* vMu;
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

  double* vPrevalence;
  if (params->useVectorPrevalence) {
    vPrevalence = params->vPrevalence;
  } else {
    vPrevalence = new double [params->nVectorSpecies];
    for (size_t i = 0; i < params->nVectorSpecies; ++i) {
      double sum = 0;
      for (size_t j = 0; j < params->nSpecies; ++j) {
        sum += gsl_vector_get(x, j) * params->theta[j] *
          params->hPrevalence[j];
      }
      vPrevalence[i] = params->bitingRate[i] * sum /
        (params->vMu[i] + params->bitingRate[i] * sum);
    }
  }

  double vectorSum = 0;
  for (size_t i = 0; i < params->nVectorSpecies; ++i) {
    vectorSum += params->vDensity[i] * vPrevalence[i];
  }

  for (size_t i = 0; i < params->nSpecies; ++i) {
    double beta = gsl_vector_get(x, i);
    double y = params->hLambda[i] - beta * params->area_convert * vectorSum /
      params->hAbundance[i];
    gsl_vector_set(f, i, y);
  }

  return GSL_SUCCESS;
}

int betafunc_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  if (params->useVectorPrevalence) {
    double vectorSum = 0;
    for (size_t i = 0; i < params->nVectorSpecies; ++i) {
      vectorSum += params->vDensity[i] * params->vPrevalence[i];
    }
    
    for (size_t i = 0; i < params->nSpecies; ++i) {
      for (size_t j = 0; j < params->nSpecies; ++j) {
        if (i == j) {
          double y = - params->area_convert * vectorSum / params->hAbundance[i];
          gsl_matrix_set(J, i, j, y);
        } else {
          gsl_matrix_set(J, i, j, 0);
        }
      }
    }
  } else {
    for (size_t i = 0; i < params->nSpecies; ++i) {
      for (size_t j = 0; j < params->nSpecies; ++j) {

        double sum = 0;

        for (size_t l = 0; l < params->nSpecies; ++l) {
          sum += gsl_vector_get(x, l) *
            params->theta[l] * params->hPrevalence[l];
        }
        
        double enumerator;
        double denominator;
        
        double y = 0;
        
        for (size_t k = 0; k < params->nVectorSpecies; ++k) {
          
          enumerator = - params->area_convert * params->bitingRate[k] *
            params->hAbundance[i] * gsl_vector_get(x, i) * params->theta[i];
          if (i == j) {
            enumerator -= params->hAbundance[i] * params->area_convert *
              (params->vMu[k] + sum) * sum;
          }
          denominator = 1 / ((params->vMu[k] + sum) * (params->vMu[k] + sum));
          
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

void betaffoiv(void *p, std::vector<double> &beta,
               size_t nSpecies)
{
  // betafunc_params* params = ((struct betafunc_params*) p);
  // std::cout << "rlambda:";
  // for (size_t i = 0; i < nSpecies; ++i) {
  //   std::cout << " " << params->hLambda[i];
  // }
  // std::cout << std::endl;

  // std::cout << "vdensity:";
  // for (size_t i = 0; i < params->nVectorSpecies; ++i) {
  //   std::cout << " " << params->vDensity[i];
  // }
  // std::cout << std::endl;

  // std::cout << "rabundance:";
  // for (size_t i = 0; i < nSpecies; ++i) {
  //   std::cout << " " << params->hAbundance[i];
  // }
  // std::cout << std::endl;
  
  // std::cout << "theta:";
  // for (size_t i = 0; i < nSpecies; ++i) {
  //   std::cout << " " << params->theta[i];
  // }
  // std::cout << std::endl;
  
  // std::cout << "biting_rate:";
  // for (size_t i = 0; i < params->nVectorSpecies; ++i) {
  //   std::cout << " " << params->bitingRate[i];
  // }
  // std::cout << std::endl;

  // std::cout << "area_convert:";
  // std::cout << " " << params->area_convert;
  // std::cout << std::endl;

  // std::cout << "vprev:";
  // if (params->useVectorPrevalence) {
  //   for (size_t i = 0; i < params->nVectorSpecies; ++i) {
  //     std::cout << " " << params->vPrevalence[i];
  //   }
  // }
  // std::cout << std::endl;
  
  // std::cout << "rprev:";
  // for (size_t i = 0; i < nSpecies; ++i) {
  //   std::cout << " " << params->hPrevalence[i];
  // }
  // std::cout << std::endl;
  
  // std::cout << "vmu:";
  // if (!params->useVectorPrevalence) {
  //   for (size_t i = 0; i < params->nVectorSpecies; ++i) {
  //     std::cout << " " << params->vMu[i];
  //   }
  // }
  // std::cout << std::endl;
  
  // const gsl_multiroot_fdfsolver_type * T =
  //   gsl_multiroot_fdfsolver_hybridsj;
  // gsl_multiroot_fdfsolver * s =
  //   gsl_multiroot_fdfsolver_alloc (T, nSpecies);
  // gsl_multiroot_function_fdf f = {&betafunc_f, &betafunc_df, &betafunc_fdf,
  //                                 nSpecies, p};

  const gsl_multiroot_fsolver_type * T =
    gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver * s =
    gsl_multiroot_fsolver_alloc (T, nSpecies);
  gsl_multiroot_function f = {&betafunc_f, nSpecies, p};

  gsl_vector* x_init = gsl_vector_alloc(nSpecies);
  for (size_t i = 0; i < nSpecies; ++i) {
    gsl_vector_set(x_init, i, 1);
  }
  
  // gsl_multiroot_fdfsolver_set(s, &f, x);
  gsl_multiroot_fsolver_set(s, &f, x_init);

  size_t iter = 0;
  // print_state (iter, s, nSpecies);

  int status;
  do {
    iter++;
    // status = gsl_multiroot_fdfsolver_iterate (s);
    status = gsl_multiroot_fsolver_iterate (s);
    // print_state (iter, s, nSpecies);

    if (status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);

  // printf ("status = %s\n", gsl_strerror (status));

  beta.resize(nSpecies);

  for (size_t i = 0; i < nSpecies; ++i) {
    beta[i] = gsl_vector_get(s->x, i);
  }
  
  // gsl_multiroot_fdfsolver_free (s);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x_init);
}
#endif
