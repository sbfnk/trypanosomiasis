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
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (header[i] == "name") {
        name = header[i];
      } else if (header[i] == "N") {
        std::istringstream s(header[i]);
        s >> N;
      } else if (header[i] == "pos_tbng") {
        if (tbng) {
          unsigned int posNG;
          std::istringstream s(header[i]);
          s >> posNG;
          M += posNG;
        }
      } else if (header[i] == "pos_tbg") {
        if (tbg) {
          unsigned int posG;
          std::istringstream s(header[i]);
          s >> posG;
          M += posG;
        }
      } else if (header[i] == "mortality") {
        std::istringstream s(header[i]);
        s >> mu;
      } else if (header[i] == "rec_rate") {
        std::istringstream s(header[i]);
        s >> gamma;
      } else if (header[i] == "abundance") {
        std::istringstream s(header[i]);
        s >> abundance;
      } else if (header[i] == "theta") {
        std::istringstream s(header[i]);
        s >> theta;
      }
    }
  }

  std::string name;
  unsigned int M, N;
  double mu, gamma, theta, abundance;
};

struct vector {

  vector(std::vector<std::string>& data, std::vector<std::string>& header,
         bool tbg = true, bool tbng = false) :
    M(0)
  {
    for (unsigned int i = 0; i < header.size(); ++i) {
       if (header[i] == "name") {
         name = header[i];
       } else if (header[i] == "N") {
         std::istringstream s(header[i]);
         s >> N;
       } else if (header[i] == "pos_tbng") {
         if (tbng) {
           unsigned int posNG;
           std::istringstream s(header[i]);
           s >> posNG;
           M += posNG;
         }
       } else if (header[i] == "pos_tbg") {
         if (tbg) {
           unsigned int posG;
           std::istringstream s(header[i]);
           s >> posG;
           M += posG;
         }
       } else if (header[i] == "mortality") {
         std::istringstream s(header[i]);
         s >> mu;
       } else if (header[i] == "rec_rate") {
         std::istringstream s(header[i]);
         s >> gamma;
       } else if (header[i] == "density") {
         std::istringstream s(header[i]);
         s >> density;
       } else if (header[i] == "biting_rate") {
        std::istringstream s(header[i]);
        s >> bitingRate;
       }
    }
  }

  std::string name;
  unsigned int M, N;
  double prevalence, mu, gamma, density, bitingRate;
};

struct param {
  param(): areaConvert(1.) {}
  param(std::vector<std::string>& data, std::vector<std::string>& header)
  {
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (header[i] == "area_convert") {
        std::istringstream s(header[i]);
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
    for (unsigned int i = 0; i < hosts.size(); ++i) {
      hPrevalence[i] = hosts[i].M / static_cast<double>(hosts[i].N);
    }
    hLambda = new double [nSpecies];
    for (unsigned int i = 0; i < hosts.size(); ++i) {
      hLambda[i] = hPrevalence[i] / (1-hPrevalence[i]) *
        (hosts[i].gamma + hosts[i].mu);
    }

    hAbundance = new double [nSpecies];
    theta = new double [nSpecies];

    vDensity = new double [nVectorSpecies];
    bitingRate = new double [nVectorSpecies];

    if (useVectorPrevalence) {
      vectorPrevalence = new double [nVectorSpecies];
    } else {
      vMu = new double [nVectorSpecies];
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
      delete [] vectorPrevalence;
    } else {
      delete [] vMu;
    }
  }
  
  unsigned int nSpecies;
  unsigned int nVectorSpecies;

  double* hLambda;
  double* vDensity;
  double* hAbundance;
  double* theta;
  double* hPrevalence;
  double* bitingRate;
  double area_convert;

  bool useVectorPrevalence;
  double* vectorPrevalence;
  double* vMu;
};

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

int betafunc(const gsl_vector * x, void * p, gsl_vector * f)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  double* vectorPrevalence;
  if (params->useVectorPrevalence) {
    vectorPrevalence = params->vectorPrevalence;
  } else {
    vectorPrevalence = new double [params->nVectorSpecies];
    for (unsigned int i = 0; i < params->nVectorSpecies; ++i) {
      double sum = 0;
      for (unsigned int j = 0; j < params->nSpecies; ++j) {
        sum += gsl_vector_get(x, j) * params->theta[j] *
          params->hPrevalence[j];
      }
      vectorPrevalence[i] = params->bitingRate[i] * sum /
        (params->vMu[i] + params->bitingRate[i] * sum);
    }
  }

  double vectorSum = 0;
  for (unsigned int i = 0; i < params->nVectorSpecies; ++i) {
    vectorSum += params->vDensity[i] * params->vectorPrevalence[i];
  }

  for (unsigned int i = 0; i < params->nSpecies; ++i) {
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
    for (unsigned int i = 0; i < params->nVectorSpecies; ++i) {
      vectorSum += params->vDensity[i] * params->vectorPrevalence[i];
    }
    
    for (unsigned int i = 0; i < params->nSpecies; ++i) {
      for (unsigned int j = 0; j < params->nSpecies; ++j) {
        if (i == j) {
          double y = - params->area_convert * vectorSum / params->hAbundance[i];
          gsl_matrix_set(J, i, j, y);
        } else {
          gsl_matrix_set(J, i, j, 0);
        }
      }
    }
  } else {
    for (unsigned int i = 0; i < params->nSpecies; ++i) {
      for (unsigned int j = 0; j < params->nSpecies; ++j) {

        double sum = 0;

        for (unsigned int j = 0; j < params->nSpecies; ++j) {
          sum += gsl_vector_get(x, j) *
            params->theta[j] * params->hPrevalence[j];
        }
        
        double enum;
        double denom;
        
        double y = 0;
        
        for (unsigned int k = 0; k < params->nVectorSpecies; ++k) {
          
        
          

          enum = - params->area_convert * params->biting_rate *
            params->hAbundance[i] * gsl_vector_get(x, i) * params->theta[i];
          if (i == j) {
            enum -= params->hAbundance[i] * params->area_convert *
              (params->vMu[k] + sum) * sum;
          }
          denom = 1 / ((params->vMu[k] + sum) * (params->vMu[k] + sum));
          
          y += enum /denom;
        }
        gsl_matrix_set(J, i, j, y);
      }
    }
  }

  return GSL_SUCCESS
}

int betafunc_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J)
{
  betafunc_f(x, p, f);
  betafunc_df(x, p, J);
}

void betaffoiv(void *betafunc_params, std::vector<double>& beta)
{
  const gsl_multiroot_fdfsolver_type * T =
    gsl_multiroot_fdfsolver_hybridsj;
  gsl_multiroot_fdfsolver * s =
    gsl_multiroot_fdfsolver_alloc (T, hosts.size());
  gsl_multiroot_function f = {&betafunc, &betafunc_df, &betafunc_fdf,
                              hosts.size(), &betafunc_params};

  gsl_vector* x = gsl_vector_alloc(hosts.size());
  for (size_t i = 0; i < hosts.size(); ++i) {
    gsl_vector_set(x, i, 1);
  }
  
  gsl_multiroot_fdfsolver_set(s, &f, x);

  print_state (iter, s);

  size_t iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate (s);
    print_state (iter, s);

    if (status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  beta.resize(hosts.size());

  for (size_t i = 0; i < hosts.size(); ++i) {
    beta[i] = gsl_vector_get(x, i);
  }
  
  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
}
#endif
