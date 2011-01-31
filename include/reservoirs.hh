#ifndef RESERVOIRS_HH
#define RESERVOIRS_HH

#include <vector>
#include <string>
#include <sstream>

struct host {

  host(std::vector<std::string>& data, std::vector<std::string>& header)
  {
    for (unsigned int i = 0; i < header.size(); ++i) {
      switch (header[i]) {
       case 'name':
        name = header[i];
        break;
       case 'N':
          std::istringstream s(header[i]);
          s >> N;
        break;
       case 'pos_tbng':
          std::istringstream s(header[i]);
          s >> posNG;
        break;
       case 'pos_tbg':
        std::istringstream s(header[i]);
        s >> posG;
        break;
       case 'mortality':
        std::istringstream s(header[i]);
        s >> mu;
        break;
       case 'rec_rate':
        std::istringstream s(header[i]);
        s >> gamma;
        break;
       case 'abundance':
        std::istringstream s(header[i]);
        s >> abundance;
        break;
       case 'theta':
        std::istringstream s(header[i]);
        s >> theta;
        break;
      }
    }
  }

  std::string name;
  unsigned int M, N;
  double mu, gamma, lambda, theta, abundance;
};

struct vector {

  vector(std::vector<std::string>& data, std::vector<std::string>& header)
  {
    for (unsigned int i = 0; i < header.size(); ++i) {
      switch (header[i]) {
       case 'name':
        name = header[i];
        break;
       case 'N':
          std::istringstream s(header[i]);
          s >> N;
        break;
       case 'pos_tbng':
          std::istringstream s(header[i]);
          s >> posNG;
        break;
       case 'pos_tbg':
        std::istringstream s(header[i]);
        s >> posG;
        break;
       case 'mortality':
        std::istringstream s(header[i]);
        s >> mu;
        break;
       case 'rec_rate':
        std::istringstream s(header[i]);
        s >> gamma;
        break;
       case 'density':
        std::istringstream s(header[i]);
        s >> density;
        break;
       case 'biting_rate':
        std::istringstream s(header[i]);
        s >> bitingRate;
        break;
      }
    }
  }

  unsigned int M, N;
  double prevalence, mu, gamma, lambda, density, bitingRate;
};

struct param {
  param(): areaConvert(1.) {}
  param(std::vector<std::string>& data, std::vector<std::string>& header)
  {
    for (unsigned int i = 0; i < header.size(); ++i) {
      switch (header[i]) {
       case 'area_convert':
          std::istringstream s(header[i]);
          s >> areaConvert;
          break;
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
                  bool uvp = false) :
    nSpecies(hosts.size()),
    nVectorSpecies(vectors.size()),
    area_convert(params.area_convert),
    useVectorPrevalence(uvp)
  {
    hLambda = new double [nSpecies];
    for (unsigned int i = 0; i < hosts.size(); ++i) {
      hLambda[i] = hosts[i].lambda;
    }

    hPrevalence = new double [nSpecies];
    for (unsigned int i = 0; i < hosts.size(); ++i) {
      hPrevalence[i] = hosts[i].lambda;
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
  
  unsigned int nSpecies;
  unsigned int nVectorSpecies;

  double* hLambda;
  double* vDensity;
  double* hAbundance;
  double* theta;
  double* hPrevalence
  double* bitingRate;
  double area_convert;

  bool useVectorPrevalence;
  double* vectorPrevalence;
  double* vMu;
};

int betafunc(gsl_vector * x, void * p, gsl_vector * f)
{
  
}

void betaffoiv(void *betafunc_params)
{
  const gsl_multiroot_fdfsolver_type * T
    = gsl_multiroot_fdfsolver_hybridsj;
  gsl_multiroot_fdfsolver * s =
    gsl_multiroot_fdfsolver_alloc (T, hosts.size());
  gsl_multiroot_fdfsolver_set
  
}
#endif
