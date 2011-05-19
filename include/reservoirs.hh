#ifndef RESERVOIRS_HH
#define RESERVOIRS_HH

#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

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
      } else if (header[i] == "mortality_low") {
        std::istringstream s(data[i]);
        s >> mu_limits.first;
      } else if (header[i] == "mortality_high") {
        std::istringstream s(data[i]);
        s >> mu_limits.second;
      } else if (header[i] == "rec_rate") {
        std::istringstream s(data[i]);
        s >> gamma;
      } else if (header[i] == "rec_rate_low") {
        std::istringstream s(data[i]);
        s >> gamma_limits.first;
      } else if (header[i] == "rec_rate_high") {
        std::istringstream s(data[i]);
        s >> gamma_limits.second;
      } else if (header[i] == "abundance") {
        std::istringstream s(data[i]);
        s >> abundance;
      } else if (header[i] == "abundance_low") {
        std::istringstream s(data[i]);
        s >> abundance_limits.first;
      } else if (header[i] == "abundance_high") {
        std::istringstream s(data[i]);
        s >> abundance_limits.second;
      } else if (header[i] == "theta") {
        std::istringstream s(data[i]);
        s >> theta;
      }
    }
  }

  std::string name;
  size_t M, N;
  double mu, gamma, theta, abundance;
  std::pair<double, double> mu_limits, gamma_limits, abundance_limits;
  std::vector<double> habitat;
  std::vector<std::pair<double, double> > habitat_limits;
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

  vector(const vector& v) :
    name(v.name),
    M(v.M),
    N(v.N),
    mu(v.mu),
    gamma(v.gamma),
    density(v.density),
    bitingRate(v.bitingRate)
  {}

  std::string name;
  size_t M, N;
  double mu, gamma, density, bitingRate;
};

struct group
{
  group() :
    theta(.0)
  {}
  
  group(size_t initMember, double theta) :
    theta(theta),
    members(std::vector<size_t>(1, initMember))
  {}
  
  double theta;
  std::vector<size_t> members;
};

// parameters of the differnetial equation system
struct betafunc_params
{

  betafunc_params(std::vector<host> const &hosts,
                  std::vector<vector> const &vectors,
                  std::vector<group> const &groups,
                  double xi, double eta) :
    hosts(hosts),
    vectors(vectors),
    groups(groups),
    xi(xi),
    eta(std::make_pair(eta, eta)),
    hPrevalence(std::vector<double>(hosts.size()))
  {
    for (size_t i = 0; i < hosts.size(); ++i) {
      hPrevalence[i] = hosts[i].M /
        static_cast<double>(hosts[i].N);
    }

    vPrevalence = vectors[0].M /
        static_cast<double>(vectors[0].N);
  }

  betafunc_params(std::vector<host> const &hosts,
                  std::vector<vector> const &vectors,
                  std::vector<group> const &groups,
                  double xi, std::pair<double> eta) :
    hosts(hosts),
    vectors(vectors),
    groups(groups),
    xi(xi),
    eta(eta),
    hPrevalence(std::vector<double>(hosts.size()))
  {
    for (size_t i = 0; i < hosts.size(); ++i) {
      hPrevalence[i] = hosts[i].M /
        static_cast<double>(hosts[i].N);
    }

    vPrevalence = vectors[0].M /
        static_cast<double>(vectors[0].N);
  }

  std::vector<host> const &hosts;
  std::vector<vector> const &vectors;
  std::vector<group> const &groups;

  double xi;

  std::pair<double, double> eta;
  
  eta::vector<double> hPrevalence;
  double vPrevalence;
};

// print current state of the fdfsolver
void print_state (size_t iter, gsl_multiroot_fdfsolver * s, size_t n)
{
  printf ("iter = %3u x =", static_cast<unsigned int>(iter));
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3f", gsl_vector_get(s->x, i));
  }
  printf("  f(x) =");
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3e", gsl_vector_get(s->f, i));
  }
  printf("\n");
}

// print current state of the fsolver
void print_state (size_t iter, gsl_multiroot_fsolver * s, size_t n)
{
  printf ("iter = %3u x =", static_cast<unsigned int>(iter));
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3f", gsl_vector_get(s->x, i));
  }
  printf("  f(x) =");
  for (size_t i = 0; i < n; ++i) {
    printf(" %.3e", gsl_vector_get(s->f, i));
  }
  printf("\n");
}

// function to find root of (host-assortative)
int betafunc_f(const gsl_vector * x, void * p, gsl_vector * f)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  double pv[params->groups.size()];
  double beta[params->hosts.size()];
  double alpha = gsl_vector_get(x, params->hosts.size() +
                                params->groups.size());
  
  double weightedVectorPrevSum = 0;
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    beta[i] = gsl_vector_get(x, i);
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    pv[j] = gsl_vector_get(x, j+params->hosts.size());
    weightedVectorPrevSum += pv[j] * params->groups[j].theta;
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    double yv = (pv[j] * params->vectors[j].mu + params->xi *
                 (pv[j] - params->vPrevalence)) /
      (1 - pv[j]);
    for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
      size_t i = params->groups[j].members[k];
      double yh = params->hPrevalence[i] / (1 - params->hPrevalence[i]) *
        (params->hosts[i].mu + params->hosts[i].gamma) -
        beta[i] * params->hosts[i].theta / params->hosts[i].abundance *
        pv[j];
      yv -= alpha * beta[i] * params->hosts[i].theta * params->hPrevalence[i] /
        params->groups[j].theta;
      gsl_vector_set(f, i, yh);
    }
    gsl_vector_set(f, j + params->hosts.size(), yv);
  }
  gsl_vector_set(f, params->hosts.size() + params->groups.size(),
                 params->vPrevalence - weightedVectorPrevSum);
    
  return GSL_SUCCESS;
}

// function to find root of (pachy area version)
int betafunc_area_f(const gsl_vector * x, void * p, gsl_vector * f)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  double pv[params->groups.size()];
  double beta[params->hosts.size()];
  double alpha = gsl_vector_get(x, params->hosts.size() +
                                params->groups.size());
  
  double weightedVectorPrevSum = 0;
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    beta[i] = gsl_vector_get(x, i);
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    pv[j] = gsl_vector_get(x, j+params->hosts.size());
    weightedVectorPrevSum += pv[j] * params->groups[j].theta;
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    double yv = (pv[j] * params->vectors[j].mu + params->xi *
                 (pv[j] - params->vPrevalence)) /
      (1 - pv[j]);
    for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
      size_t i = params->groups[j].members[k];
      double yh = params->hPrevalence[i] / (1 - params->hPrevalence[i]) *
        (params->hosts[i].mu + params->hosts[i].gamma) -
        beta[i] * params->hosts[i].theta / params->hosts[i].abundance *
        pv[j];
      yv -= alpha * beta[i] * params->hosts[i].theta * params->hPrevalence[i] /
        params->groups[j].theta;
      gsl_vector_set(f, i, yh);
    }
    gsl_vector_set(f, j + params->hosts.size(), yv);
  }
  gsl_vector_set(f, params->hosts.size() + params->groups.size(),
                 params->vPrevalence - weightedVectorPrevSum);
    
  return GSL_SUCCESS;
}

// function to find root of (derivative)
int betafunc_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  double pv[params->groups.size()];
  double beta[params->hosts.size()];
  double alpha = gsl_vector_get(x, params->hosts.size() +
                                params->groups.size());
  
  for (size_t i = 0; i < params->hosts.size() + params->groups.size() + 1;
       ++i) {
    for (size_t j = 0; j < params->hosts.size() + params->groups.size() + 1;
         ++j) {
      gsl_matrix_set(J, i, j, 0);
    }
  }
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    beta[i] = gsl_vector_get(x, i);
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    pv[j] = gsl_vector_get(x, j+params->hosts.size());
  }

  for (size_t j = 0; j < params->groups.size(); ++j) {

    double dgjda = .0;

    for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
      size_t i = params->groups[j].members[k];

      double dfidbi = params->hosts[i].theta / params->hosts[i].abundance *
        pv[j];
      gsl_matrix_set(J, i, i, dfidbi);
      
      double dfidpj = beta[i] * params->hosts[i].theta /
        params->hosts[i].abundance;
      gsl_matrix_set(J, i, j + params->hosts.size(), dfidpj);

      double dgjdbi = -alpha * params->hosts[i].theta *
        params->hPrevalence[i] / params->groups[j].theta;
      gsl_matrix_set(J, j + params->hosts.size(), i, dgjdbi);

      double dgjdpj = -(params->vectors[0].mu + params->xi) /
        pow(1-pv[j], 2);
      gsl_matrix_set(J, j + params->hosts.size(), j + params->hosts.size(),
                     dgjdpj);
      
      dgjda -= beta[i] * params->hosts[i].theta * params->hPrevalence[i];
    }
    dgjda /= params->groups[j].theta;
    gsl_matrix_set(J, j + params->hosts.size(), params->hosts.size() +
                   params->groups.size(), dgjda);

    double dhdpj = params->groups[j].theta;
    gsl_matrix_set(J, params->hosts.size() + params->groups.size(),
                   j + params->hosts.size(), dhdpj);
  }      

  // for (size_t i = 0; i < params->hosts.size() + params->groups.size() + 1; ++i) {
  //   for (size_t j = 0; j < params->hosts.size() + params->groups.size() + 1; ++j) {
  //     std::cout << "J(" << i << "," << j << ") = " << gsl_matrix_get(J, i, j)
  //               << std::endl;
  //   }
  // }
  
  return GSL_SUCCESS;
}

// function to find root of (2nd derivative)
int betafunc_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J)
{
  betafunc_f(x, p, f);
  betafunc_df(x, p, J);

  return GSL_SUCCESS;
}

// find beta (and alpha and p^v_i( from forces of infection
int betaffoiv(void *p, std::vector<double> &vars,
              bool jac = false, unsigned int verbose = 0)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  bool area = false;
  if (params->xi == 0 && (params->eta.first == 0 || params->eta.second == 0)) {
    area = true;
  }

  if (verbose) {
    std::cout << "vdensity:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].density;
    }
    std::cout << std::endl;

    std::cout << "rabundance:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i].abundance;
    }
    std::cout << std::endl;
  
    std::cout << "theta:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i].theta;
    }
    std::cout << std::endl;
  
    std::cout << "biting_rate:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].bitingRate;
    }
    std::cout << std::endl;

    std::cout << "rprev:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hPrevalence[i];
    }
    std::cout << std::endl;
  
    std::cout << "vmu:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].mu;
    }
    std::cout << std::endl;

    std::cout << "xi: " << params->xi << std::endl;
    std::cout << "eta: " << params->eta.first << "," << params->eta.second << std::endl;
  }

  size_t nvars;
  
  if (area) {
    nvars = params->hosts.size() * params->hosts[0].habitat.size() + 2;
  } else {
    nvars = params->hosts.size() + params->groups.size() + 1;
  }

  const gsl_multiroot_fdfsolver_type * Tdf =
    gsl_multiroot_fdfsolver_hybridsj;
  gsl_multiroot_fdfsolver * sdf =
    gsl_multiroot_fdfsolver_alloc (Tdf, nvars);
  gsl_multiroot_function_fdf fdf = {&betafunc_f, &betafunc_df, &betafunc_fdf,
                                    nvars, p};

  const gsl_multiroot_fsolver_type * T =
    gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver * s =
     gsl_multiroot_fsolver_alloc (T, nvars);
  gsl_multiroot_function f;
  if (area) {
    f = {&betafunc_area_f, nvars, p};
  } else {
    f = {&betafunc_f, nvars, p};
  }

  gsl_vector* x_init = gsl_vector_alloc(nvars);
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    gsl_vector_set(x_init, i, 1);
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    gsl_vector_set(x_init, j + params->hosts.size(), params->vPrevalence);
  }
  gsl_vector_set(x_init, params->hosts.size() + params->groups.size(), 1);

  gsl_multiroot_fdfsolver_set(sdf, &fdf, x_init);
  gsl_multiroot_fsolver_set(s, &f, x_init);

  size_t iter = 0;
  if (verbose > 0) {
    if (jac) {
      print_state (iter, sdf, nvars);
    } else {
      print_state (iter, s, nvars);
    }
  }

  int status;
  do {
    iter++;
    if (jac) {
      status = gsl_multiroot_fdfsolver_iterate (sdf);
    } else {
      status = gsl_multiroot_fsolver_iterate (s);
    }

    if (verbose > 0) {
      if (jac) {
        print_state (iter, sdf, nvars);
      } else {
        print_state (iter, s, nvars);
      }
    }

    if (status)
      break;

    if (jac) {
      status = gsl_multiroot_test_residual (sdf->f, 1e-7);
    } else {
      status = gsl_multiroot_test_residual (s->f, 1e-7);
    }
  } while (status == GSL_CONTINUE); // && iter < 10000);

  if (verbose > 0) {
    printf ("status = %s\n", gsl_strerror (status));
  }

  vars.resize(nvars);
  gsl_vector* sol;

  if (jac) {
    sol = sdf->x;
  } else {
    sol = s->x;
  }

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    vars[i] = gsl_vector_get(sol, i);
  }
  for (size_t j = 0; j < params->groups.size(); ++j) {
    vars[j+params->hosts.size()] =
      gsl_vector_get(sol, j+params->hosts.size());
  }
  for (size_t i = 0; i < 1; ++i) {
    vars[i+params->hosts.size() + params->groups.size()] =
      gsl_vector_get(sol, i+params->hosts.size()+params->groups.size());
  }

  gsl_multiroot_fdfsolver_free (sdf);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x_init);

  return status;
 }

 #endif
