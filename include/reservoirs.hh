// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file reservoirs.hh
  \brief Header file for reservoir analysis based on prevalence data

  This includes the function to be solved to determine parameter values
*/
#ifndef RESERVOIRS_HH
#define RESERVOIRS_HH

#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

typedef typename std::pair<double, std::pair<double, double> > Parameter;

namespace po = boost::program_options;

/*! Base class Class for (host/vector/system) parameters */
class ParamContainer
{
public:

  /*! \brief Constructor.
    
  This is where classes derived from Model define the parameters and command
  line options.
  */
  ParamContainer(std::string optionName = "") :
    options(po::options_description("\nOptions ("+optionName+")"))
  {;}
  //! Destructor -- this is a virtual class
  virtual ~ParamContainer() = 0; 

  void ReadTable(std::vector<std::string> const &data,
                 std::vector<std::string> const &header);
  void ReadParams(po::variables_map const &vm);

  void addParam(std::string option,
                std::string description,
                Parameter* param);

  //! Accessor for options
  const po::options_description& getOptions() const
  { return options; }

protected:
  /*! \brief Parameters
    
    A map of command line options to the species parameters
  */
  std::map<std::string, Parameter*> params;
  //! Command line options
  po::options_description options;

  std::string name;
};
  
/*! Base class Class for (host/vector) parameters with habitats */
class HabitatContainer :
  public ParamContainer
{
public:

  HabitatContainer(std::string optionName = ""):
    ParamContainer(optionName) {;}
  //! Destructor -- this is a virtual class
  virtual ~HabitatContainer() = 0; 

  void Normalise();
  void ReadTable(std::vector<std::string> const &data,
                                   std::vector<std::string> const &header);
  void ReadParams(po::variables_map const &vm);
private:
  std::vector<Parameter> habitat;
};
  
/*! Class for host species and its data/parameters */
class Host :
  public HabitatContainer
{

  Host();
public:
  Parameter M, N, mu, gamma, f, n;
};

/*! Class for vector species and its data/parameters */
class Vector :
  public HabitatContainer
{

  Vector();
public:
  Parameter M, N, mu, tau, b, xi;
};

struct Group
{
  Group() :
    f(.0)
  {}
  
  Group(size_t initMember, double f) :
    f(f),
    members(std::vector<size_t>(1, initMember))
  {}
  
  double f;
  std::vector<size_t> members;
};

/*! Class for system-wide data/parameters */
class GlobalParams :
  public ParamContainer
{
  GlobalParams();
  
  void ReadParams(po::variables_map const &vm);

public:
  bool estimateXi; //!< estimate xi or alpha
};

/*! parameters of the differential equation system */
struct betafunc_params
{

  betafunc_params(std::vector<Host> const &hosts,
                  std::vector<Vector> const &vector,
                  std::vector<Group> const &groups,
                  GlobalParams const &global,
                  std::string overlapType);

  std::vector<Host> const &hosts;
  std::vector<Vector> const &vectors;
  std::vector<Group> const &groups;

  GlobalParams const &global;

  std::vector<double> hPrevalence; //!< auxiliary variable so the parasite
                                   //!prevalence in hosts does not have to be
                                   //!calculated every time
  std::vector<double> vPrevalence; //!< auxiliary variable so the parasite
                                   //!prevalence in vectors does not have to be
                                   //!calculated every time 
  std::vector<std::vector<double> > habitatOverlap;
};

/*! \brief Read data from table.
   
Reads parameters and measured data from table

*/
void ParamContainer::ReadTable(std::vector<std::string> const &data,
                               std::vector<std::string> const &header)
{
  for (size_t i = 0; i < header.size(); ++i) {
    for (std::map<std::string, Parameter*>::iterator it = params.begin();
         it != params.end(); it++) {
      if (header[i] == "name") {
        name = header[i];
      } else {
        if ((header[i] == it->first)) {
          std::istringstream s(data[i]);
          s >> it->second->first;
        } else if (header[i].substr(0, it->first.length() + 1) ==
                   (it->first + "_low")) {
          std::istringstream s(data[i]);
          s >> it->second->second.first;
        } else if (header[i].substr(0, it->first.length() + 1) ==
                   (it->first + "_high")) {
          std::istringstream s(data[i]);
          s >> it->second->second.second;
        }
      }
    }
  }
}
  
/*! \brief Read command line parameters.
   
This should be called after the table has been read and command line parameters
of the model haven been assigned. It initialises the model parameter variables
with the values found in the command line parameters.

\param[in] vm The map of command line parameters
*/
void ParamContainer::ReadParams(po::variables_map const &vm) {
  // if we have a name, we add a prefix to parameter names
  std::string paramPrefix = "";
  if (name != "") paramPrefix = name + "-";
  for (std::map<std::string, Parameter*>::iterator it = params.begin();
       it != params.end(); it++) {
    if (vm.count(paramPrefix + it->first)) {
      // command line parameter has been specified, assign to model variable
      it->second->first = vm[paramPrefix + it->first].as<double>();
    } else if (vm.count(paramPrefix + it->first + "_low")){
      it->second->second.first = vm[paramPrefix + it->first].as<double>();
    } else if (vm.count(paramPrefix + it->first + "_high")){
      it->second->second.second = vm[paramPrefix + it->first].as<double>();
    }
  }
}


void ParamContainer::addParam(std::string option,
                              std::string description,
                              Parameter* param)
{
  params.insert(std::make_pair(option, param));
  options.add_options()
    (option.c_str(), po::value<double>(), description.c_str());
}

void HabitatContainer::Normalise()
{
  double habitatSum = .0; //!< sum of all habitat contributions (for
                          //!normalisation) 
  for (size_t i = 0; i < habitat.size(); ++i) {
    habitatSum += habitat[i].first;
  }
  for (std::vector<Parameter>::iterator it = habitat.begin();
       it != habitat.end(); it++) {
    (*it).first /= habitatSum;
    (*it).second.first /= habitatSum;
    (*it).second.second /= habitatSum;
  }
}

void HabitatContainer::ReadTable(std::vector<std::string> const &data,
                                 std::vector<std::string> const &header)
{
  ParamContainer::ReadTable(data, header);

  for (size_t i = 0; i < header.size(); ++i) {
    if (header[i].substr(0,1) == "X") {
      std::istringstream s(data[i]);
      size_t index;
      size_t numPos = header[i].find_first_of("0123456789");
      std::istringstream numStr(header[i].substr(numPos));
      numStr >> index;
      while (index > habitat.size()) {
        habitat.push_back(std::make_pair(.0, std::make_pair(.0, .0)));
      }
      if (header[i].find("_low_") != std::string::npos) {
        s >> habitat[index-1].second.first;
      } else if (header[i].find("_high_") != std::string::npos) {
        s >> habitat[index-1].second.second;
      } else if (numPos == 1) {
        s >> habitat[index-1].first;
      }
    }
  }

  Normalise();
}

void HabitatContainer::ReadParams(po::variables_map const &vm) 
{
  for (size_t i = 0; i < habitat.size(); ++i) {
    std::stringstream ss;
    ss << i;
    if (vm.count(name + "-X" + ss.str())) {
      // command line parameter has been specified, assign to model variable
      habitat[i].first = vm[name + "-X" + ss.str()].as<double>();
    } else if (vm.count(name + "-X" + ss.str() + "_low")){
      habitat[i].second.first = vm[name + "-X" + ss.str() + "_low"].as<double>();
    } else if (vm.count(name + "-X" + ss.str() + "_high")){
      habitat[i].second.second = vm[name + "-X" + ss.str() + "_low"].as<double>();
    }
  }

  Normalise();
}

void GlobalParams::ReadParams(po::variables_map const &vm)
{
  ParamContainer::ReadParams(vm);
  if (vm.count("estimate-xi")) {
    estimateXi = true;
  }
}

Host::Host() :
  HabitatContainer("hosts") 
{
  addParam("N", "Population size", &N);
  addParam("M", "Number infected", &M);
  addParam("mu", "Mortality rate", &mu);
  addParam("gamma", "Recovery rate", &gamma);
  addParam("n", "Abundance", &n);
  addParam("f", "Biting preference", &f);
}
  
Vector::Vector() :
  HabitatContainer("vectors") 
{
  addParam("N", "Population size", &N);
  addParam("M", "Number infected", &M);
  addParam("mu", "Mortality rate", &mu);
  addParam("tau", "Biting rate", &tau);
  addParam("b", "Susceptibility", &b);
  addParam("xi", "Correlation", &xi);
}
  
GlobalParams::GlobalParams() :
  ParamContainer("global"),
  estimateXi(false)
{
  options.add_options()
    ("estimate-xi", "Estimate xi (with vector susceptibility set)");
}

betafunc_params::betafunc_params(std::vector<Host> const &hosts,
                                 std::vector<Vector> const &vector,
                                 std::vector<Group> const &groups,
                                 GlobalParams const &global,
                                 std::string overlapType) :
  hosts(hosts),
  vectors(vectors),
  groups(groups),
  global(global),
  hPrevalence(std::vector<double>(hosts.size())),
  vPrevalence(std::vector<double>(vectors.size()))
{
  for (size_t i = 0; i < hosts.size(); ++i) {
    hPrevalence[i] = hosts[i].M / hosts[i].N;
  }
    
  for (size_t v = 0; v < vectors.size(); ++v) {
    vPrevalence[v] = vectors[v].M / vectors[v].N;
  }
    
  //!< amount of overlap between host habitats
  habitatOverlap = std::vector<std::vector<double> >
    (groups.size(), std::vector<double>(groups.size(), .0));
    
  if (overlapType == "b" ||
      overlapType == "f") {
    std::vector<double> normaliseSum(hosts.size(), .0);
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        for (size_t k = 0; k < groups[j].members.size(); ++k) {
          size_t i = groups[j].members[k];
          for (size_t n = 0; n < groups[m].members.size(); ++n) {
            size_t l = groups[m].members[n];
            // std::cout << i << " " << l << " " << hosts[0].habitat.size()
            //           << std::endl;
            for (size_t o = 0; o < hosts[0].habitat.size(); ++o) {
              // std::cout << "A " << i << " " << l << " " << o
              //           << " " << hosts[i].habitat[o] << " "
              //           << hosts[l].habitat[o] << std::endl;
              if (hosts[i].habitat[o] > 0 &&
                  hosts[l].habitat[o] > 0) {
                if (vm["habitat"].as<std::string>() == "b") {
                  habitatOverlap[j][m] = 1;
                  habitatOverlap[m][j] = 1;
                } else {
                  double overlap = hosts[i].habitat[o] * hosts[l].habitat[o];
                  habitatOverlap[j][m] += overlap;
                  habitatOverlap[m][j] += overlap;
                }
              }
            }
          }
        }
        normaliseSum[j] += habitatOverlap[j][m] * groups[m].f;
        normaliseSum[m] += habitatOverlap[m][j] * groups[j].f;
        // std::cout << "Y " << j << " " << m << " " << habitatOverlap[j][m]
        //           << " " << groups[m].f << " " << groups[j].f
        //           << std::endl;
      }
    }
    
    // normalise
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        habitatOverlap[j][m] *= groups[j].f / normaliseSum[j];
        habitatOverlap[m][j] *= groups[m].f / normaliseSum[m];
        // std::cout << "X " << j << " " << m << " "
        //           << normaliseSum[j] << " " << habitatOverlap[j][m]
        //           << std::endl;
      }
    }
  } else if (overlapType == "n") {
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        habitatOverlap[j][m] = groups[j].f;
        habitatOverlap[m][j] = groups[m].f;
      }
    }
  } else {
    std::cerr << "WARNING: Ignoring unknown habitat option "
              << overlapType
              << std::endl;
  }
}

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

  double bhost[params->hosts.size()];
  double bvector[params->vectors.size()];
  double pv[params->vectors.size()][params->groups.size()];
  double xi[params->vectors.size()];
  
  double weightedVectorPrevSum[params->vectors.size()];
  std::vector<std::vector<double> > incomingVectorSum
    (params->vectors.size(), std::vector<double>(params->groups.size(), .0);
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    bhost[i] = gsl_vector_get(x, i);
  }
  if (estimateXi) {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      bvector[v] = params->vectors[v].b;
      xi[v] = gsl_vector_get(params->hosts.size());
    }
  } else {
    for (size_t v = 0; i < params->vectors.size(); ++v) {
      bvector[v] = gsl_vector_get(x, v + params->hosts.size());
    }
    xi[v] = params->vectors[v].xi;
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      pv[v][j] =
        gsl_vector_get(x, j + v * params->groups.size() +
                       params->hosts.size() + params->vectors.size());
    }
  }
  double yv[params->vectors.size()];
  double yh[params->hosts.size()];
  for (size_t j = 0; j < params->groups.size(); ++j) {
    double enumerator = .0;
    double denominator = .0;
    for (size_t l = 0; l < params->groups.size(); ++l) {
      for (size_t v = 0; v < params->vectors.size(); ++v) {
        enumerator += pv[v][l] * params->habitatOverlap[j][l] *
          params->groups[l].f;
        denominator += params->habitatOverlap[j][l] *
          params->groups[v][l].f;
        // std::cout << j << " " << l << " " << incomingVectorSum[j] << " "
        //           << pv[l] << " " << params->habitatOverlap[j][l]
        //           << " " << params->groups[l].f << std::endl;
      }
      incomingVectorSum[v][j] = enumerator / denominator;
      // std::cout << j << " " << incomingVectorSum[j] << std::endl;
    }
    yv[v] = (pv[v][j] * params->vectors[v].mu + params->xi *
             (pv[v][j] - incomingVectorSum[v][j])) /
      (1 - pv[v][j]);
  }
  for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
    size_t i = params->groups[j].members[k];
    yh[i] = params->hPrevalence[i] / (1 - params->hPrevalence[i]) *
      (params->hosts[i].mu + params->hosts[i].gamma);
  }
  
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        yh[i] -= bhost[i] * params->vectors[v].tau * params->hosts[i].f /
          params->hosts[i].n * pv[v][j];
        yv[v] -= bvector[v] * params->vectors[v].tau * params->hosts[i].f *
        params->hPrevalence[i] / params->groups[j].f;
        // std::cout << k << " " << i << " " << j << " " << yh << " " << yv << " "
        //           << params->groups[j].f << " " << alpha
        //           << " " << bhost[i] << " " << params->hosts[i].f
        //           << " " << params->hPrevalence[i] << " " << pv[j]
        //           << " " << params->vectors[j].mu << " "
        //           << params->xi << " " << incomingVectorSum[j] << std::endl;
      }
    }
  }
    
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        gsl_vector_set(f, i, yh);
      }
      gsl_vector_set(f, j + v * params->groups.size() +
                     params->hosts.size() + params->vectors.size(), yv[v]);
      weightedVectorPrevSum[v] += pv[v][j] * params->groups[j].f;
    }
    gsl_vector_set(f, params->hosts.size() +
                   params->groups.size() * params->vectors.size() +
                   params->vectors.size(),
                   params->vPrevalence[v] - weightedVectorPrevSum[v]);
  }

    
  return GSL_SUCCESS;
}

// // function to find root of (derivative)
// int betafunc_df(const gsl_vector * x, void * p, gsl_matrix * J)
// {
//   betafunc_params* params = ((struct betafunc_params*) p);

//   double pv[params->groups.size()];
//   double bhost[params->hosts.size()];
//   double alpha = gsl_vector_get(x, params->hosts.size() +
//                                 params->groups.size());
  
//   for (size_t i = 0; i < params->hosts.size() + params->groups.size() + 1;
//        ++i) {
//     for (size_t j = 0; j < params->hosts.size() + params->groups.size() + 1;
//          ++j) {
//       gsl_matrix_set(J, i, j, 0);
//     }
//   }
//   for (size_t i = 0; i < params->hosts.size(); ++i) {
//     bhost[i] = gsl_vector_get(x, i);
//   }
//   for (size_t j = 0; j < params->groups.size(); ++j) {
//     pv[j] = gsl_vector_get(x, j+params->hosts.size());
//   }

//   for (size_t j = 0; j < params->groups.size(); ++j) {

//     double dgjda = .0;

//     for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
//       size_t i = params->groups[j].members[k];

//       double dfidbi = params->hosts[i].f / params->hosts[i].n *
//         pv[j];
//       gsl_matrix_set(J, i, i, dfidbi);
      
//       double dfidpj = bhost[i] * params->hosts[i].f /
//         params->hosts[i].n;
//       gsl_matrix_set(J, i, j + params->hosts.size(), dfidpj);

//       double dgjdbi = -alpha * params->hosts[i].f *
//         params->hPrevalence[i] / params->groups[j].f;
//       gsl_matrix_set(J, j + params->hosts.size(), i, dgjdbi);

//       double dgjdpj = -(params->vectors[0].mu + params->xi) /
//         pow(1-pv[j], 2);
//       gsl_matrix_set(J, j + params->hosts.size(), j + params->hosts.size(),
//                      dgjdpj);
      
//       dgjda -= bhost[i] * params->hosts[i].f * params->hPrevalence[i];
//     }
//     dgjda /= params->groups[j].f;
//     gsl_matrix_set(J, j + params->hosts.size(), params->hosts.size() +
//                    params->groups.size(), dgjda);

//     double dhdpj = params->groups[j].f;
//     gsl_matrix_set(J, params->hosts.size() + params->groups.size(),
//                    j + params->hosts.size(), dhdpj);
//   }      

//   // for (size_t i = 0; i < params->hosts.size() + params->groups.size() + 1; ++i) {
//   //   for (size_t j = 0; j < params->hosts.size() + params->groups.size() + 1; ++j) {
//   //     std::cout << "J(" << i << "," << j << ") = " << gsl_matrix_get(J, i, j)
//   //               << std::endl;
//   //   }
//   // }
  
//   return GSL_SUCCESS;
// }

// // function to find root of and derivative
// int betafunc_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J)
// {
//   betafunc_f(x, p, f);
//   betafunc_df(x, p, J);

//   return GSL_SUCCESS;
// }

// find beta (and alpha and p^v_i) from forces of infection
int betaffoiv(void *p, std::vector<double> &vars,
              bool jac = false, unsigned int verbose = 0)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  if (verbose) {
    std::cout << "vdensity:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].density;
    }
    std::cout << std::endl;

    std::cout << "rabundance:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i].n;
    }
    std::cout << std::endl;
  
    std::cout << "biting:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i].f;
    }
    std::cout << std::endl;
  
    std::cout << "biting_rate:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].tau;
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

    std::cout << "xi:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i].xi;
    }
    std::cout << std::endl;
  }

  size_t nvars =
    params->hosts.size() + params->groups.size() * params->vectors.size() +
    params->vectors.size() + 1;

  // beta + p_V + alpha
  
  // const gsl_multiroot_fdfsolver_type * Tdf;
  // gsl_multiroot_fdfsolver * sdf;
  // gsl_multiroot_function_fdf fdf;
  // if (jac) {
  //   Tdf = gsl_multiroot_fdfsolver_hybridsj;
  //   sdf = gsl_multiroot_fdfsolver_alloc (Tdf, nvars);
  //   fdf.f = &betafunc_f;
  //   fdf.df = &betafunc_df;
  //   fdf.fdf = &betafunc_fdf;
  //   fdf.n = nvars;
  //   fdf.params = p;
  // }

  const gsl_multiroot_fsolver_type * T =
    gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver * s =
     gsl_multiroot_fsolver_alloc (T, nvars);
  gsl_multiroot_function f;

  f.f = &betafunc_f;

  f.n = nvars;
  f.params = p;

  gsl_vector* x_init = gsl_vector_alloc(nvars);
  for (size_t i = 0; i < params->hosts.size(); ++i) {
    gsl_vector_set(x_init, i, 1);
  }
  if (estimateXi) {
    gsl_vector_set(x_init, params->hosts.size(), .5);
  } else {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      gsl_vector_set(x_init, params->hosts.size() + v, 1);
    }
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      gsl_vector_set(x_init, j + v * params->groups.size() +
                     params->hosts.size() + params->vectors.size(),
                     params->vPrevalence);
    }
  }

  if (jac) {
    gsl_multiroot_fdfsolver_set(sdf, &fdf, x_init);
  }
  gsl_multiroot_fsolver_set(s, &f, x_init);

  size_t iter = 0;
  if (verbose > 0) {
    // if (jac) {
    //   print_state (iter, sdf, nvars);
    // } else {
      print_state (iter, s, nvars);
    // }
  }

  int status;
  do {
    iter++;
    // if (jac) {
    //   status = gsl_multiroot_fdfsolver_iterate (sdf);
    // } else {
      status = gsl_multiroot_fsolver_iterate (s);
    // }

    if (verbose > 0) {
      // if (jac) {
      //   print_state (iter, sdf, nvars);
      // } else {
        print_state (iter, s, nvars);
      // }
    }

    if (status)
      break;

    // if (jac) {
    //   status = gsl_multiroot_test_residual (sdf->f, 1e-7);
    // } else {
      status = gsl_multiroot_test_residual (s->f, 1e-7);
    // }
  } while (status == GSL_CONTINUE); // && iter < 10000);

  if (verbose > 0) {
    printf ("status = %s\n", gsl_strerror (status));
  }

  vars.resize(nvars);
  gsl_vector* sol;

  // if (jac) {
  //   sol = sdf->x;
  // } else {
    sol = s->x;
  // }

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    vars[i] = gsl_vector_get(sol, i);
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    vars[v+params->hosts.size()] =
      gsl_vector_get(sol, v+params->hosts.size());
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      vars[j + v * params->groups.size() + params->hosts.size() +
           params->vectors.size()] =
        gsl_vector_get(sol, j+params->hosts.size());
    }
  }
  
  // if (jac) {
  //   gsl_multiroot_fdfsolver_free (sdf);
  // }
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x_init);

  return status;
}

#endif
