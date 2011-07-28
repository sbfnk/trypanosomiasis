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
    options(new po::options_description("\nOptions ("+optionName+")"))
  {;}
  //! Destructor
  virtual ~ParamContainer()
  {;}

  void ReadTable(std::vector<std::string> const &data,
                 std::vector<std::string> const &header);
  void ReadParams(po::variables_map const &vm);

  //! Accessor for options
  const po::options_description* getOptions() const
  { return options; }

  //! Accessor for name
  const std::string getName() const
  { return name; }

protected:
  struct ParamInfo {
    ParamInfo(std::string o, std::string d, Parameter* p) :
      option(o), description(d), param(p) 
    {;}
    
    std::string option;
    std::string description;
    Parameter* param;
  };
    
  /*! \brief Parameters
    
    A map of command line options to the species parameters
  */
  std::vector<ParamInfo> params;
  //! Command line options
  po::options_description* options;

  std::string name;
};
  
/*! Base class Class for (host/vector) parameters with habitats */
class HabitatContainer :
  public ParamContainer
{
public:

  HabitatContainer(std::string optionName = "") :
    ParamContainer(optionName) {;}
  //! Destructor
  ~HabitatContainer() {;}

  void Normalise();
  void ReadTable(std::vector<std::string> const &data,
                 std::vector<std::string> const &header);
  void ReadParams(po::variables_map const &vm);
  std::vector<Parameter> habitat;

};
  
/*! Class for host species and its data/parameters */
class Host :
  public HabitatContainer
{
public:
  Host();
  Parameter M, N, mu, gamma, f, n;
};

/*! Class for vector species and its data/parameters */
class Vector :
  public HabitatContainer
{
public:
  Vector();
  
  void ReadParams(po::variables_map const &vm) {
    HabitatContainer::ReadParams(vm);
  }
    
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
public:
  GlobalParams();
  void ReadParams(po::variables_map const &vm);
  std::string habType;
  bool estimateXi; //!< estimate xi or alpha
};

/*! parameters of the differential equation system */
struct betafunc_params
{

  betafunc_params(std::vector<Host*> const &hosts,
                  std::vector<Vector*> const &vector,
                  std::vector<Group> const &groups,
                  GlobalParams const *global);

  std::vector<Host*> const &hosts;
  std::vector<Vector*> const &vectors;
  std::vector<Group> const &groups;

  GlobalParams const *global;

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
    for (std::vector<ParamInfo>::iterator it = params.begin();
         it != params.end(); it++) {
      if (header[i] == "name") {
        name = data[i];
      } else {
        if ((header[i] == it->option)) {
          std::istringstream s(data[i]);
          s >> it->param->first;
        } else if (header[i].substr(0, it->option.length() + 1) ==
                   (it->option + "_low")) {
          std::istringstream s(data[i]);

          s >> it->param->second.first;
        } else if (header[i].substr(0, it->option.length() + 1) ==
                   (it->option + "_high")) {
          std::istringstream s(data[i]);
          s >> it->param->second.second;
        }
      }
    }
  }

  if (name.length() > 0) {
    for (std::vector<ParamInfo>::iterator it = params.begin();
         it != params.end(); it++) {
      options->add_options()
        ((name+"-"+it->option).c_str(), po::value<double>(),
         it->description.c_str());
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
  for (std::vector<ParamInfo>::iterator it = params.begin();
       it != params.end(); it++) {
    if (vm.count(paramPrefix + it->option)) {
      // command line parameter has been specified, assign to model variable
      it->param->first = vm[paramPrefix + it->option].as<double>();
    } else if (vm.count(paramPrefix + it->option + "_low")){
      it->param->second.first = vm[paramPrefix + it->option].as<double>();
    } else if (vm.count(paramPrefix + it->option + "_high")){
      it->param->second.second = vm[paramPrefix + it->option].as<double>();
    }
  }
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
  ParamContainer::ReadParams(vm);
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
  habType = vm["habitat"].as<std::string>();
  if (!(habType == "b" || habType == "f")) {
    if (habType != "n") {
      std::cerr << "WARNING: Overlap type " << habType << "unrecognised, "
                << "not considering habitats" << std::endl;
    }
    habType = "n";
  }
}

Host::Host() :
  HabitatContainer("host")
{
  params.push_back(ParamInfo("N", "Population size", &N));
  params.push_back(ParamInfo("M", "Number infected", &M));
  params.push_back(ParamInfo("mu", "Mortality rate", &mu));
  params.push_back(ParamInfo("gamma", "Recovery rate", &gamma));
  params.push_back(ParamInfo("n", "Abundance", &n));
  params.push_back(ParamInfo("f", "Biting preference", &f));
}
  
Vector::Vector() :
  HabitatContainer("vector")
{
  params.push_back(ParamInfo("N", "Population size", &N));
  params.push_back(ParamInfo("M", "Number infected", &M));
  params.push_back(ParamInfo("mu", "Mortality rate", &mu));
  params.push_back(ParamInfo("tau", "Biting rate", &tau));
  params.push_back(ParamInfo("b", "Susceptibility", &b));
  params.push_back(ParamInfo("xi", "Correlation", &xi));
}
  
GlobalParams::GlobalParams() :
  ParamContainer("global"),
  estimateXi(false)
{
  options->add_options()
    ("estimate-xi", "Estimate xi (with vector susceptibility set)")
    ("habitat,a", po::value<std::string>()->default_value("n"),
     "type of habitat overlap (n=none,b=binary, f=fractional)")
    ;
}

betafunc_params::betafunc_params(std::vector<Host*> const &hosts,
                                 std::vector<Vector*> const &vectors,
                                 std::vector<Group> const &groups,
                                 GlobalParams const *global) :
  hosts(hosts),
  vectors(vectors),
  groups(groups),
  global(global),
  hPrevalence(std::vector<double>(hosts.size())),
  vPrevalence(std::vector<double>(vectors.size()))
{
  for (size_t i = 0; i < hosts.size(); ++i) {
    hPrevalence[i] = hosts[i]->M.first / hosts[i]->N.first;
  }
    
  for (size_t v = 0; v < vectors.size(); ++v) {
    vPrevalence[v] = vectors[v]->M.first / vectors[v]->N.first;
  }
    
  //!< amount of overlap between host habitats
  habitatOverlap = std::vector<std::vector<double> >
    (groups.size(), std::vector<double>(groups.size(), .0));
    
  if (global->habType == "b" ||
      global->habType == "f") {
    std::vector<double> normaliseSum(hosts.size(), .0);
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        for (size_t k = 0; k < groups[j].members.size(); ++k) {
          size_t i = groups[j].members[k];
          for (size_t n = 0; n < groups[m].members.size(); ++n) {
            size_t l = groups[m].members[n];
            // std::cout << i << " " << l << " " << hosts[0]->habitat.size()
            //           << std::endl;
            for (size_t o = 0; o < hosts[0]->habitat.size(); ++o) {
              // std::cout << "A " << i << " " << l << " " << o
              //           << " " << hosts[i]->habitat[o] << " "
              //           << hosts[l]->habitat[o] << std::endl;
              if (hosts[i]->habitat[o].first > 0 &&
                  hosts[l]->habitat[o].first > 0) {
                if (global->habType == "b") {
                  habitatOverlap[j][m] = 1;
                  habitatOverlap[m][j] = 1;
                } else {
                  double overlap =
                    hosts[i]->habitat[o].first *
                    hosts[l]->habitat[o].first;
                  habitatOverlap[j][m] += overlap;
                  habitatOverlap[m][j] += overlap;
                }
              }
            }
          }
        }
        // normaliseSum[j] += habitatOverlap[j][m] * groups[m].f;
        // normaliseSum[m] += habitatOverlap[m][j] * groups[j].f;
        // std::cout << "Y " << j << " " << m << " " << habitatOverlap[j][m]
        //           << " " << groups[m].f << " " << groups[j].f
        //           << std::endl;
      }
    }
    
    // normalise
    // for (size_t j = 0; j < groups.size(); ++j) {
    //   for (size_t m = j; m < groups.size(); ++m) {
    //     habitatOverlap[j][m] *= groups[j].f / normaliseSum[j];
    //     habitatOverlap[m][j] *= groups[m].f / normaliseSum[m];
    //     // std::cout << "X " << j << " " << m << " "
    //     //           << normaliseSum[j] << " " << habitatOverlap[j][m]
    //     //           << std::endl;
    //   }
    // }
  } else {
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        // habitatOverlap[j][m] = groups[j].f;
        // habitatOverlap[m][j] = groups[m].f;
        habitatOverlap[j][m] = 1;
        habitatOverlap[m][j] = 1;
      }
    }
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

  std::vector<double> bhost(params->hosts.size(), .0);
  std::vector<double> bvector(params->vectors.size(), .0);
  std::vector<std::vector<double> > pv
    (params->groups.size(), std::vector<double>(params->vectors.size(), .0));
  std::vector<double> xi(params->vectors.size(), .0);
  
  std::vector<double> weightedVectorPrevSum(params->vectors.size(), .0);
  std::vector<std::vector<double> > incomingVectorSum
    (params->vectors.size(), std::vector<double>(params->groups.size(), .0));

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    bhost[i] = gsl_vector_get(x, i);
  }
  if (params->global->estimateXi) {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      bvector[v] = params->vectors[v]->b.first;
      xi[v] = gsl_vector_get(x, v + params->hosts.size());
    }
  } else {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      bvector[v] = gsl_vector_get(x, v + params->hosts.size());
      xi[v] = params->vectors[v]->xi.first;
    }
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      pv[j][v] =
        gsl_vector_get(x, j + v * params->groups.size() +
                       params->hosts.size() + params->vectors.size());
    }
  }
  double yv[params->groups.size()][params->vectors.size()];
  double yh[params->hosts.size()];
  for (size_t j = 0; j < params->groups.size(); ++j) {
    double enumerator = .0;
    double denominator = .0;
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      // for (size_t l = 0; l < params->groups.size(); ++l) {
      //   enumerator += pv[l][v] * params->habitatOverlap[j][l] *
      //     params->groups[l].f;
      //   denominator += params->habitatOverlap[j][l] *
      //     params->groups[l].f;
      // }
      // incomingVectorSum[v][j] = enumerator / denominator;
      // std::cout << pv[j][v] << " " << incomingVectorSum[v][j] << std::endl;
      yv[j][v] = (pv[j][v] * params->vectors[v]->mu.first + xi[v] *
                  (pv[j][v] - params->vPrevalence[v]) /
        (1 - pv[j][v]);
    }
    for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
      size_t i = params->groups[j].members[k];
      yh[i] = params->hPrevalence[i] / (1 - params->hPrevalence[i]) *
        (params->hosts[i]->mu.first + params->hosts[i]->gamma.first);
    }
  }
  
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        yh[i] -= bhost[i] * params->vectors[v]->tau.first * params->hosts[i]->f.first /
          params->hosts[i]->n.first * pv[j][v];
        yv[j][v] -= bvector[v] * params->vectors[v]->tau.first *
          params->hosts[i]->f.first * params->hPrevalence[i] / params->groups[j].f;
        // std::cout << k << " " << i << " " << j << " " << yh << " " << yv << " "
        //           << params->groups[j].f << " " << alpha
        //           << " " << bhost[i] << " " << params->hosts[i]->f
        //           << " " << params->hPrevalence[i] << " " << pv[j]
        //           << " " << params->vectors[j]->mu << " "
        //           << params->xi << " " << incomingVectorSum[j] << std::endl;
      }
    }
  }
    
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        gsl_vector_set(f, i, yh[i]);
      }
      gsl_vector_set(f, j + v * params->groups.size() +
                     params->hosts.size(), yv[j][v]);
      weightedVectorPrevSum[v] += pv[j][v] * params->groups[j].f;
    }
    gsl_vector_set(f, params->hosts.size() +
                   params->groups.size() * params->vectors.size() + v,
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

//       double dfidbi = params->hosts[i]->f / params->hosts[i]->n *
//         pv[j];
//       gsl_matrix_set(J, i, i, dfidbi);
      
//       double dfidpj = bhost[i] * params->hosts[i]->f /
//         params->hosts[i]->n;
//       gsl_matrix_set(J, i, j + params->hosts.size(), dfidpj);

//       double dgjdbi = -alpha * params->hosts[i]->f *
//         params->hPrevalence[i] / params->groups[j].f;
//       gsl_matrix_set(J, j + params->hosts.size(), i, dgjdbi);

//       double dgjdpj = -(params->vectors[0]->mu + params->xi) /
//         pow(1-pv[j], 2);
//       gsl_matrix_set(J, j + params->hosts.size(), j + params->hosts.size(),
//                      dgjdpj);
      
//       dgjda -= bhost[i] * params->hosts[i]->f * params->hPrevalence[i];
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
    std::cout << "rabundance:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->n.first;
    }
    std::cout << std::endl;
  
    std::cout << "biting:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->f.first;
    }
    std::cout << std::endl;
  
    std::cout << "biting_rate:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i]->tau.first;
    }
    std::cout << std::endl;

    std::cout << "rprev:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hPrevalence[i];
    }
    std::cout << std::endl;
  
    std::cout << "vprev:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vPrevalence[i];
    }
    std::cout << std::endl;
  
    std::cout << "vmu:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i]->mu.first;
    }
    std::cout << std::endl;

    if (params->global->estimateXi) {
      std::cout << "bvector:";
      for (size_t i = 0; i < params->vectors.size(); ++i) {
        std::cout << " " << params->vectors[i]->b.first;
      }
      std::cout << std::endl;
    } else {  
      std::cout << "xi:";
      for (size_t i = 0; i < params->vectors.size(); ++i) {
        std::cout << " " << params->vectors[i]->xi.first;
      }
      std::cout << std::endl;
    }
  }

  size_t nvars =
    params->hosts.size() + params->groups.size() * params->vectors.size() +
    params->vectors.size();

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
  if (params->global->estimateXi) {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      gsl_vector_set(x_init, params->hosts.size() + v, .5);
    }
  } else {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      gsl_vector_set(x_init, params->hosts.size() + v, 1);
    }
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      gsl_vector_set(x_init, j + v * params->groups.size() +
                     params->hosts.size() + params->vectors.size(),
                     params->vPrevalence[v]);
    }
  }

  // if (jac) {
  //   gsl_multiroot_fdfsolver_set(sdf, &fdf, x_init);
  // }
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
        gsl_vector_get(sol, j + v * params->groups.size() +
                       params->hosts.size() + params->vectors.size());
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
