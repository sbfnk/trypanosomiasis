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

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

enum SamplingType {Normal, Linear, Log};

struct Parameter
{
  Parameter() :
    sampling(Normal), limits(std::make_pair(-1.,-1.)) {setMean(0);}
  Parameter(double v) :
    sampling(Normal), limits(std::make_pair(-1.,-1.)) {setMean(v);}
  Parameter(double v, std::pair<double, double> l) :
    sampling(Normal), limits(l) {setMean(v);}
             
  SamplingType sampling;
  std::pair<double, double> limits;
  double value;

  void setMean(double m) { mean = m; value = mean;}
  double getMean() const { return mean; }

private:
  
  double mean;
};

namespace po = boost::program_options;

/*! Base class Class for (host/vector/system) parameters */
class ParamContainer
{
public:

  struct ParamInfo {
    ParamInfo(std::string o, std::string d, Parameter* p) :
      option(o), description(d), param(p) {;}
    
    std::string option;
    std::string description;
    Parameter* param;
  };

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
  
  //! Accessor for params
  const std::vector<ParamInfo>& getParams() const
  { return params; }

protected:
    
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

  void NormaliseHabitats();
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
  Parameter M, N, mu, gamma, b, f, n;
  Parameter x0; // initial number of infected hosts
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
    
  Parameter M, N, mu, tau, b, xi, alpha;
  std::vector<double> groupPrev;
};

struct Group
{
  Group() :
    f(.0)
  {}
  
  Group(size_t initMember) :
    f(.0),
    members(std::vector<size_t>(1, initMember))
  {}
  
  double f;
  std::vector<size_t> members;
  int x0; // initial number of infected vectors in group
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
  bool xiTau; //!< set xi to tau
  bool teneralOnly; //!< only consider teneral flies to be infected
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

struct Event {
  int add;
  int remove;
  double rate;
};

int getSeed();

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
        std::replace(name.begin(), name.end(), ' ', '.');
      } else {
        if ((header[i] == it->option)) {
          if (data[i].length() > 0) {
            std::istringstream s(data[i]);
            double m;
            s >> m;
            it->param->setMean(m);
          }
        } else if (header[i] == (it->option + "_low")) {
          if (data[i].length() > 0) {
            std::istringstream s(data[i]);
            s >> it->param->limits.first;
          }
        } else if (header[i] == (it->option + "_high")) {
          if (data[i].length() > 0) {
            std::istringstream s(data[i]);
            s >> it->param->limits.second;
          }
        } else if (header[i] == (it->option + "_sampling")) {
          if (data[i].length() > 0) {
            if (data[i] == "l") {
              it->param->sampling = Linear;
            } else if (data[i] == "o") {
              it->param->sampling = Log;
            } else if (data[i] == "p") {
              it->param->sampling = Normal;
            } else if (data[i] != "") {
              std::cerr << "WARNING: Unknown sampling type for " << name << ": "
                        << data[i] << std::endl;
            }
          }
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
      it->param->setMean(vm[paramPrefix + it->option].as<double>());
    } else if (vm.count(paramPrefix + it->option + "_low")){
      it->param->limits.first = vm[paramPrefix + it->option].as<double>();
    } else if (vm.count(paramPrefix + it->option + "_high")){
      it->param->limits.second = vm[paramPrefix + it->option].as<double>();
    }
  }
}

/*! \brief Normalise habitat contributions
   
Normalises all habitat contributions to be 
  
\param[in] vm The map of command line parameters
*/
void HabitatContainer::NormaliseHabitats()
{
  double habitatSum = .0; //!< sum of all habitat contributions (for
                          //!normalisation) 
  for (size_t i = 0; i < habitat.size(); ++i) {
    habitatSum += habitat[i].value;
  }
  for (std::vector<Parameter>::iterator it = habitat.begin();
       it != habitat.end(); it++) {
    (*it).value /= habitatSum;
    (*it).limits.first /= habitatSum;
    (*it).limits.second /= habitatSum;
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
      if (numPos == 1); {
        std::istringstream numStr(header[i].substr(numPos));
        numStr >> index;
        std::stringstream paramStr;
        paramStr << "X" << index;
        while (index > habitat.size()) {
          habitat.push_back(Parameter(.0, std::make_pair(.0, .0)));
        }
        if (header[i] == paramStr.str()) {
          double m;
          s >> m;
          habitat[index-1].setMean(m);
        } else if (header[i] == (paramStr.str() + "_low")) {
          s >> habitat[index-1].limits.first;
        } else if (header[i] == (paramStr.str() + "_high")) {
          s >> habitat[index-1].limits.second;
        }
      }
    }
  }

  for (size_t i = 0; i < habitat.size(); ++i) {
    if (habitat[i].value > 0) {
      std::stringstream paramStr;
      paramStr << "X" << i;
      std::stringstream describStr;
      describStr << "Density in habitat " << i;
      params.push_back(ParamInfo(paramStr.str(), describStr.str(), &(habitat[i])));
    }
  }
}

void HabitatContainer::ReadParams(po::variables_map const &vm) 
{
  ParamContainer::ReadParams(vm);
  for (size_t i = 0; i < habitat.size(); ++i) {
    std::stringstream ss;
    ss << i;
    if (vm.count(name + "-X" + ss.str())) {
      // command line parameter has been specified, assign to model variable
      habitat[i].setMean(vm[name + "-X" + ss.str()].as<double>());
    } else if (vm.count(name + "-X" + ss.str() + "_low")){
      habitat[i].limits.first = vm[name + "-X" + ss.str() + "_low"].as<double>();
    } else if (vm.count(name + "-X" + ss.str() + "_high")){
      habitat[i].limits.second = vm[name + "-X" + ss.str() + "_low"].as<double>();
    }
  }
}

void GlobalParams::ReadParams(po::variables_map const &vm)
{
  ParamContainer::ReadParams(vm);
  if (vm.count("estimate-xi")) {
    estimateXi = true;
  }
  if (vm.count("teneral-only")) {
    teneralOnly = true;
  }
  if (vm.count("xi-tau")) {
    xiTau = true;
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
  HabitatContainer("host"),
  x0(-1)
{
  params.push_back(ParamInfo("N", "Population size", &N));
  params.push_back(ParamInfo("M", "Number infected", &M));
  params.push_back(ParamInfo("mu", "Mortality rate", &mu));
  params.push_back(ParamInfo("gamma", "Recovery rate", &gamma));
  params.push_back(ParamInfo("b", "Susceptibility", &b));
  params.push_back(ParamInfo("n", "Abundance", &n));
  params.push_back(ParamInfo("f", "Biting preference", &f));
  params.push_back(ParamInfo("x0", "Initial number of infected", &x0));
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
  params.push_back(ParamInfo("alpha", "Rate of incubation", &alpha));
}
  
GlobalParams::GlobalParams() :
  ParamContainer("global"),
  estimateXi(false),
  xiTau(false),
  teneralOnly(false)
{
  options->add_options()
    ("estimate-xi,x", "Estimate xi (with vector susceptibility set)")
    ("xi-tau", "set xi to tau")
    ("teneral-only,y", "Only consider teneral flies")
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
    hPrevalence[i] = hosts[i]->M.value / hosts[i]->N.value;
  }
    
  for (size_t v = 0; v < vectors.size(); ++v) {
    vPrevalence[v] = vectors[v]->M.value / vectors[v]->N.value;
  }
    
  //!< amount of overlap between host habitats
  habitatOverlap = std::vector<std::vector<double> >
    (groups.size(), std::vector<double>(groups.size(), .0));

  size_t nSpecies = 0;
  for (size_t j = 0; j < groups.size(); ++j) {
    nSpecies += groups[j].members.size();
  }

  std::vector<double> habitat_densities(hosts[0]->habitat.size(), .0);
  std::vector<double> host_densities(nSpecies, .0);

  if (global->habType == "f") {
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t k = 0; k < groups[j].members.size(); ++k) {
        size_t i = groups[j].members[k];
        for (size_t o = 0; o < hosts[0]->habitat.size(); ++o) {
          // std::cout << "Host " << i << " Habitat " << o << " " << hosts[i]->habitat[o].value << std::endl;
          host_densities[i] += hosts[i]->habitat[o].value;
          habitat_densities[o] += hosts[i]->habitat[o].value;
        }
      }
    }
    // for (size_t j = 0; j < groups.size(); ++j) {
    //   for (size_t k = 0; k < groups[j].members.size(); ++k) {
    //     size_t i = groups[j].members[k];
    //     std::cout << "Host " << i << " " << host_densities[i] << std::endl;
    //     for (size_t o = 0; o < hosts[0]->habitat.size(); ++o) {
    //       std::cout << "Habitat " << o << " " << habitat_densities[o] << std::endl;
    //     }
    //   }
    // }
  }

  if (global->habType == "b" ||
      global->habType == "f") {
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
        for (size_t k = 0; k < groups[j].members.size(); ++k) {
          size_t i = groups[j].members[k];
          for (size_t n = 0; n < groups[m].members.size(); ++n) {
            size_t l = groups[m].members[n];
            for (size_t o = 0; o < hosts[0]->habitat.size(); ++o) {
              if (global->habType == "b") {
                habitatOverlap[j][m] = 1;
                habitatOverlap[m][j] = 1;
              } else {
                habitatOverlap[j][m] +=
                  hosts[i]->habitat[o].value / host_densities[i] *
                  hosts[l]->habitat[o].value / host_densities[l];
                if (j != m) {
                  habitatOverlap[m][j] += 
                    hosts[l]->habitat[o].value / host_densities[l] *
                    hosts[i]->habitat[o].value / host_densities[i];
                }
              }
            }
          }
        }
      }
    }

//    for (size_t j = 0; j < groups.size(); ++j) {
//      for (size_t m = 0; m < groups.size(); ++m) {
//        std::cout << j << " " << m << " " << habitatOverlap[j][m] 
//                  << std::endl;
//      }
//    }
                  
    
    // for (size_t j = 0; j < groups.size(); ++j) {
    //   for (size_t k = 0; k < groups[j].members.size(); ++k) {
    //     size_t i = groups[j].members[k];
    //     hosts[i]->NormaliseHabitats();
    //   }
    // }
  } else {
    for (size_t j = 0; j < groups.size(); ++j) {
      for (size_t m = j; m < groups.size(); ++m) {
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
    (params->groups.size(), std::vector<double>(params->vectors.size(), .0));

  for (size_t i = 0; i < params->hosts.size(); ++i) {
    bhost[i] = gsl_vector_get(x, i);
    // std::cout << "bhost[" << i << "]: " << bhost[i] << std::endl;
  }
  if (params->global->estimateXi) {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      bvector[v] = params->vectors[v]->b.value;
      xi[v] = gsl_vector_get(x, v + params->hosts.size());
      // std::cout << "bvector[" << v << "]: " << bvector[v] << std::endl;
      // std::cout << "xi[" << v << "]: " << xi[v] << std::endl;
    }
  } else {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      bvector[v] = gsl_vector_get(x, v + params->hosts.size());
      xi[v] = params->vectors[v]->xi.value;
      // std::cout << "bvector[" << v << "]: " << bvector[v] << std::endl;
      // std::cout << "xi[" << v << "]: " << xi[v] << std::endl;
    }
  }
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      pv[j][v] =
        gsl_vector_get(x, j + v * params->groups.size() +
                       params->hosts.size() + params->vectors.size());
      // std::cout << "pv[" << j << "," << v << "]: " << pv[j][v] << std::endl;
    }
  }

  std::vector<std::vector<double> > yv
    (params->groups.size(), std::vector<double>(params->vectors.size(), .0));
  std::vector<double> yh(params->hosts.size(), .0);

  for (size_t j = 0; j < params->groups.size(); ++j) {
    for (size_t v = 0; v < params->vectors.size(); ++v) {
      double enumerator = .0;
      double denominator = .0;
      for (size_t k = 0; k < params->groups.size(); ++k) {
        enumerator +=
          params->groups[k].f * params->habitatOverlap[j][k] * pv[k][v];
        denominator +=
          params->groups[k].f * params->habitatOverlap[j][k];
      }
      if (params->global->teneralOnly) {
        yv[j][v] = (params->vectors[v]->mu.value +
                    params->vectors[v]->tau.value) *
          (pv[j][v] + xi[v] / params->vectors[v]->mu.value *
           (pv[j][v] - enumerator/denominator));
      } else {
        yv[j][v] = (pv[j][v] * params->vectors[v]->mu.value + xi[v] *
                    (pv[j][v] - enumerator/denominator)) /
          (1 - pv[j][v]);
      }
    }
    for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
      size_t i = params->groups[j].members[k];
       yh[i] = params->hPrevalence[i] / (1 - params->hPrevalence[i]) *
        (params->hosts[i]->mu.value + params->hosts[i]->gamma.value);
    }
  }
  
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        double i_av;
        if (params->vectors[v]->alpha.value > 0) {
          i_av = params->vectors[v]->alpha.value /
            (params->vectors[v]->alpha.value + params->vectors[v]->mu.value) *
            ((params->vectors[v]->alpha.value + params->vectors[v]->mu.value) /
             (params->vectors[v]->alpha.value + params->vectors[v]->mu.value + xi[v]) *
             pv[j][v] +
             xi[v] /
             (params->vectors[v]->alpha.value + params->vectors[v]->mu.value + xi[v]) *
             params->vPrevalence[v]);
        } else {
          i_av = pv[j][v];
        }
        yh[i] -= bhost[i] * params->vectors[v]->tau.value *
          params->hosts[i]->f.value / params->hosts[i]->n.value * i_av;
        yv[j][v] -= bvector[v] * params->vectors[v]->tau.value *
          params->hosts[i]->f.value * params->hPrevalence[i] /
          params->groups[j].f;
        // std::cout << k << " " << i << " " << j << " " << params->vectors[v]->tau.value
        //           << " " << params->groups[j].f << " " << bvector[v]
        //           << " " << bhost[i] << " " << params->hosts[i]->f.value
        //           << " " << params->hPrevalence[i] << " " << pv[j][v]
        //           << " " << params->vectors[v]->mu.value << " "
        //           << xi[v] << " " << yv[j][v] << " " << yh[i]
        //           << " " << i_av << std::endl;
      }
    }
  }
    
  for (size_t v = 0; v < params->vectors.size(); ++v) {
    for (size_t j = 0; j < params->groups.size(); ++j) {
      for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
        size_t i = params->groups[j].members[k];
        gsl_vector_set(f, i, yh[i]);
        // std::cout << "yh[" << i << "]: " << yh[i] << std::endl;
      }
      gsl_vector_set(f, j + v * params->groups.size() +
                     params->hosts.size(), yv[j][v]);
      // std::cout << "yv[" << j << "," << v << "]: " << yv[j][v] << std::endl;
      weightedVectorPrevSum[v] += pv[j][v] * params->groups[j].f;
    }
    gsl_vector_set(f, params->hosts.size() +
                   params->groups.size() * params->vectors.size() + v,
                   params->vPrevalence[v] - weightedVectorPrevSum[v]);
    // std::cout << "y[" << v << "]: "
    //           << (params->vPrevalence[v] - weightedVectorPrevSum[v])
    //           << std::endl;
  }

    
  return GSL_SUCCESS;
}

// function to find root of (derivative)
int betafunc_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  // betafunc_params* params = ((struct betafunc_params*) p);

  // std::vector<double> bhost(params->hosts.size(), .0);
  // std::vector<double> bvector(params->vectors.size(), .0);
  // std::vector<std::vector<double> > pv
  //   (params->groups.size(), std::vector<double>(params->vectors.size(), .0));
  // std::vector<double> xi(params->vectors.size(), .0);
  
  // for (size_t i = 0; i < params->hosts.size() +
  //        params->vectors.size() * params->groups.size() + params->vectors.size();
  //      ++i) {
  //   for (size_t j = 0; j < params->hosts.size() + 
  //          params->vectors.size() * params->groups.size() + params->vectors.size();
  //        ++j) {
  //     gsl_matrix_set(J, i, j, 0);
  //   }
  // }
  // for (size_t i = 0; i < params->hosts.size(); ++i) {
  //   bhost[i] = gsl_vector_get(x, i);
  //   std::cout << "bhost[" << i << "]: " << bhost[i] << std::endl;
  // }
  // if (params->global->estimateXi) {
  //   for (size_t v = 0; v < params->vectors.size(); ++v) {
  //     bvector[v] = params->vectors[v]->b.value;
  //     xi[v] = gsl_vector_get(x, v + params->hosts.size());
  //     std::cout << "bvector[" << v << "]: " << bvector[v] << std::endl;
  //     std::cout << "xi[" << v << "]: " << xi[v] << std::endl;
  //   }
  // } else {
  //   for (size_t v = 0; v < params->vectors.size(); ++v) {
  //     bvector[v] = gsl_vector_get(x, v + params->hosts.size());
  //     xi[v] = params->vectors[v]->xi.value;
  //     std::cout << "bvector[" << v << "]: " << bvector[v] << std::endl;
  //     std::cout << "xi[" << v << "]: " << xi[v] << std::endl;
  //   }
  // }
  // for (size_t v = 0; v < params->vectors.size(); ++v) {
  //   for (size_t j = 0; j < params->groups.size(); ++j) {
  //     pv[j][v] =
  //       gsl_vector_get(x, j + v * params->groups.size() +
  //                      params->hosts.size() + params->vectors.size());
  //     std::cout << "pv[" << j << "," << v << "]: " << pv[j][v] << std::endl;
  //   }
  // }

  // for (size_t j = 0; j < params->groups.size(); ++j) {
  //   for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
  //     size_t i = params->groups[j].members[k];
  //     double dfidbi = 0;
  //     double dfidpvjv = 0;
  //     for (size_t v = 0; v < params->vectors.size(); ++v) {
  //       dfidbi -= params->vectors[v]->tau.value *
  //         params->hosts[i]->f.value / params->groups[j].f /
  //         params->hosts[i]->n.value * pv[j][v];
  //       dfidpvjv = -bhost[i] * params->vectors[v]->tau.value *
  //         params->hosts[i]->f.value / params->groups[j].f /
  //         params->hosts[i]->n.value;
  //       gsl_matrix_set(J, i, j + v * params->groups.size() +
  //                      params->hosts.size() + params->vectors.size(),
  //                      dfidpvjv);
  //     }
  //     gsl_matrix_set(J, i, i, dfidbi);
  //   }
  // }

  // for (size_t v = 0; v < params->vectors.size(); ++v) {
  //   for (size_t j = 0; j < params->groups.size(); ++j) {
  //     double enumerator = .0;
  //     double denominator = .0;
  //     for (size_t l = 0; l < params->groups.size(); ++l) {
  //       enumerator += pv[l][v] * params->habitatOverlap[j][l] *
  //         params->groups[l].f;
  //       denominator += params->habitatOverlap[j][l] *
  //         params->groups[l].f;
  //     }
  //     if (params->global->estimateXi) {
  //       double dgjvdxiv = (pv[j][v] - enumerator/denominator) / (1 - pv[j][v]);
  //       gsl_matrix_set(J, j + v * params->groups.size() + params->hosts.size(),
  //                      v + params->hosts.size(),
  //                      dgjvdxiv);
  //     } else {
  //       double dgjvdbv = 0;
  //       for (size_t k = 0; k < params->groups[j].members.size(); ++k) {
  //         size_t i = params->groups[j].members[k];
  //         dgjvdbv -= params->vectors[v]->tau.value *
  //           params->hosts[i]->f.value / params->groups[j].f *
  //           params->hPrevalence[i];
  //       }
  //       gsl_matrix_set(J, j + v * params->groups.size() + params->hosts.size(),
  //                      v + params->hosts.size(),
  //                      dgjvdbv);
  //     }
  //     double dgjvdpvjv = 1 / pow(1-pv[j][v], 2) *
  //       (params->vectors[v]->mu.value +
  //        (xi[j] * (pv[j][v] - 1) * params->groups[j].f - enumerator) / denominator);
  //     gsl_matrix_set(J, j + v * params->groups.size() + params->hosts.size(),
  //                    j + v * params->groups.size() + params->hosts.size() +
  //                    params->vectors.size(),
  //                    dgjvdpvjv);
  //     for (size_t l = 0; l < params->groups.size(); ++l) {
  //       if (j != l) {
  //         double dgjvdpvlv = -xi[v] / (1 - pv[j][v]) * params->habitatOverlap[j][l] /
  //           denominator;
  //         gsl_matrix_set(J, j + v * params->groups.size() + params->hosts.size(),
  //                        l + v * params->groups.size() + params->hosts.size() +
  //                        params->vectors.size(),
  //                        dgjvdpvlv);
  //       }
  //     }
  //   }
  // }

  // for (size_t v = 0; v < params->vectors.size(); ++v) {
  //   for (size_t j = 0; j < params->groups.size(); ++j) {
  //     double dhvdpjv = -params->groups[j].f;
  //         gsl_matrix_set(J, j + v * params->groups.size() + params->hosts.size() +
  //                        params->vectors.size(),
  //                        j + v * params->groups.size() + params->hosts.size() +
  //                        params->vectors.size(),
  //                        dhvdpjv);
  //   }
  // }

  // for (size_t i = 0; i < params->hosts.size() + params->vectors.size() +
  //        params->groups.size()*params->vectors.size(); ++i) {
  //   for (size_t j = 0; j < params->hosts.size() + params->vectors.size() +
  //          params->groups.size()*params->vectors.size(); ++j) {
  //     std::cout << "J(" << i << "," << j << ") = " << gsl_matrix_get(J, i, j)
  //               << std::endl;
  //   }
  // }
  
  return GSL_SUCCESS;
}

// function to find root of and derivative
int betafunc_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J)
{
  betafunc_f(x, p, f);
  betafunc_df(x, p, J);

  return GSL_SUCCESS;
}

// find beta (and alpha and p^v_i) from forces of infection
int betaffoiv(void *p, std::vector<double> &vars,
              bool jac = false, unsigned int verbose = 0)
{
  betafunc_params* params = ((struct betafunc_params*) p);

  if (verbose) {
    std::cout << "N:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->N.value;
    }
    std::cout << std::endl;

    std::cout << "M:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->M.value;
    }
    std::cout << std::endl;

    std::cout << "rabundance:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->n.value;
    }
    std::cout << std::endl;
  
    std::cout << "mu:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->mu.value;
    }
    std::cout << std::endl;
  
    std::cout << "gamma:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->gamma.value;
    }
    std::cout << std::endl;
  
    std::cout << "biting:";
    for (size_t i = 0; i < params->hosts.size(); ++i) {
      std::cout << " " << params->hosts[i]->f.value;
    }
    std::cout << std::endl;
  
    std::cout << "biting_rate:";
    for (size_t i = 0; i < params->vectors.size(); ++i) {
      std::cout << " " << params->vectors[i]->tau.value;
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
      std::cout << " " << params->vectors[i]->mu.value;
    }
    std::cout << std::endl;

    if (params->global->estimateXi) {
      std::cout << "bvector:";
      for (size_t i = 0; i < params->vectors.size(); ++i) {
        std::cout << " " << params->vectors[i]->b.value;
      }
      std::cout << std::endl;
    } else {  
      std::cout << "xi:";
      for (size_t i = 0; i < params->vectors.size(); ++i) {
        std::cout << " " << params->vectors[i]->xi.value;
      }
      std::cout << std::endl;
    }
  }

  size_t nvars =
    params->hosts.size() + params->groups.size() * params->vectors.size() +
    params->vectors.size();

  const gsl_multiroot_fdfsolver_type * Tdf;
  gsl_multiroot_fdfsolver * sdf;
  gsl_multiroot_function_fdf fdf;
  if (jac) {
    Tdf = gsl_multiroot_fdfsolver_hybridsj;
    sdf = gsl_multiroot_fdfsolver_alloc (Tdf, nvars);
    fdf.f = &betafunc_f;
    fdf.df = &betafunc_df;
    fdf.fdf = &betafunc_fdf;
    fdf.n = nvars;
    fdf.params = p;
  }

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
                     0.1);
                     // params->vPrevalence[v]);
    }
  }

  if (jac) {
    gsl_multiroot_fdfsolver_set(sdf, &fdf, x_init);
  }
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
  
  if (jac) {
    gsl_multiroot_fdfsolver_free (sdf);
  }
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x_init);

  return status;
}

int getSeed()
{
  std::ifstream rand("/dev/urandom");
  char tmp[sizeof(int)];
  rand.read(tmp,sizeof(int));
  rand.close();
  int* number = reinterpret_cast<int*>(tmp);
  return (*number);
}

#endif
