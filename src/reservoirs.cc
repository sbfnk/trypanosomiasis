// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file reservoirs.cc
  \brief Main program for reservoir analysis based on prevalence data

  This takes prevalence data of some disease for a number of host species, as
  well as data on abundance, vector biting preference, mortality and recovery
  rates to estimate the relative contribution of each host species to the
  next-generation matrix, and, consequently, to the basic reproductive number
  R_0. This is done assuming the prevalences measured reflect an equilibrium of
  the multi-host system. Some parts of the code are specific to analysing data
  for African Trypanosomiasis, but this can be easily generalised.
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <armadillo>
#include <math.h>
#include <sys/time.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/algorithm/string.hpp>

#include "reservoirs.hh"
#include "ihs.hh"

/*! Main program. */
int main(int argc, char* argv[])
{

  std::string dataFile; //!< name of the file holding host parameters
  std::string vectorFile; //!< name of the file holding vector parameters
  std::string outFile; //!< output file

  bool jacobian = false; //!< use jacobian?
  size_t lhsSamples = 0; //!< number of samples to be used in linear hypercube
                         //!< sampling 
  bool calcBetas = true; //!< calculate the betas (or use analytical formula in
                          //!< simple model?)

  size_t samples = 0; //!< number of samples to be used from negative binomial
                      //!< distribution around measured prevalence

  unsigned int verbose = 0; //!< be verbose

  std::vector<Group> groups; //!< composition of groups
  std::vector<size_t> hostGroups; //!< per-host group membership

  std::string addColumn; //!< entry of custom columns
  std::string addHeader; //!< header of custom columns
  size_t attempts = 1; //!< number of attempts for solving equations

  bool sim = false; // simulate?
  double recordStep = 1.;
  double tmax = 100.;
  size_t humanmax = 0;

  // main options
  po::options_description main_options
    ("Usage: reservoirs [options]... \n\nOptions");

  main_options.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("data-file,f", po::value<std::string>()->default_value("data/bipindi.csv"),
     "data file")
    ("vector-file,e", po::value<std::string>()->
     default_value("data/bipindi_vector.csv"),
     "vector file")
    ("output-file,o", po::value<std::string>()->default_value("-"),
     "output file (\"-\" for output to stdout")
    ("groups,g", po::value<std::string>(),
     "groups (semicolon-separated list of comma-separated hosts")
    ("lhs,l", po::value<size_t>()->default_value(0),
     "number of latin hypercube samples (0 for no sampling)")
    ("samples,s", po::value<size_t>()->default_value(0),
     "number of samples (0 for no sampling)")
    ("normal,r", 
     "assume normally distributed parameters (as opposed to linear ones)")
    ("jacobian,j", 
     "use jacobian")
    ("calculate,u", 
     "calulate parameters, don't estimate")
    ("noheader,d", 
     "do not print header")
    ("print,p", 
     "print results")
    ("addcolumn,c", po::value<std::string>(),
     "custom column(s) to add")
    ("addheader,C", po::value<std::string>(),
     "header for custom columns")
    ("attempts,m", po::value<size_t>()->default_value(1),
     "number of solving attempts")
    ("sim",
     "simulate")
    ("timestep,t", po::value<double>()->default_value(1.),
     "timestep of recording data")
    ("tmax", po::value<double>()->default_value(100),
     "how long to run")
    ("humanmax", po::value<size_t>(),
     "human level to stop at")
    ;

  // read options
  po::variables_map vm;
    
  try {
    po::store(po::command_line_parser(argc, argv).options(main_options).
              allow_unregistered().run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);

  // help option
  if (vm.count("help")) {
    std::cout << main_options << std::endl;
    return 0;
  }

  // verbose option
  if (vm.count("verbose")) {
    verbose = 1;
  }
  if (vm.count("very-verbose")) {
    verbose = 2;
  }
  
  // host data file, necessary option
  if (vm.count("data-file")) {
    dataFile = vm["data-file"].as<std::string>();
  } else {
    std::cerr << "Error: must specify data file" << std::endl;
    return 1;
  }

  // vector data file, necessary option
  if (vm.count("vector-file")) {
    vectorFile = vm["vector-file"].as<std::string>();
  } else {
    std::cerr << "Error: must specify vector file" << std::endl;
    return 1;
  }

  samples = vm["samples"].as<size_t>();

  // vector data file, necessary option if more than one sample is generated
  if (vm.count("output-file")) {
    outFile = vm["output-file"].as<std::string>();
  } else {
    outFile = "-";
  }

  if (vm.count("jacobian")) {
    jacobian = true;
  }

  if (vm.count("lhs")) {
    lhsSamples = vm["lhs"].as<size_t>();
  }

  if (vm.count("calculate")) {
    calcBetas = false;
  }

  if (vm.count("sim")) {
    sim = true;
    recordStep = vm["timestep"].as<double>();
    tmax = vm["tmax"].as<double>();
    if (vm.count("humanmax")) {
      humanmax = vm["humanmax"].as<size_t>();
    }
  }

  if (vm.count("addcolumn")) {
    addColumn = vm["addcolumn"].as<std::string>();
  }

  if (vm.count("addheader")) {
    addHeader = vm["addheader"].as<std::string>();
  }

  attempts = vm["attempts"].as<size_t>();

  std::vector<Host*> hosts; // vector of hosts
  std::vector<Vector*> vectors; // vector of vectors

  // tokenizer for reading csv file
  typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
  boost::escaped_list_separator<char> sep('\\', ',', '\"');

  std::vector<std::string> headings;
  bool firstLine = true;

  GlobalParams* global = new GlobalParams();
  main_options.add(*(global->getOptions()));

  // ********************** read data file ********************
  
  std::ifstream in(dataFile.c_str());
  if (!in.is_open()) {
    std::cerr << "Could not open " << dataFile << std::endl;
    return 1;
  }

  std::string line;
  while (std::getline(in,line))  {
    std::vector<std::string> lineVector;
    Tokenizer tok(line, sep);
    lineVector.assign(tok.begin(), tok.end());

    if (firstLine) {
      headings = lineVector;
      firstLine = false;
    } else {
      Host* newHost = new Host();
      newHost->ReadTable(lineVector, headings);
      main_options.add(*(newHost->getOptions()));
      hosts.push_back(newHost);
    }
  }
  in.close();

  // ************** read group composition option ****************
  // for the moment, works only with one vector species
    
  if (vm.count("groups")) {
    std::string groupString = vm["groups"].as<std::string>();
    std::vector<std::string> result;
    boost::algorithm::split(result, groupString, boost::is_any_of(";"));
    for (std::vector<std::string>::iterator it = result.begin();
         it != result.end(); it++) {
      std::vector<std::string> subres;
      boost::algorithm::split(subres, *it, boost::is_any_of(","));
      Group newGroup;
      for (std::vector<std::string>::iterator it2 = subres.begin();
           it2 != subres.end(); it2++) {
        size_t s;
        std::istringstream myStream(*it2);
        myStream >> s;
        newGroup.members.push_back(s);
      }
      groups.push_back(newGroup);
    }
  } else {
    groups.push_back(Group());
    for (size_t i = 0; i < hosts.size(); ++i) {
      groups[0].members.push_back(i);
    }
  }

  hostGroups.resize(hosts.size());
  for (size_t i = 0; i < groups.size(); ++i) {
    for (size_t j = 0; j < groups[i].members.size(); ++j) {
      hostGroups[groups[i].members[j]] = i;
    }
  }

  po::options_description group_options("\nInitial number of infected vectors");
  for (size_t i = 0; i < groups.size(); ++i) {
    std::stringstream group_option;
    group_option << "g" << i << "-x0";
    group_options.add_options()
      (group_option.str().c_str(), po::value<double>()->default_value(-1),
       "Number of infected vectors in group "+i)
      ;
  }
  main_options.add(group_options);
  
  if (verbose) {
    std::cout << "Number of groups: " << groups.size() << std::endl;
    for (size_t i = 0; i < groups.size(); ++i) {
      std::cout << "  group #" << i << ":";
      for (size_t j = 0; j < groups[i].members.size(); ++j) {
        std::cout << " " << groups[i].members[j];
      }
      std::cout << std::endl;
    }
  }
    
  // ********************** read vector file ********************
  
  in.open(vectorFile.c_str());
  firstLine = true;
  
  if (!in.is_open()) {
    std::cerr << "Could not open " << vectorFile << std::endl;
    return 1;
  }

  while (std::getline(in,line))  {
    std::vector<std::string> lineVector;
    Tokenizer tok(line, sep);
    lineVector.assign(tok.begin(), tok.end());

    if (firstLine) {
      headings = lineVector;
      firstLine = false;
    } else {
      Vector* newVector = new Vector();
      newVector->ReadTable(lineVector, headings);
      main_options.add(*(newVector->getOptions()));
      newVector->groupPrev.resize(groups.size(), .0);
      vectors.push_back(newVector);
    }
  }
  in.close();

  // ************** read additional parameters ****************

  try {
    po::store(po::command_line_parser(argc, argv).
              options(main_options).run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::cout << main_options << std::endl;
    return 1;
  }
  po::notify(vm);

  global->ReadParams(vm);
  for (size_t i = 0; i < hosts.size(); ++i) {
    hosts[i]->ReadParams(vm);
  }
  for (size_t i = 0; i < vectors.size(); ++i) {
    vectors[i]->ReadParams(vm);
  }

  // ********************* estimate betas *********************

  std::streambuf * buf;
  std::ofstream of;
  if (outFile == "-") {
    buf = std::cout.rdbuf();
  } else {
    of.open(outFile.c_str());
    buf = of.rdbuf();
  }
  
  std::ostream out(buf);

  int seed = getSeed();
  int *x = 0;
  
  boost::mt19937 gen(seed);
  
  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    randGen(gen, boost::uniform_real<> (0,1));

  if (!(vm.count("noheader") ||
        (outFile == "-" && vm.count("print")) ||
        sim
      )) {
    if (samples > 0) {
      out << "n";
      if (addHeader.length() > 0) {
        out << ",";
      }
    }
    if (addHeader.length() > 0) {
      out << addHeader;
    }

    for (size_t j = 0; j < hosts.size(); ++j) {
      out << ",\"" << hosts[j]->getName() << "_prev\"";
      if (lhsSamples > 0) {
        for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
          out << ",\"" << hosts[j]->getName() << "_"
              << hosts[j]->getParams()[k].option << "\"";
        }
      }
    }
    for (size_t v = 0; v < vectors.size(); ++v) {
      out << ",\"" << vectors[v]->getName() << "_prev\"";
      if (groups.size() > 1) {
        for (size_t i = 0; i < groups.size(); ++i) {
          out << ",\"" << vectors[v]->getName() << "_" << i+1 << "_prev\"";
        }
      }
      if (lhsSamples > 0) {
        for (size_t k = 0; k < vectors[v]->getParams().size(); ++k) {
          out << ",\"" << vectors[v]->getName() << "_"
              << vectors[v]->getParams()[k].option << "\"";
        }
      }
    }
    for (size_t j = 0; j < hosts.size(); ++j) {
      out << ",\"" << hosts[j]->getName() << "\"";
    }
    out << ",\"domestic\",\"wildlife\",\"domestic+wildlife\","
        << "\"human+domestic\",\"human+wildlife\",\"R0\"";
    out << std::endl;
  }

  boost::math::normal_distribution<> stdNormal
    (0,1);
  size_t i = 0;
  size_t sample = 0;
  do {

    size_t nParams = 0;
    for (size_t j = 0; j < hosts.size(); ++j) {
      for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
        if (hosts[j]->getParams()[k].param->limits.first >= 0 &&
            hosts[j]->getParams()[k].param->limits.second >= 0) {
          ++nParams;
        }
      }
    }
    for (size_t j = 0; j < vectors.size(); ++j) {
      for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
        if (vectors[j]->getParams()[k].param->limits.first >= 0 &&
            vectors[j]->getParams()[k].param->limits.second >= 0) {
          ++nParams;
        }
      }
    }

    if (lhsSamples > 0) {
      if (i % lhsSamples == 0) {
        // generate latin hypercube samples
        x = (int*) malloc (nParams * lhsSamples * sizeof(int));
        ihs(nParams, lhsSamples, 5, &seed, x);
        sample = 0;
      }
      
      for (size_t j = 0; j < hosts.size(); ++j) {
        for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
          // std::cout << hosts[j]->getName() << " "
          //           << hosts[j]->getParams()[k].option << " " 
          //           << hosts[j]->getParams()[k].param->getMean() << " "
          //           << hosts[j]->getParams()[k].param->value << " "
          //           << hosts[j]->getParams()[k].param->limits.first << " " 
          //           << hosts[j]->getParams()[k].param->limits.second
          //           << std::endl;
          if (hosts[j]->getParams()[k].param->limits.first >= 0 &&
              hosts[j]->getParams()[k].param->limits.second >= 0) {
            if (hosts[j]->getParams()[k].param->sampling == Linear) {
              hosts[j]->getParams()[k].param->value =
                hosts[j]->getParams()[k].param->limits.first +
                (hosts[j]->getParams()[k].param->limits.second -
                 hosts[j]->getParams()[k].param->limits.first) *
                (x[sample] - 1) / static_cast<double>(lhsSamples - 1);
            } else if (hosts[j]->getParams()[k].param->sampling == Log) {
              hosts[j]->getParams()[k].param->value =
                hosts[j]->getParams()[k].param->limits.first +
                (hosts[j]->getParams()[k].param->limits.second -
                 hosts[j]->getParams()[k].param->limits.first) *
                (exp(10* (x[sample] - 1) / static_cast<double>(lhsSamples - 1)) - 1) /
                (exp(10) - 1);
            } else if (hosts[j]->getParams()[k].param->sampling == Normal) {
              double sd = fabs(hosts[j]->getParams()[k].param->limits.second -
                              hosts[j]->getParams()[k].param->limits.first) / 2.;
              hosts[j]->getParams()[k].param->value =
                hosts[j]->getParams()[k].param->getMean() +
                sd * quantile(stdNormal,
                              x[sample] / static_cast<double>(lhsSamples + 1) +
                              cdf(stdNormal,
                                  -hosts[j]->getParams()[k].param->getMean() /
                                  sd) *
                              (1 - x[sample] /
                               static_cast<double>(lhsSamples + 1)));
              // std::cout << (x[sample] / static_cast<double>(lhsSamples + 1)) << " "
              //           << cdf(stdNormal,
              //                     -hosts[j]->getParams()[k].param->getMean() /
              //                  sd) << " "
              //           << quantile(stdNormal,
              //                 x[sample] / static_cast<double>(lhsSamples + 1) +
              //                 cdf(stdNormal,
              //                     -hosts[j]->getParams()[k].param->getMean() /
              //                     sd) *
              //                 (1 - x[sample] /
              //                  static_cast<double>(lhsSamples + 1))) << " "
              //           << hosts[j]->getParams()[k].param->value << " "
              //           <<std::endl;
            }
            ++sample;
          }
        }
      }
      for (size_t j = 0; j < vectors.size(); ++j) {
        for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
          // std::cout << vectors[j]->getName() << " "
          //           << vectors[j]->getParams()[k].option << " " 
          //           << vectors[j]->getParams()[k].param->getMean() << " "
          //           << vectors[j]->getParams()[k].param->value << " "
          //           << vectors[j]->getParams()[k].param->limits.first << " " 
          //           << vectors[j]->getParams()[k].param->limits.second
          //           << std::endl;
          if (vectors[j]->getParams()[k].param->limits.first >= 0 &&
              vectors[j]->getParams()[k].param->limits.second >= 0) {
            if (vectors[j]->getParams()[k].param->sampling == Linear) {
              vectors[j]->getParams()[k].param->value =
                vectors[j]->getParams()[k].param->limits.first +
                (vectors[j]->getParams()[k].param->limits.second -
                 vectors[j]->getParams()[k].param->limits.first) *
                (x[sample] - 1) / static_cast<double>(lhsSamples - 1);
            } else if (vectors[j]->getParams()[k].param->sampling == Log) {
              vectors[j]->getParams()[k].param->value =
                vectors[j]->getParams()[k].param->limits.first +
                (vectors[j]->getParams()[k].param->limits.second -
                 vectors[j]->getParams()[k].param->limits.first) *
                (exp(10* (x[sample] - 1) / static_cast<double>(lhsSamples - 1)) - 1) /
                (exp(10) - 1);
            } else if (vectors[j]->getParams()[k].param->sampling == Normal) {
              double sd = fabs(vectors[j]->getParams()[k].param->limits.second -
                               vectors[j]->getParams()[k].param->limits.first) /
                2.;
              vectors[j]->getParams()[k].param->value =
                vectors[j]->getParams()[k].param->getMean() +
                sd * quantile(stdNormal,
                              x[sample] / static_cast<double>(lhsSamples + 1) +
                              cdf(stdNormal,
                                  -vectors[j]->getParams()[k].param->getMean() /
                                  sd) *
                              (1 - x[sample] /
                               static_cast<double>(lhsSamples + 1)));
            }
            ++sample;
          }
        }
      }
    }

    // normalise biting
    double bitingSum = .0;
    for (std::vector<Host*>::iterator it = hosts.begin();
         it != hosts.end(); it++) {
      bitingSum += (*it)->f.value;
    }
    for (std::vector<Host*>::iterator it = hosts.begin();
         it != hosts.end(); it++) {
      (*it)->f.value /= bitingSum;
    }
    for (std::vector<Group>::iterator it = groups.begin();
         it != groups.end(); it++) {
      it->f = .0;
      for (std::vector<size_t>::iterator it2 = it->members.begin();
           it2 != it->members.end(); it2++) {
        it->f += hosts[*it2]->f.value;
      }
    }

    betafunc_params p (hosts, vectors, groups, global);
  
    if (samples > 0) {
      // generate beta distributions of prevalence
      for (size_t j = 0; j < hosts.size(); ++j) {
        boost::math::beta_distribution<> distribution
          (hosts[j]->M.value+1, hosts[j]->N.value-hosts[j]->M.value+1);
        p.hPrevalence[j] = quantile(distribution, randGen());
      }
      for (size_t j = 0; j < vectors.size(); ++j) {
        boost::math::beta_distribution<> distribution
          (vectors[j]->M.value+1, vectors[j]->N.value-vectors[j]->M.value+1);
        p.vPrevalence[j] = quantile(distribution, randGen());
      }
    }
    
    std::stringstream outLine;
    if (samples > 0) {
      outLine << i;
      if (addColumn.length() > 0) {
        outLine << ",";
      }
    }
    if (addColumn.length() > 0) {
      outLine << addColumn;
    }
    
    std::vector<double> contrib (hosts.size());
    double domestic;
    double domwild;
    double wildlife;
    double humdom;
    double humwild;
    double R0;

    int status = GSL_SUCCESS;
    
    arma::mat K;
    arma::mat S;
    arma::mat T;

    size_t matrixSize = hosts.size() + groups.size() * vectors.size();
    for (size_t v = 0; v < vectors.size(); ++v) {
      if (vectors[v]->alpha.value > 0) {
        matrixSize += groups.size();
      }
    }
    
    K.set_size(matrixSize, matrixSize);
    S.set_size(matrixSize, matrixSize);
    T.set_size(matrixSize, matrixSize);

    std::vector<double> vars; // beta^*, p^v_i and alpha, the variables
  
    if (calcBetas) { // estimate betas

      // find betas, p^v_is and alpha
      size_t nAttempts = 0;
      do {
        status = betaffoiv(&p, vars, jacobian, verbose);
        ++nAttempts;
      } while (status != GSL_SUCCESS && nAttempts < attempts);

      if (status == GSL_SUCCESS) {

        for (size_t i = 0; i < hosts.size(); ++i) {
          hosts[i]->b.value = vars[i];
        }
        if (global->estimateXi) {
          for (size_t v = 0; v < vectors.size(); ++v) {
            vectors[v]->xi.value = vars[v + hosts.size()];
          }
        } else {
          for (size_t v = 0; v < vectors.size(); ++v) {
            vectors[v]->b.value = vars[v + hosts.size()];
          }
        }
        for (size_t v = 0; v < vectors.size(); ++v) {
          vectors[v]->groupPrev.resize(groups.size(), .0);
          for (size_t g = 0; g < groups.size(); ++g) {
            vectors[v]->groupPrev[g] =
              vars[g+v * groups.size() + hosts.size() + vectors.size()];
          }
        }

        // compose NGM
        S.zeros();
        T.zeros();

        size_t matrixCounter = 0;
        for (size_t v = 0; v < vectors.size(); ++v) {
          size_t vectorStart = matrixCounter;
          for (size_t j = 0; j < groups.size(); ++j) {
            double denominator = .0;
            for (size_t l = 0; l < groups.size(); ++l) {
              denominator += p.habitatOverlap[j][l] *
                groups[l].f;
            }
            S(matrixCounter, matrixCounter) =
              S(matrixCounter, matrixCounter) -
              vectors[v]->mu.value - vectors[v]->xi.value;
            if (vectors[v]->alpha.value > 0) {
              S(matrixCounter, matrixCounter) =
                S(matrixCounter, matrixCounter) -
                vectors[v]->alpha.value;
              S(matrixCounter + 1, matrixCounter) =
                S(matrixCounter + 1, matrixCounter) +
                vectors[v]->alpha.value;
              S(matrixCounter + 1, matrixCounter + 1) =
                S(matrixCounter + 1, matrixCounter + 1) -
                vectors[v]->mu.value - vectors[v]->xi.value;
            }
            size_t matrixCounter2 = vectorStart;
            for (size_t k = 0; k < groups.size(); ++k) {
              S(matrixCounter,matrixCounter2) =
                S(matrixCounter,matrixCounter2) +
                vectors[v]->xi.value * groups[j].f *
                p.habitatOverlap[j][k] / denominator;
              if (vectors[v]->alpha.value > 0) {
                S(matrixCounter + 1, matrixCounter2 + 1) =
                  S(matrixCounter + 1, matrixCounter2 + 1) +
                  vectors[v]->xi.value * groups[j].f *
                  p.habitatOverlap[j][k] / denominator;
                ++matrixCounter2;
              }
              ++matrixCounter2;
            }
            for (size_t k = 0; k < groups[j].members.size(); ++k) {
              size_t i = groups[j].members[k];
              T(matrixCounter, matrixSize - hosts.size() + i) =
                vectors[v]->b.value * vectors[v]->tau.value * hosts[i]->f.value /
                hosts[i]->n.value;
            }
            if (vectors[v]->alpha.value > 0) {
              ++matrixCounter;
            }
            for (size_t k = 0; k < groups[j].members.size(); ++k) {
              size_t i = groups[j].members[k];
              T(matrixSize - hosts.size() + i, matrixCounter) =
                hosts[i]->b.value * vectors[v]->tau.value * hosts[i]->f.value /
                groups[j].f;
            }
            ++matrixCounter;
          }
        }
        for (size_t i = 0; i < hosts.size(); ++i) {
          S(matrixSize - hosts.size() + i,
            matrixSize - hosts.size() + i) =
            S(matrixSize - hosts.size() + i,
              matrixSize - hosts.size() + i) -
            hosts[i]->gamma.value - hosts[i]->mu.value;
        }

        if (verbose >= 2) {
          T.print("T:");
          S.print("S:");
          arma::mat I = -inv(S);
          I.print("inv(S):");
        }

        K = - T * inv(S);

        if (verbose) {
          K.print("K:");
        }
      
        arma::cx_vec eigval;
        arma::cx_mat eigvec;
      
        arma::eig_gen(eigval, eigvec, K);
        arma::colvec colr0 = arma::real(eigval);
        R0 = arma::max(colr0);
        if (verbose) {
          for (size_t i = 0; i < hosts.size() + groups.size(); ++i) {
            std::cout << eigval(i) << std::endl;
          }
        }

        arma::mat P;
        arma::mat KP;
        P.copy_size(K);
        KP.copy_size(K);
        P.zeros();
        for (size_t i = 0; i < matrixSize - hosts.size(); ++i) {
          P(i,i) = 1;
        }
        arma::mat tempP;
        tempP.copy_size(P);

        for (size_t i = 0; i < hosts.size(); ++i) {
          tempP = P;
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
          KP = tempP * K;
          if (verbose >= 2) {
            std::stringstream s;
            s << "KP[" << i << "]";
            KP.print(s.str());
          }
          arma::eig_gen(eigval, eigvec, KP);
          arma::colvec col = arma::real(eigval);
          contrib[i] = arma::max(col);
        }

        tempP = P;
        for (size_t i = 1; i < 4; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec coldomestic = arma::real(eigval);
        domestic = arma::max(coldomestic);
      
        tempP = P;
        for (size_t i = 1; i < 12; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec coldomwild = arma::real(eigval);
        domwild = arma::max(coldomwild);
      
        tempP = P;
        for (size_t i = 4; i < 12; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec colwild = arma::real(eigval);
        wildlife = arma::max(colwild);

        tempP = P;
        for (size_t i = 0; i < 4; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec colhumdom = arma::real(eigval);
        humdom = arma::max(colhumdom);
        
        tempP = P;
        for (size_t i = 0; i < 1; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        for (size_t i = 4; i < 12; ++i) {
          tempP(i + matrixSize - hosts.size(),
                i + matrixSize - hosts.size()) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec colhumwild = arma::real(eigval);
        humwild = arma::max(colhumwild);
      } else {
        std::cerr << "ERROR: Could not find solution" << std::endl;
      }
    } else {

      std::vector<double> hostContrib;

      double hostSum = .0;
      
      for (size_t i = 0; i < hosts.size(); ++i) {
        hostSum += hosts[i]->f.value * p.hPrevalence[i];
      }

      for (size_t i = 0; i < hosts.size(); ++i) {
        hosts[i]->b.value = p.hPrevalence[i] / (1-p.hPrevalence[i]) *
          hosts[i]->n.value *
          (hosts[i]->gamma.value + hosts[i]->mu.value) /
          (vectors[0]->tau.value * hosts[i]->f.value * p.vPrevalence[0]);
      }

      vectors[0]->b.value = 
        p.vPrevalence[0] / (1-p.vPrevalence[0]) *
        vectors[0]->mu.value / vectors[0]->tau.value / hostSum;
      
      double hostContribSum = 0;
      for (size_t i = 0; i < hosts.size(); ++i) {
        double newContrib =
          1 / ((1 - p.hPrevalence[i]) * (1 - p.vPrevalence[0])) *
          p.hPrevalence[i] * hosts[i]->f.value / hostSum;
        hostContrib.push_back(newContrib);
        hostContribSum += newContrib;
      }
      if (verbose) {
        std::cout << "Measured vector prevalence: " 
                  << vectors[0]->M.value/static_cast<double>(vectors[0]->N.value)
                  << std::endl;
      }

      domestic = .0;
      domwild = .0;
      wildlife = .0;
      humdom = .0;
      humwild = .0;
      R0 = .0;

      for (size_t i = 0; i < hosts.size(); ++i) {
        R0 += hostContrib[i];
        if (i == 0) {
          humdom += hostContrib[i];
          humwild += hostContrib[i];
        } else if (i < 4) {
          domestic += hostContrib[i];
          domwild += hostContrib[i];
          humdom += hostContrib[i];
        } else {
          wildlife += hostContrib[i];
          domwild += hostContrib[i];
          humwild += hostContrib[i];
        }
      
        contrib[i] = sqrt(hostContrib[i]);
        if (verbose) {
          std::cout << "hosts[" << i << "].gamma=" << hosts[i]->gamma.value
                    << ", hosts[" << i << "].mu=" << hosts[i]->mu.value
                    << ", hosts[" << i << "].n="
                    << hosts[i]->n.value << std::endl;
        }
      }
      R0 = sqrt(R0);
      domestic = sqrt(domestic);
      domwild = sqrt(domwild);
      wildlife = sqrt(wildlife);
      humdom = sqrt(humdom);
      humwild = sqrt(humwild);
    }
  
    // output

    if (verbose) {
      for (size_t i = 0; i < hosts.size(); ++i) {
        std::cout << "\\hat{b^h_" << i << "}=" << hosts[i]->b.value
                  << std::endl;
      }
      std::cout << "b_hat <-" << std::endl;
      std::cout << "  c(" << hosts[0]->b.value;
      for (size_t i = 1; i < hosts.size(); ++i) {
        std::cout << ", " << hosts[i]->b.value;
      }
      std::cout << ")" << std::endl;
      std::cout << std::endl;
      for (size_t v = 0; v < vectors.size(); ++v) {
        if (global->estimateXi) {
          std::cout << "xi[" << v << "]" << "=" << vectors[v]->xi.value
                    << std::endl;;
        } else {
          std::cout << "b^v_" << v << "=" << vectors[v]->b.value
                    << std::endl;
        }
      }
      if (global->estimateXi) {
        std::cout << "xi <- c(" << vectors[0]->xi.value;
        for (size_t v = 1; v < vectors.size(); ++v) {
          std::cout << ", " << vectors[i]->xi.value;
        }
        std::cout << ")" << std::endl;
      } else {
        std::cout << "bv <- c(" << vectors[0]->b.value;
        for (size_t v = 1; v < vectors.size(); ++v) {
          std::cout << ", " << vectors[v]->b.value;
        }
      }
      std::cout << ")" << std::endl;
      
      for (size_t v = 0; v < vectors.size(); ++v) {
        for (size_t g = 0; g < groups.size(); ++g) {
          std::cout << "p_{" << g << "," << v << "}="
                    << vectors[v]->groupPrev[g]
                    << std::endl;
        }
      }
      std::cout << "pv <- " << std::endl;
      std::cout << "  c(" << vars[hosts.size() + vectors.size()];
      bool first = true;
      for (size_t v = 0; v < vectors.size(); ++v) {
        for (size_t j = 0; j < groups.size(); ++j) {
          if (first) {
            first = false;
          } else {
            std::cout << ", " << ;
          }
        }
      }
      std::cout << ")" << std::endl;
    }
    
    for (size_t j = 0; j < hosts.size(); ++j) {
      outLine << "," << p.hPrevalence[j];
      if (lhsSamples > 0) {
        for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
          outLine << "," << hosts[j]->getParams()[k].param->value;
        }
      }
    }
    
    for (size_t v = 0; v < vectors.size(); ++v) {
      outLine << "," << p.vPrevalence[v];
      if (groups.size() > 1) {
        for (size_t g = 0; g < groups.size(); ++g) {
          outLine << "," << vectors[v]->groupPrev[g];
        }
      }
      if (lhsSamples > 0) {
        for (size_t k = 0; k < vectors[v]->getParams().size(); ++k) {
          outLine << "," << vectors[v]->getParams()[k].param->value;
        }
      }
    }
    
    if (status == GSL_SUCCESS) {
      if (vm.count("print")) {
        std::cout << "\nNGM contributions:" << std::endl;
      }
      
      for (size_t j = 0; j < groups.size(); ++j) {
        for (size_t k = 0; k < groups[j].members.size(); ++k) {
          size_t i = groups[j].members[k];
          if (vm.count("print")) {
            std::cout << hosts[i]->getName() << ": " << contrib[i] << std::endl;
          }
          outLine << "," << contrib[i];
        }
      }
    
      if (vm.count("print")) {
        std::cout << "\nR0: " << R0 << std::endl;
      
        std::cout << "\nDomestic cycle: " << domestic << std::endl;
        std::cout << "Wildlife cycle: " << wildlife << std::endl;
        std::cout << "Domestic+wildlife: " << domwild << std::endl;
        std::cout << "Human+domestic: " << humdom << std::endl;
        std::cout << "Human+wildlife: " << humwild << std::endl;
      }
      outLine << "," << domestic << "," << wildlife << "," << domwild << ","
              << humdom << "," << humwild << "," << R0;
      if (!((outFile == "-" && vm.count("print")) || sim)) {
        out << outLine.str() << std::endl;
      }
    }
      
    ++i;
    if (lhsSamples > 0 && (i % lhsSamples == 0)) {
      free(x);
    }

    if (sim) {
      double t = 0;  // Current time
      double tn = 0; // Next time point at which the state of the system will be
                     // recorded

      size_t nvars = hosts.size();
      // for (size_t v = 0; v < vectors.size(); ++v) {
      size_t vectorTypes = 1 + 
        (vectors[0]->alpha.value > 0 ? 1 : 0) +
        (global->teneralOnly ? 1 : 0);
      nvars += groups.size() * vectorTypes;
      // }
                                
      std::vector<double> data(nvars, 0.);

      size_t index = 0;

      if (!(vm.count("noheader"))) {
        out << "time";
      }
    
      // initialise values    
      for (size_t h = 0; h < hosts.size(); ++h) {
        data[index] = (hosts[h]->x0.value >= 0 ?
                       hosts[h]->x0.value :
                       hosts[h]->M.value);
        ++index;
        if (!(vm.count("noheader"))) {
          out << ",I" << (h+1);
        }
      }

      // for the moment, only work with one vector
      for (size_t g = 0; g < groups.size(); ++g) {
        std::stringstream group_option;
        group_option << "g" << g << "-x0";
        double x0 = vm[group_option.str().c_str()].as<double>();
        if (x0 < 0) {
          x0 = vars[g + hosts.size() + vectors.size()] * groups[g].f *
            vectors[0]->N.value;
        } else {
          for (size_t i = 0; i < groups[g].members.size(); ++i) {
            data[groups[g].members[i]] = x0;
          }
        }

        if (global->teneralOnly) {
          data[index] = vectors[0]->mu.value / (vectors[0]->tau.value + vectors[0]->mu.value) *
            groups[g].f * vectors[0]->N.value;
          ++index;
          if (!(vm.count("noheader"))) {
            out << ",Tv" << (g + 1);
          }
        }

        if (vectors[0]->alpha.value > 0) {
          double inf =
            vectors[0]->alpha.value/(vectors[0]->alpha.value +
                                     vectors[0]->mu.value) *
            ((vectors[0]->alpha.value + vectors[0]->mu.value) /
             (vectors[0]->alpha.value + vectors[0]->mu.value +
              vectors[0]->xi.value) * vars[g + hosts.size() + vectors.size()] +
             (vectors[0]->xi.value) /
             (vectors[0]->alpha.value + vectors[0]->mu.value +
              vectors[0]->xi.value) *
             vectors[0]->M.value / vectors[0]->N.value);

          double incub = ((vectors[0]->mu.value + vectors[0]->xi.value) * inf -
                          vectors[0]->xi.value * vectors[0]->M.value /
                          vectors[0]->N.value) / vectors[0]->alpha.value;

          data[index] = round(x0 * incub / (inf + incub));
          ++index;
          if (!(vm.count("noheader"))) {
            out << ",Cv" << (g + 1);
          }
          data[index] = round(x0 * inf / (inf + incub));
          ++index;
          if (!(vm.count("noheader"))) {
            out << ",Iv" << (g + 1);
          }
        } else {
          data[index] = round(x0);
          ++index;
          if (!(vm.count("noheader"))) {
            out << ",Iv" << (g + 1);
          }
        }
      }
      if (!(vm.count("noheader"))) {
        out << std::endl;
      }
      
      // run gillespie
      do {
        if (t >= tn) {
          out << t;
          for (size_t i = 0; i < data.size(); ++i) {
            out << "," << data[i];
          }
          out << std::endl;
          tn += recordStep;
        }

        double eventRateSum = .0;
        std::vector<Event> events;

        for (size_t h = 0; h < hosts.size(); ++h) {
          Event newEvent;

          // susceptible host infected by vector
          newEvent.add = h;
          newEvent.remove = -1;
          newEvent.rate = vars[h] / hosts[h]->n.value *
            hosts[h]->f.value / groups[hostGroups[h]].f * 
            (vectors[0]->tau.value *
             data[hosts.size() + hostGroups[h] * vectorTypes +
                  vectorTypes - 1] /
             (vectors[0]->N.value)) * (hosts[h]->N.value - data[h]);
          events.push_back(newEvent);
          if (verbose >= 2) {
            std::cout << "Event added: infection; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }
          eventRateSum += newEvent.rate;

          // infected host losing infectiousness
          newEvent.add = -1;
          newEvent.remove = h;
          newEvent.rate = (hosts[h]->mu.value + hosts[h]->gamma.value) * data[h];
          events.push_back(newEvent);
          eventRateSum += newEvent.rate;
          if (verbose >= 2) {
            std::cout << "Event added: infectiousness loss; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }
        }

        for (size_t g = 0; g < groups.size(); ++g) {
          Event newEvent;

          double foi = 0;
          for (size_t i = 0; i < groups[g].members.size(); ++i) {
            size_t h = groups[g].members[i];
            foi += hosts[h]->f.value * data[h] / hosts[h]->N.value;
          }
          foi *= vectors[0]->tau.value / groups[g].f;
          double susc;
          if (global->teneralOnly) {
            susc = data[hosts.size() + g * vectorTypes];
          } else {
            susc = vectors[0]->N.value * groups[g].f;
            for (size_t i = 0; i < vectorTypes; ++i) {
              susc -= data[hosts.size() + g * vectorTypes + i];
            }
          }
          foi *= susc;

          // vector infected by host
          if (global->teneralOnly) {
            newEvent.add = hosts.size() + g * vectorTypes + 1;
            newEvent.remove = hosts.size() + g * vectorTypes;
          } else {
            newEvent.add = hosts.size() + g * vectorTypes;
            newEvent.remove = -1;
          }
          newEvent.rate = vars[hosts.size()] * foi;
          events.push_back(newEvent);
          eventRateSum += newEvent.rate;
          if (verbose >= 2) {
            std::cout << "Event added: vector infection; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }

          // infection maturing in incubating vector
          if (vectors[0]->alpha.value > 0) {
            newEvent.add = hosts.size() + g * vectorTypes + vectorTypes - 1;
            newEvent.remove = hosts.size() + g * vectorTypes + vectorTypes - 2;
            newEvent.rate =
              vectors[0]->alpha.value *
              data[hosts.size() + g * vectorTypes + vectorTypes - 2];
            events.push_back(newEvent);
            eventRateSum += newEvent.rate;
          if (verbose >= 2) {
            std::cout << "Event added: infection maturing; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }
          }
            
          // incubating vector death
          if (vectors[0]->alpha.value > 0) {
            if (global->teneralOnly) {
              newEvent.add = hosts.size() + g * vectorTypes;
            } else {
              newEvent.add = -1;
            }
            newEvent.remove = hosts.size() + g * vectorTypes + vectorTypes - 2;
            newEvent.rate =
              vectors[0]->mu.value *
              data[hosts.size() + g * vectorTypes + vectorTypes - 2];
            events.push_back(newEvent);
            eventRateSum += newEvent.rate;
          if (verbose >= 2) {
            std::cout << "Event added: incubating death; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }
          }

          // infected vector death
          if (global->teneralOnly) {
            newEvent.add = hosts.size() + g * vectorTypes;
          } else {
            newEvent.add = -1;
          }
          newEvent.remove = hosts.size() + g * vectorTypes + vectorTypes - 1;
          newEvent.rate = 
            vectors[0]->mu.value *
            data[hosts.size() + g * vectorTypes + vectorTypes - 1];
          events.push_back(newEvent);
          eventRateSum += newEvent.rate;
          if (verbose >= 2) {
            std::cout << "Event added: infected death; add " << newEvent.add << ", remove "
                      << newEvent.remove << ", rate " << newEvent.rate << std::endl;
          }

          // non-infected vector death
          if (global->teneralOnly) {
            newEvent.add = hosts.size() + g * vectorTypes;
            newEvent.remove = -1;
            newEvent.rate = vectors[0]->N.value;
            for (size_t i = 1; i < vectorTypes; ++i) {
              newEvent.rate -= data[hosts.size() + g * vectorTypes + i];
            }
            newEvent.rate *= vectors[0]->mu.value;
            events.push_back(newEvent);
            eventRateSum += newEvent.rate;
            if (verbose >= 2) {
              std::cout << "Event added: noninfectious death; add " << newEvent.add << ", remove "
                        << newEvent.remove << ", rate " << newEvent.rate << std::endl;
            }

            // non-infectious bite on host
            newEvent.add = -1;
            newEvent.remove = hosts.size() + g * vectorTypes;
            newEvent.rate = (1 - vars[hosts.size()]) * susc;
            events.push_back(newEvent);
            eventRateSum += newEvent.rate;
            if (verbose >= 2) {
              std::cout << "Event added: noninfectious bite; add " << newEvent.add << ", remove "
                        << newEvent.remove << ", rate " << newEvent.rate << std::endl;
            }
          }

          // host switches
          for (size_t k = 0; k < groups.size(); ++k) {
            if (k != g) {
              if (global->teneralOnly) {
                newEvent.add = hosts.size() + k * vectorTypes;
                newEvent.remove = hosts.size() + g * vectorTypes;
                newEvent.rate =
                  vectors[0]->xi.value * groups[k].f *
                  data[hosts.size() + g * vectorTypes];
                events.push_back(newEvent);
                eventRateSum += newEvent.rate;
                if (verbose >= 2) {
                  std::cout << "Event added: host switch; add " << newEvent.add << ", remove "
                            << newEvent.remove << ", rate " << newEvent.rate << std::endl;
                }
              }

              if (vectors[0]->alpha.value > 0) {
                newEvent.add = hosts.size() + k * vectorTypes + vectorTypes - 2;
                newEvent.remove = hosts.size() + g * vectorTypes + vectorTypes - 2;
                newEvent.rate =
                  vectors[0]->xi.value * groups[k].f *
                  data[hosts.size() + g * vectorTypes + vectorTypes - 2];
                events.push_back(newEvent);
                eventRateSum += newEvent.rate;
                if (verbose >= 2) {
                  std::cout << "Event added: host switch; add " << newEvent.add << ", remove "
                            << newEvent.remove << ", rate " << newEvent.rate << std::endl;
                }
              }
                
              newEvent.add = hosts.size() + k * vectorTypes + vectorTypes - 1;
              newEvent.remove = hosts.size() + g * vectorTypes + vectorTypes - 1;
              newEvent.rate =
                vectors[0]->xi.value * groups[k].f *
                data[hosts.size() + g * vectorTypes + vectorTypes - 1];
              events.push_back(newEvent);
              eventRateSum += newEvent.rate;
              if (verbose >= 2) {
                std::cout << "Event added: host switch; add " << newEvent.add << ", remove "
                          << newEvent.remove << ", rate " << newEvent.rate << std::endl;
              }
            }
          }
        }

        if (eventRateSum > 0) {

          // determine timestep and update time
          double r = -log(randGen())/eventRateSum;
          t += r;

          // determine event
          r = randGen() * eventRateSum;
          if (verbose >= 2) {
            std::cout << "random: " << r << " out of " << eventRateSum
                      << std::endl;
          }

          size_t e = 0;
          do {
            r -= events[e].rate;
            if (verbose >= 2) {
              std::cout << "Event: " << e << " step " << r << std::endl;
            }
            if (r <= 0) {
              if (verbose >= 2) {
                std::cout << "chosen" << std::endl;
              }
              if (events[e].add >= 0) {
                ++data[events[e].add];
                if (verbose >= 2) {
                  std::cout << "adding a " << events[e].add << std::endl;
                }
              }
              if (events[e].remove >= 0) {
                --data[events[e].remove];
                if (verbose >= 2) {
                  std::cout << "removing a " << events[e].remove << std::endl;
                }
              }
            }
            ++e;
          } while (r > 0);
        } else {
          t = INFINITY;
        }
      } while (t < tmax && (humanmax == 0 || data[0] < humanmax));
      out << t;
      for (size_t i = 0; i < data.size(); ++i) {
        out << "," << data[i];
      }
      out << std::endl;
    }
  } while (i < samples);
  
  of.close();

  for (std::vector<Host*>::iterator it = hosts.begin();
       it != hosts.end(); it++) {
    delete *it;
  }
  for (std::vector<Vector*>::iterator it = vectors.begin();
       it != vectors.end(); it++) {
    delete *it;
  }
  if (global) delete global;
}
   
//------------------------------------------------------------
