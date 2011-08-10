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

  std::string firstColumn = "0"; //!< entry of first column in data output --
                                 //!this can be controlled as a command line
                                 //!parameter if samples = 0 (only one line of
                                 //!data output)
  std::string firstHeader = "n"; //!< header of first column (normally n for
                                 //!sample number)
  size_t attempts = 1; //!< number of attempts for solving equations

  // main options
  po::options_description main_options
    ("Usage: reservoirs [options]... \n\nOptions");

  main_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message")
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
     "number of samples for latin hypercube sampling")
    ("samples,n", po::value<size_t>()->default_value(0),
     "number of samples (0 for no sampling)")
    ("jacobian,j", 
     "use jacobian")
    ("simple,i", 
     "use simple model")
    ("random,r", 
     "assume random mixing")
    ("noheader,d", 
     "do not print header")
    ("print,p", 
     "print results")
    ("firstcolumn,c", po::value<std::string>(),
     "entry of first column")
    ("firstheader,C", po::value<std::string>(),
     "header of first column")
    ("attempts,t", po::value<size_t>()->default_value(1),
     "number of solving attempts")
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

  if (vm.count("longhelp")) {
    Host *host = new Host();
    Vector *vector = new Vector();
    GlobalParams *global = new GlobalParams();
    std::cout << main_options << *(global->getOptions())
              << *(host->getOptions()) << *(vector->getOptions()) << std::endl;
    delete host;
    delete vector;
    delete global;
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

  if (vm.count("simple")) {
    calcBetas = false;
  }

  if (vm.count("firstcolumn")) {
    if (samples == 0) {
      firstColumn = vm["firstcolumn"].as<std::string>();
    } else {
      std::cerr << "WARNING: --firstcolumn ignored because samples > 0"
                << std::endl;
    }
  }

  if (vm.count("firstheader")) {
    if (samples == 0) {
      firstHeader = vm["firstheader"].as<std::string>();
    } else {
      std::cerr << "WARNING: --firstheader ignored because samples > 0"
                << std::endl;
    }
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

  // ************** read group composition option ****************
    
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
    if (vm.count("random")) { // random mixing
      groups.push_back(Group());
      for (size_t i = 0; i < hosts.size(); ++i) {
        groups[0].members.push_back(i);
      }
    } else { // species-specific mixing
      for (size_t i = 0; i < hosts.size(); ++i) {
        groups.push_back(Group(i));
      }
    }
  }

  // ********************* estimate betas *********************

  betafunc_params p (hosts, vectors, groups, global);
  
  std::streambuf * buf;
  std::ofstream of;
  if (outFile == "-") {
    buf = std::cout.rdbuf();
  } else {
    of.open(outFile.c_str());
    buf = of.rdbuf();
  }
  
  std::ostream out(buf);

  int seed = get_seed();
  int *x = 0;
  
  std::vector<boost::math::beta_distribution<>*> distributions;
  
  // generate beta distributions of prevalence
  for (size_t j = 0; j < hosts.size(); ++j) {
    distributions.push_back
      (new boost::math::beta_distribution<>
       (hosts[j]->M.first+1, hosts[j]->N.first-hosts[j]->M.first+1));
  }
  for (size_t j = 0; j < vectors.size(); ++j) {
    distributions.push_back
      (new boost::math::beta_distribution<>
       (vectors[j]->M.first+1, vectors[j]->N.first-vectors[j]->M.first+1));
  }
  
  boost::mt19937 gen(seed);
  
  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    randGen(gen, boost::uniform_real<> (0,1));

  if (!(vm.count("noheader") || (outFile == "-" && vm.count("print")))) {
    if (samples == 0) {
      out << firstHeader;
    } else {
      out << "n";
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
    for (size_t j = 0; j < vectors.size(); ++j) {
      out << ",\"" << vectors[j]->getName() << "_prev\"";
      if (lhsSamples > 0) {
        for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
          out << ",\"" << vectors[j]->getName() << "_"
              << vectors[j]->getParams()[k].option << "\"";
        }
      }
    }
    for (size_t j = 0; j < hosts.size(); ++j) {
      out << ",\"" << hosts[j]->getName() << "\"";
    }
    out << ",\"domestic\",\"wildlife\",\"domestic+wildlife\",\"R0\"";
    out << std::endl;
  }

  size_t i = 0;
  do {

    size_t nParams = 0;
    for (size_t j = 0; j < hosts.size(); ++j) {
      for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
        if (hosts[j]->getParams()[k].param->second.first >= 0 &&
            hosts[j]->getParams()[k].param->second.second >= 0) {
          ++nParams;
        }
      }
    }
    for (size_t j = 0; j < vectors.size(); ++j) {
      for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
        if (vectors[j]->getParams()[k].param->second.first >= 0 &&
            vectors[j]->getParams()[k].param->second.second >= 0) {
          ++nParams;
        }
      }
    }
    
    if (lhsSamples > 0) {
      if (i % lhsSamples == 0) {
        // generate latin hypercube samples
        x = (int*) malloc (nParams * lhsSamples * sizeof(int));
        ihs(nParams, lhsSamples, 5, &seed, x);
      }
      
      size_t sample = 0;
      for (size_t j = 0; j < hosts.size(); ++j) {
        for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
          if (hosts[j]->getParams()[k].param->second.first >= 0 &&
              hosts[j]->getParams()[k].param->second.second >= 0) {
            hosts[j]->getParams()[k].param->first =
              hosts[j]->getParams()[k].param->second.first +
              (hosts[j]->getParams()[k].param->second.second -
               hosts[j]->getParams()[k].param->second.first) *
              x[sample] /
              static_cast<double>(lhsSamples);
            ++sample;
          }
        }
      }
      for (size_t j = 0; j < vectors.size(); ++j) {
        for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
          if (vectors[j]->getParams()[k].param->second.first >= 0 &&
              vectors[j]->getParams()[k].param->second.second >= 0) {
            vectors[j]->getParams()[k].param->first =
              vectors[j]->getParams()[k].param->second.first +
              (vectors[j]->getParams()[k].param->second.second -
               vectors[j]->getParams()[k].param->second.first) *
              x[sample] /
              static_cast<double>(lhsSamples);
            ++sample;
          }
        }
      }
    }

    if (samples > 0) {
      for (size_t j = 0; j < hosts.size(); ++j) {
        p.hPrevalence[j] = quantile(*distributions[j], randGen());
      }
      for (size_t j = 0; j < vectors.size(); ++j) {
        p.vPrevalence[j] = quantile(*distributions[j + hosts.size()],
                                    randGen());
      }
    }

    // normalise biting
    double bitingSum = .0;
    for (std::vector<Host*>::iterator it = hosts.begin();
         it != hosts.end(); it++) {
      bitingSum += (*it)->f.first;
    }
    for (std::vector<Host*>::iterator it = hosts.begin();
         it != hosts.end(); it++) {
      (*it)->f.first /= bitingSum;
    }
    for (std::vector<Group>::iterator it = groups.begin();
         it != groups.end(); it++) {
      it->f = .0;
      for (std::vector<size_t>::iterator it2 = it->members.begin();
           it2 != it->members.end(); it2++) {
        it->f += hosts[*it2]->f.first;
      }
    }

    std::stringstream outLine;
    if (samples == 0) {
      outLine << firstColumn;
    } else {
      outLine << i;
    }
    
    std::vector<double> contrib (hosts.size());
    double domestic;
    double domwild;
    double wildlife;
    double R0;

    int status = GSL_SUCCESS;
    
    arma::mat K;
    arma::mat S;
    arma::mat T;

    size_t matrixSize = hosts.size() + groups.size() * vectors.size();
    for (size_t v = 0; v < vectors.size(); ++v) {
      if (vectors[v]->alpha.first > 0) {
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
          hosts[i]->b.first = vars[i];
        }
        if (global->estimateXi) {
          for (size_t v = 0; v < vectors.size(); ++v) {
            vectors[v]->xi.first = vars[v + hosts.size()];
          }
        } else {
          for (size_t v = 0; v < vectors.size(); ++v) {
            vectors[v]->b.first = vars[v + hosts.size()];
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
              vectors[v]->mu.first - vectors[v]->xi.first;
            if (vectors[v]->alpha.first > 0) {
              S(matrixCounter, matrixCounter) =
                S(matrixCounter, matrixCounter) -
                vectors[v]->alpha.first;
              S(matrixCounter + 1, matrixCounter) =
                S(matrixCounter + 1, matrixCounter) +
                vectors[v]->alpha.first;
              S(matrixCounter + 1, matrixCounter + 1) =
                S(matrixCounter + 1, matrixCounter + 1) -
                vectors[v]->mu.first - vectors[v]->xi.first;
            }
            size_t matrixCounter2 = vectorStart;
            for (size_t k = 0; k < groups.size(); ++k) {
              S(matrixCounter,matrixCounter2) =
                S(matrixCounter,matrixCounter2) +
                vectors[v]->xi.first * groups[j].f *
                p.habitatOverlap[j][k] / denominator;
              if (vectors[v]->alpha.first > 0) {
                S(matrixCounter + 1, matrixCounter2 + 1) =
                  S(matrixCounter + 1, matrixCounter2 + 1) +
                  vectors[v]->xi.first * groups[j].f *
                  p.habitatOverlap[j][k] / denominator;
                ++matrixCounter2;
              }
              ++matrixCounter2;
            }
            for (size_t k = 0; k < groups[j].members.size(); ++k) {
              size_t i = groups[j].members[k];
              T(matrixCounter, matrixSize - hosts.size() + i) =
                vectors[v]->b.first * vectors[v]->tau.first * hosts[i]->f.first /
                hosts[i]->n.first;
            }
            if (vectors[v]->alpha.first > 0) {
              ++matrixCounter;
            }
            for (size_t k = 0; k < groups[j].members.size(); ++k) {
              size_t i = groups[j].members[k];
              T(matrixSize - hosts.size() + i, matrixCounter) =
                hosts[i]->b.first * vectors[v]->tau.first * hosts[i]->f.first /
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
            hosts[i]->gamma.first - hosts[i]->mu.first;
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
      } else {
        std::cerr << "ERROR: Could not find solution" << std::endl;
      }
    } else {

      std::vector<double> hostContrib;
      vars.resize(hosts.size() + groups.size() + 1);

      double hostSum = .0;
      
      for (size_t i = 0; i < hosts.size(); ++i) {
        hostSum += hosts[i]->f.first * p.hPrevalence[i];
      }
      
      for (size_t i = 0; i < hosts.size(); ++i) {
        vars[i] = p.hPrevalence[i] / (1-p.hPrevalence[i]) * hosts[i]->n.first *
          (hosts[i]->gamma.first + hosts[i]->mu.first) /
          (vectors[0]->tau.first * hosts[i]->f.first * p.vPrevalence[0]);
      }

      vars[hosts.size() + groups.size()] = p.vPrevalence[0];

      vars[hosts.size()] =
        p.vPrevalence[0] / (1-p.vPrevalence[0]) *
        vectors[0]->mu.first / vectors[0]->tau.first / hostSum;
      
      double hostContribSum = 0;
      for (size_t i = 0; i < hosts.size(); ++i) {
        double newContrib =
          1 / ((1 - p.hPrevalence[i]) * (1 - p.vPrevalence[0])) *
          p.hPrevalence[i] * hosts[i]->f.first / hostSum;
        hostContrib.push_back(newContrib);
        hostContribSum += newContrib;
      }
      if (verbose) {
        std::cout << "Measured vector prevalence: " 
                  << vectors[0]->M.first/static_cast<double>(vectors[0]->N.first)
                  << std::endl;
      }

      domestic = .0;
      domwild = .0;
      wildlife = .0;
      R0 = .0;

      for (size_t i = 0; i < hosts.size(); ++i) {
        R0 += hostContrib[i];
        if (i == 0) {
        } else if (i < 4) {
          domestic += hostContrib[i];
          domwild += hostContrib[i];
        } else {
          wildlife += hostContrib[i];
          domwild += hostContrib[i];
        }
      
        contrib[i] = sqrt(hostContrib[i]);
        if (verbose) {
          std::cout << "hosts[" << i << "].gamma=" << hosts[i]->gamma.first
                    << ", hosts[" << i << "].mu=" << hosts[i]->mu.first
                    << ", hosts[" << i << "].n="
                    << hosts[i]->n.first << std::endl;
        }
      }
      R0 = sqrt(R0);
      domestic = sqrt(domestic);
      domwild = sqrt(domwild);
      wildlife = sqrt(wildlife);
    }
  
    // output

    if (verbose) {
      for (size_t i = 0; i < hosts.size(); ++i) {
        std::cout << "\\hat{b^h_" << i << "}=" << vars[i] << std::endl;
      }
      std::cout << std::endl;
      for (size_t v = 0; v < vectors.size(); ++v) {
        if (global->estimateXi) {
          std::cout << "xi[" << v << "]";
        } else {
          std::cout << "b^v_" << v;
        }
        std::cout << "=" << vars[v + hosts.size()] << std::endl;
      }
      for (size_t v = 0; v < vectors.size(); ++v) {
        for (size_t j = 0; j < groups.size(); ++j) {
          std::cout << "p_{" << j << "," << v << "}="
                    << vars[j+v * groups.size() + hosts.size() + vectors.size()]
                    << std::endl;
        }
      }
    }
    
    for (size_t j = 0; j < hosts.size(); ++j) {
      outLine << "," << p.hPrevalence[j];
      if (lhsSamples > 0) {
        for (size_t k = 0; k < hosts[j]->getParams().size(); ++k) {
          outLine << "," << hosts[j]->getParams()[k].param->first;
        }
      }
    }
    
    for (size_t j = 0; j < vectors.size(); ++j) {
      outLine << "," << p.vPrevalence[j];
      if (lhsSamples > 0) {
        for (size_t k = 0; k < vectors[j]->getParams().size(); ++k) {
          outLine << "," << vectors[j]->getParams()[k].param->first;
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
      }
      outLine << "," << domestic << "," << wildlife << "," << domwild << "," << R0;
      if (!(outFile == "-" && vm.count("print"))) {
        out << outLine.str() << std::endl;
      }
    }
      
    ++i;
    if (lhsSamples > 0 && (i % lhsSamples == 0)) {
      free(x);
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
