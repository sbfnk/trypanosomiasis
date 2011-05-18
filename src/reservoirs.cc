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

namespace po = boost::program_options;

/*! Main program. */
int main(int argc, char* argv[])
{

  std::string dataFile; //!< name of the file holding host parameters
  std::string vectorFile; //!< name of the file holding vector parameters
  std::string outFile; //!< output file

  bool gambiense = false; //!< include gambiense data?
  bool nonGambiense = false; //!< include nongambiense data?
  bool jacobian = false; //!< use jacobian?
  size_t lhsSamples = 0; //!< number of samples to be used in linear hypercube
                         //!< sampling 
  bool calcBetas = true; //!< calculate the betas (or use analytical formula in
                          //!< simple model?)

  size_t samples = 0; //!< number of samples to be used from negative binomial
                      //!< distribution around measured prevalence

  double xi; //!< xi (assortativity)
  double eta; //!< eta (habitat preference)

  unsigned int verbose = 0; //!< be verbose

  std::vector<group> groups; //!< composition of groups

  std::string firstColumn = "0";
  std::string firstHeader = "n";
  size_t attempts = 1;

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
    ("species,s", po::value<std::string>()->default_value("g"),
     "trypanosome species to consider (g=gambiense, n=nongambiense)")
    ("groups,g", po::value<std::string>(),
     "groups (semicolon-separated list of comma-separated hosts")
    ("xi,x", po::value<double>()->default_value(.0),
     "assortativity parameter xi")
    ("eta,a", po::value<double>()->default_value(.0),
     "habitat preference xi")
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

  // one of the two trypanosome species must be specified
  gambiense = (vm["species"].as<std::string>().find_first_of("g") !=
               std::string::npos);
  nonGambiense = (vm["species"].as<std::string>().find_first_of("n") !=
               std::string::npos);
  if (!(gambiense || nonGambiense)) {
    std::cerr << "Error: must include some trypanosome species" << std::endl;
    return 1;
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

  xi = vm["xi"].as<double>();

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

  std::vector<host> hosts; // vector of hosts
  std::vector<vector> vectors; // vector of vectors

  // tokenizer for reading csv file
  typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
  boost::escaped_list_separator<char> sep('\\', ',', '\"');

  std::vector<std::string> headings;
  bool firstLine = true;

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
      host newHost(lineVector, headings, gambiense, nonGambiense);
      hosts.push_back(newHost);
    }
  }
  in.close();

  po::options_description host_options
    ("\nHost options");

  for (size_t i = 0; i < hosts.size(); ++i) {
    host_options.add_options()
      ((hosts[i].name+"-M").c_str(), po::value<size_t>());
    host_options.add_options()
      ((hosts[i].name+"-N").c_str(), po::value<size_t>());
    host_options.add_options()
      ((hosts[i].name+"-mu").c_str(), po::value<double>());
    host_options.add_options()
      ((hosts[i].name+"-gamma").c_str(), po::value<double>());
    host_options.add_options()
      ((hosts[i].name+"-theta").c_str(), po::value<double>());
    host_options.add_options()
      ((hosts[i].name+"-abundance").c_str(), po::value<double>());
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
      vector newVector(lineVector, headings, gambiense, nonGambiense);
      vectors.push_back(newVector);
    }
  }
  in.close();

  po::options_description vector_options
    ("\nVector options");
  
  for (size_t i = 0; i < vectors.size(); ++i) {
    vector_options.add_options()
      ((vectors[i].name+"-M").c_str(), po::value<size_t>());
    vector_options.add_options()
      ((vectors[i].name+"-N").c_str(), po::value<size_t>());
    vector_options.add_options()
      ((vectors[i].name+"-mu").c_str(), po::value<double>());
    vector_options.add_options()
      ((vectors[i].name+"-gamma").c_str(), po::value<double>());
    vector_options.add_options()
      ((vectors[i].name+"-biting-rate").c_str(), po::value<double>());
    vector_options.add_options()
      ((vectors[i].name+"-density").c_str(), po::value<double>());
  }

  // ************** read additional parameters ****************

  po::options_description long_options;
  long_options.add(main_options).
    add(host_options).add(vector_options);

  if (vm.count("longhelp")) {
    std::cout << long_options << std::endl;
    return 0;
  }
  
  try {
    po::store(po::command_line_parser(argc, argv).
              options(long_options).run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);

  for (size_t i = 0; i < hosts.size(); ++i) {
    if (vm.count((hosts[i].name+"-M").c_str())) {
      hosts[i].M = vm[(hosts[i].name+"-M").c_str()].as<size_t>();
    }
    if (vm.count((hosts[i].name+"-N").c_str())) {
      hosts[i].N = vm[(hosts[i].name+"-N").c_str()].as<size_t>(); 
    }
    if (vm.count((hosts[i].name+"-mu").c_str())) {
      hosts[i].mu = vm[(hosts[i].name+"-mu").c_str()].as<double>(); 
    }
    if (vm.count((hosts[i].name+"-gamma").c_str())) {
      hosts[i].gamma = vm[(hosts[i].name+"-gamma").c_str()].as<double>(); 
    }
    if (vm.count((hosts[i].name+"-theta").c_str())) {
      hosts[i].theta = vm[(hosts[i].name+"-theta").c_str()].as<double>(); 
    }
    if (vm.count((hosts[i].name+"-abundance").c_str())) {
      hosts[i].abundance = vm[(hosts[i].name+"-abundance").c_str()].as<double>();
    }
  }
  for (size_t i = 0; i < vectors.size(); ++i) {
    if (vm.count((vectors[i].name+"-M").c_str())) {
      vectors[i].M = vm[(vectors[i].name+"-M").c_str()].as<size_t>(); 
    }
    if (vm.count((vectors[i].name+"-N").c_str())) {
      vectors[i].N = vm[(vectors[i].name+"-N").c_str()].as<size_t>(); 
    }
    if (vm.count((vectors[i].name+"-mu").c_str())) {
      vectors[i].mu = vm[(vectors[i].name+"-mu").c_str()].as<double>(); 
    }
    if (vm.count((vectors[i].name+"-gamma").c_str())) {
      vectors[i].gamma = vm[(vectors[i].name+"-gamma").c_str()].as<double>(); 
    }
    if (vm.count((vectors[i].name+"-biting-rate").c_str())) {
      vectors[i].bitingRate = vm[(vectors[i].name+"-biting-rate").c_str()].as<double>();
    }
    if (vm.count((vectors[i].name+"-density").c_str())) {
      vectors[i].density = vm[(vectors[i].name+"-density").c_str()].as<double>(); 
    }
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
      group newGroup;
      for (std::vector<std::string>::iterator it2 = subres.begin();
           it2 != subres.end(); it2++) {
        size_t s;
        std::istringstream myStream(*it2);
        myStream >> s;
        newGroup.members.push_back(s);
        newGroup.theta += hosts[s].theta;
      }
      groups.push_back(newGroup);
    }
  } else {
    if (vm.count("random")) { // random mixing
      groups.push_back(group());
      groups[0].theta = 1;
      for (size_t i = 0; i < hosts.size(); ++i) {
        groups[0].members.push_back(i);
      }
    } else { // species-specific mixing
      for (size_t i = 0; i < hosts.size(); ++i) {
        groups.push_back(group(i, hosts[i].theta));
      }
    }
  }

  // create a copy of the vector variable for each group
  for (size_t i = 1; i < groups.size(); ++i) {
    vector newVector(vectors[0]);
    vectors.push_back(newVector);
  }
  
  // ********************* estimate betas *********************
  
  arma::mat K(hosts.size() + groups.size(), hosts.size() + groups.size());
  arma::mat S(hosts.size() + groups.size(), hosts.size() + groups.size());
  arma::mat T(hosts.size() + groups.size(), hosts.size() + groups.size());
  
  betafunc_params p (hosts, vectors, groups, xi);

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
       (hosts[j].M+1, hosts[j].N-hosts[j].M+1));
  }
  for (size_t j = 0; j < vectors.size(); ++j) {
    distributions.push_back
      (new boost::math::beta_distribution<>
       (vectors[j].M+1, vectors[j].N-vectors[j].M+1));
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
      out << ",\"" << hosts[j].name << "_prev\"";
      if (lhsSamples > 0) {
        out << ",\"" << hosts[j].name << "_abundance\"";
        out << ",\"" << hosts[j].name << "_mu\"";
        out << ",\"" << hosts[j].name << "_gamma\"";
      }
    }
    for (size_t j = 0; j < 1; ++j) {
      out << ",\"" << vectors[0].name << "_prev\"";
    }
    for (size_t j = 0; j < hosts.size(); ++j) {
      out << ",\"" << hosts[j].name << "\"";
    }
    out << ",\"domestic\",\"wildlife\",\"domestic+wildlife\",\"R0\"";
    out << std::endl;
  }

  size_t i = 0;
  std::stringstream outLine;
  do {

    if (lhsSamples > 0) {
      if (i % lhsSamples == 0) {
        // generate latin hypercube samples
        x = (int*) malloc (3 * hosts.size() * lhsSamples * sizeof(int));
        ihs(3*hosts.size(), lhsSamples, 5, &seed, x);
      }
      
      for (size_t j = 0; j < hosts.size(); ++j) {
        hosts[j].abundance = hosts[j].abundance_limits.first +
          (hosts[j].abundance_limits.second -
           hosts[j].abundance_limits.first) *
          x[i % lhsSamples * hosts.size() * 3 + j] /
          static_cast<double>(lhsSamples);
        hosts[j].mu = hosts[j].mu_limits.first +
          (hosts[j].mu_limits.second -
           hosts[j].mu_limits.first) *
          x[i % lhsSamples * hosts.size() * 3 + j + 1] /
          static_cast<double>(lhsSamples);
        hosts[j].gamma = hosts[j].gamma_limits.first +
          (hosts[j].gamma_limits.second -
           hosts[j].gamma_limits.first) *
          x[i % lhsSamples * hosts.size() * 3 + j + 2] /
          static_cast<double>(lhsSamples);
      }
    }

    if (samples > 0) {
      for (size_t j = 0; j < hosts.size(); ++j) {
        p.hPrevalence[j] = quantile(*distributions[j], randGen());
      }
      // to be corrected for more vectors
      // for (size_t j = 0; j < vectors.size(); ++j) {
      p.vPrevalence = quantile(*distributions[hosts.size()],
                               randGen());
    }

    if (samples == 0) {
      outLine << firstColumn;
    } else {
      outLine << i;
    }
    
    for (size_t j = 0; j < hosts.size(); ++j) {
      outLine << "," << p.hPrevalence[j];
      if (lhsSamples > 0) {
        outLine << "," << hosts[j].abundance;
        outLine << "," << hosts[j].mu;
        outLine << "," << hosts[j].gamma;
      }
    }
    outLine << "," << p.vPrevalence;

    std::vector<double> contrib (hosts.size());
    double domestic;
    double domwild;
    double wildlife;
    double R0;

    int status = GSL_SUCCESS;
    
    if (calcBetas) { // estimate betas
      std::vector<double> vars; // beta^*, p^v_i and alpha, the variables

      // find betas, p^v_is and alpha
      size_t nAttempts = 0;
      do {
        status = betaffoiv(&p, vars, jacobian, verbose);
        ++nAttempts;
      } while (status != GSL_SUCCESS && nAttempts < attempts);

      if (status == GSL_SUCCESS) {
        if (verbose) {
          for (size_t i = 0; i < hosts.size(); ++i) {
            std::cout << "beta^*[" << i << "]=" << vars[i] << std::endl;
          }
          for (size_t i = 0; i < groups.size(); ++i) {
            std::cout << "p^v_i[" << i << "]=" << vars[i+hosts.size()]
                      << std::endl;
          }
          for (size_t i = 0; i < 1; ++i) {
            std::cout << "alpha[" << i << "]=" << vars[i + hosts.size() +
                                                       groups.size()];
          }
          std::cout << std::endl;
        }
      
        // compose NGM
        S.zeros();
        T.zeros();
      
        for (size_t j = 0; j < groups.size(); ++j) {
          S(j,j) = S(j,j) - vectors[0].mu - xi;
          for (size_t k = 0; k < groups.size(); ++k) {
            S(j,k) = S(j,k) + xi * groups[k].theta;
          }
          for (size_t k = 0; k < groups[j].members.size(); ++k) {
            size_t i = groups[j].members[k];
            T(j,i + groups.size()) = vars[hosts.size() + groups.size()] *
              vars[i] * hosts[i].theta / hosts[i].abundance;
            T(i + groups.size(),j) = vars[i] * hosts[i].theta / groups[j].theta;
          }
        }
        for (size_t i = 0; i < hosts.size(); ++i) {
          S(i + groups.size(),i + groups.size()) =
            S(i + groups.size(),i + groups.size()) -
            hosts[i].gamma - hosts[i].mu;
        }
      
        K = - T * inv(S);

        if (verbose) {
          T.print("T:");
          S.print("S:");
          arma::mat I = inv(S);
          I.print("inv(S):");
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

        arma::mat P(hosts.size() + groups.size(), hosts.size() + groups.size());
        arma::mat KP(hosts.size() + groups.size(), hosts.size() + groups.size());
        P.zeros();
        for (size_t i = 0; i < groups.size(); ++i) {
          P(i,i) = 1;
        }
        arma::mat tempP(hosts.size() + groups.size(),
                        hosts.size() + groups.size());

        for (size_t i = groups.size(); i < groups.size() + hosts.size(); ++i) {
          tempP = P;
          tempP(i,i) = 1;
          KP = tempP * K;
          arma::eig_gen(eigval, eigvec, KP);
          arma::colvec col = arma::real(eigval);
          contrib[i-groups.size()] = arma::max(col);
        }

        tempP = P;
        for (size_t i = groups.size() + 1; i < groups.size() + 4; ++i) {
          tempP(i,i) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec coldomestic = arma::real(eigval);
        domestic = arma::max(coldomestic);
      
        tempP = P;
        for (size_t i = groups.size() + 1; i < groups.size() + 12; ++i) {
          tempP(i,i) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec coldomwild = arma::real(eigval);
        domwild = arma::max(coldomwild);
      
        tempP = P;
        for (size_t i = groups.size() + 4; i < groups.size() + 12; ++i) {
          tempP(i,i) = 1;
        }
        KP = tempP * K;
        arma::eig_gen(eigval, eigvec, KP);
        arma::colvec colwild = arma::real(eigval);
        wildlife = arma::max(colwild);
      } else {
        std::cerr << "ERROR: Could not find solution" << std::endl;
      } 
    } else {

      std::vector<double> K;
      std::vector<double> hostContrib;
      double hostContribSum = 0;
      for (size_t i = 0; i < hosts.size(); ++i) {
        double newContrib =
          pow(p.hPrevalence[i], 2) / (1 - p.hPrevalence[i]) *
          hosts[i].abundance * (hosts[i].gamma + hosts[i].mu);
        if (verbose) {
          std::cout << "Host contrib " << i << ": " << newContrib << std::endl;
        }
        hostContrib.push_back(newContrib);
        hostContribSum += newContrib;
      }
      if (verbose) {
        std::cout << "Measured vector prevalence: " 
                  << vectors[0].M/vectors[0].N << std::endl;
      }
      for (size_t i = 0; i < hosts.size(); ++i) {
        K.push_back(1/((1-p.vPrevalence)*(1-p.hPrevalence[i])) *
                    hostContrib[i] / hostContribSum);
      }

      domestic = .0;
      domwild = .0;
      wildlife = .0;
      R0 = .0;

      for (size_t i = 0; i < hosts.size(); ++i) {
        R0 += K[i];
        if (i == 0) {
        } else if (i < 4) {
          domestic += K[i];
          domwild += K[i];
        } else {
          wildlife += K[i];
          domwild += K[i];
        }
      
        contrib[i] = sqrt(K[i]);
        if (verbose) {
          std::cout << "hosts[" << i << "].gamma=" << hosts[i].gamma
                    << ", hosts[" << i << "].mu=" << hosts[i].mu
                    << ", vectors[0].density=" << vectors[0].density
                    << ", hosts[" << i << "].abundance="
                    << hosts[i].abundance << std::endl;
        }
      }
      R0 = sqrt(R0);
      domestic = sqrt(domestic);
      domwild = sqrt(domwild);
      wildlife = sqrt(wildlife);
    }

    if (status == GSL_SUCCESS) {
      if (vm.count("print")) {
        std::cout << "\nNGM contributions:" << std::endl;
      }
      
      for (size_t j = 0; j < groups.size(); ++j) {
        for (size_t k = 0; k < groups[j].members.size(); ++k) {
          size_t i = groups[j].members[k];
          if (vm.count("print")) {
            std::cout << hosts[i].name << ": " << contrib[i] << std::endl;
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
}
   
//------------------------------------------------------------
