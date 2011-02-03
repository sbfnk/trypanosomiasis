#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <sys/time.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "reservoirs.hh"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{

  std::string dataFile;
  std::string vectorFile;
  std::string paramsFile;
  std::string outFile;

  bool gambiense = false;
  bool nonGambiense = false;
  bool vectorPrevalence = false;
  bool bitingPreference = false;

  size_t samples = 0;

  unsigned int verbose = 0;

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
    ("params-file,p", po::value<std::string>()->
     default_value("params/params.csv"),
     "params file")
    ("output-file,o", po::value<std::string>()->default_value("reservoirs.dat"),
     "output file")
    ("species,s", po::value<std::string>()->default_value("g"),
     "trypanosome species to consider (g=gambiense, n=nongambiense)")
    ("convert-area,c", 
     "factor to convert area to prevalence")
    ("vector-prevalence,r", 
     "consider vector prevalence")
    ("biting-preference,b", 
     "consider biting preference")
    ("samples,m", po::value<size_t>()->default_value(0),
     "number of samples")
    ;

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

  if (vm.count("help")) {
    std::cout << main_options << std::endl;
    return 0;
  }

  if (vm.count("verbose")) {
    verbose = 1;
  }
  if (vm.count("very-verbose")) {
    verbose = 2;
  }
  
  if (vm.count("data-file")) {
    dataFile = vm["data-file"].as<std::string>();
  } else {
    std::cerr << "Error: must specify data file" << std::endl;
    return 1;
  }

  if (vm.count("vector-file")) {
    vectorFile = vm["vector-file"].as<std::string>();
  } else {
    std::cerr << "Error: must specify vector file" << std::endl;
    return 1;
  }

  if (vm.count("params-file")) {
    paramsFile = vm["params-file"].as<std::string>();
  } else {
    std::cerr << "Error: must specify params file" << std::endl;
    return 1;
  }

  samples = vm["samples"].as<size_t>();

  if (vm.count("output-file")) {
    outFile = vm["output-file"].as<std::string>();
  } else {
    if (samples > 0) {
      std::cerr << "Error: must specify output file for samples >0"
                << std::endl;
      return 1;
    }
  }

  gambiense = (vm["species"].as<std::string>().find_first_of("g") !=
               std::string::npos);
  nonGambiense = (vm["species"].as<std::string>().find_first_of("n") !=
               std::string::npos);
  if (!(gambiense || nonGambiense)) {
    std::cerr << "Error: must include some trypanosome species" << std::endl;
    return 1;
  }

  if (vm.count("vector-prevalence")) {
    vectorPrevalence = true;
  }

  if (vm.count("biting-preference")) {
    bitingPreference = true;
  }

  std::vector<host> hosts;
  std::vector<vector> vectors;
  param params;

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

  // ********************** read params file ********************
  
  in.open(paramsFile.c_str());
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
      params = param(lineVector, headings);
    }
  }
  in.close();

  po::options_description param_options
    ("\nParameter options");

  param_options.add_options()
    ("area-convert", po::value<double>());

  if (!vm.count("convert-area")) {
    params.areaConvert = 1.;
  }

  if (vm.count("longhelp")) {
    std::cout << main_options << host_options << vector_options
              << param_options;
  }

  // ************** read additional parameters ****************

  po::options_description long_options;
  long_options.add(main_options).
    add(host_options).add(vector_options).
    add(param_options);
  
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
  
  if (vm.count("area-convert")) {
    params.areaConvert = vm["area-convert"].as<double>();
  }
  
  // ********************* estimate betas *********************

  betafunc_params p (hosts, vectors, params, vectorPrevalence,
                     bitingPreference);
  std::vector<double> beta;
  
  if (samples == 0) {

    betaffoiv(&p, beta);

    double r0 = 0;

    double domestic = 0;
    double wildLife = 0;
    double animal = 0;
  
    std::cout << "\nNGM contributions: \n";
    for (size_t i = 0; i < hosts.size(); ++i) {
      double K = pow(beta[i] * hosts[i].theta * vectors[0].bitingRate, 2) * 
        params.areaConvert * vectors[0].density / hosts[i].abundance /
        ((hosts[i].gamma + hosts[i].mu) * vectors[0].mu);
      r0 += K;
      if (i == 0) {
      } else if (i < 4) {
        domestic += K;
        animal += K;
      } else {
        wildLife += K;
        animal += K;
      }
      
      double contrib = sqrt(K);
      // std::cout << "beta[" << i << "]=" << beta[i] << ", hosts[" << i
      //           << "].gamma=" << hosts[i].gamma << ", hosts[" << i << "].mu="
      //           << hosts[i].mu << ", params.areaConvert="
      //           << params.areaConvert << ", vectors[0]="
      //           << vectors[0].density << ", hosts[" << i << "].abundance="
      //           << hosts[i].abundance << std::endl;
      std::cout << hosts[i].name << ": " << contrib << std::endl;

    }
    std::cout << "\nR0: " << sqrt(r0) << std::endl;

    std::cout << "\nDomestic cycle: " << sqrt(domestic) << std::endl;
    std::cout << "Wildlife cycle: " << sqrt(wildLife) << std::endl;
    std::cout << "Domestic+wildlife: " << sqrt(animal) << std::endl;
  } else {
    std::ofstream out(outFile.c_str());
    std::vector<boost::math::beta_distribution<>*> distributions;

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

    unsigned int seed;
    struct timeval tv;
    gettimeofday(&tv, 0);
    seed = tv.tv_sec + tv.tv_usec;
    boost::mt19937 gen(seed);
    
    boost::variate_generator<boost::mt19937, boost::uniform_real<> >
      randGen(gen, boost::uniform_real<> (0,1));

    for (size_t i = 0; i < samples; ++i) {
      out << i;
      for (size_t j = 0; j < hosts.size(); ++j) {
        p.hPrevalence[j] = quantile(*distributions[j], randGen());
        out << " " << p.hPrevalence[j];
      }
      if (p.useVectorPrevalence) {
        for (size_t j = 0; j < vectors.size(); ++j) {
          p.vPrevalence[j] = quantile(*distributions[j+hosts.size()],
                                      randGen());
          out << " " << p.vPrevalence[j];
        }
      }
        

      betaffoiv(&p, beta);
      for (size_t j = 0; j < beta.size(); ++j) {
        out << " " << beta[j];
      }
      out << std::endl;
    }
    
    out.close();
  }
}
   
//------------------------------------------------------------
