#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include "reservoirs.hh"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{

  std::string dataFile();
  std::string vectorFile();
  std::string paramsFile();

  bool gambiense = false;
  bool nonGambiense = false;
  bool vectorPrevalence = false;

  param params;

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
    ("params-file,p", po::value<std::string>()->
     default_value("params/params.csv"),
     "params file")
    ("species,s", po::value<std::string>()->
     default_value("g")
     "trypanosome species to consider (g=gambiense, n=nongambiense)")
    ("area-convert,a", 
     "factor to convert area to prevalence")
    ("vector-prevalence,r", 
     "consider vector prevalence")
    ;

  po::variables_map vm;
    
  try {
    po::store(po::command_line_parser(argc, argv).options(main_options).
              positional(file_option).run(), vm); 
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);

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

  // ********************** read vector file ********************
  
  in = std::ifstream(vectorFile.c_str());
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

  // ********************** read params file ********************
  
  in = std::ifstream(vectorFile.c_str());
  firstLine = true;
  
  std::ifstream in(vectorFile.c_str());
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

  if (!vm.count("area-convert")) {
    params.areaConvert = 1.;
  }
  
  // ********************* estimate betas *********************

  betafunc_params p(hosts, vectors, params, vectorPrevalence);

  std::vector<double> beta;

  betaffoiv(p, beta);

  std::cout << "\nNGM contributions: \n";
  for (size_t i = 0; i < hosts.size(); ++i) {
    double K = beta[i] * beta[i] / (hosts[i].gamma + hosts[i].mu) *
      area_convert * vector[0].density / hosts[i].abundance;
    double contrib = sqrt(K*K);
    std::cout << hosts[i].name << ": " << K << std::endl;
  }
}
   
//------------------------------------------------------------
