#ifndef ODE_IO_UTILS_H
#define ODE_IO_UTILS_H

//------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <string>
#include <exception>
#include <boost/program_options.hpp>

//#include "ode_solver.hh"

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

std::string find_ode_type(int argc, char* argv[])
{
  std::string type("no type specified");
  
  int i;
  for (i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"--ode-type=mf")) {
      type = "mf";
      break;
    } else if (!strcmp(argv[i],"--ode-type=pa")) {
      type = "pa";
      break;
    }
  }
  
  return type;
}

//------------------------------------------------------------

po::options_description* generate_main_options()
{
  po::options_description* opt =
    new po::options_description("Usage: ./ode_solver.x -o ode_params_file -m model_params_file\n                      -g graph_params_file [options]...\n\nMain options");

  opt->add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("ode-type",po::value<std::string>(),
     "the specific model to solve [mf/pa]")
    ("ode-params-file,o",po::value<std::string>(),
     "file containing ode parameters")
    ("model-params-file,m",po::value<std::string>(),
     "file containing model parameters")
    ("graph-params-file,g",po::value<std::string>(),
     "file containing graph parameters")
    ("no-gp-file",
     "do not write gnuplot script params file")
    ("check-convergence",
     "check for convergence and stop before Tmax")
    ;
  
  return opt;
}

//------------------------------------------------------------

po::options_description* generate_ode_options()
{
  po::options_description* opt =
    new po::options_description("ODE parameters");
  
  opt->add_options()
    ("file-id", po::value<std::string>(),
     "file id for output files")
    ("ic-file", po::value<std::string>(),
     "initial conditions file")
    ("tmax", po::value<double>(),
     "stopping time")
    ("dt", po::value<double>(),
     "size of initial time step")
    ("nsave", po::value<unsigned int>(),
     "save solution every nsave time steps")
    ("step-algo", po::value<std::string>(),
     "name of stepping algorithm type")
    ("abs-tol", po::value<double>(),
     "absolute tolerance")
    ("rel-tol", po::value<double>(),
     "relative tolerance")
    ;

  return opt;
}

//------------------------------------------------------------

bool parse_comm_line_args(int argc, char* argv[],
                          po::options_description& all_options,
                          po::variables_map& vm)
{
  
  std::vector<std::string> unregistered;
  
  try {
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(all_options).allow_unregistered().run();
    po::store(parsed, vm);
    unregistered =
      po::collect_unrecognized(parsed.options, po::exclude_positional);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line arguments: " << e.what()
              << std::endl;
    return 1;
  }
  
  po::notify(vm);
  
  for (std::vector<std::string>::iterator it = unregistered.begin();
       it != unregistered.end(); it++) {
    std::cerr << "WARNING: ignoring unknown option " << *it << std::endl;
  }

  return 0;
  
}

//------------------------------------------------------------

bool parse_ode_args_file(po::options_description& all_ode_options,
                         po::variables_map& vm)
{

  if (vm.count("ode-params-file")) {
    std::ifstream ifs(vm["ode-params-file"].as<std::string>().c_str());
    try {
      po::store(po::parse_config_file(ifs, all_ode_options), vm);
    }
    catch (std::exception& e) {
      std::cerr << "Error parsing params file: " << e.what()
                << std::endl;
      return 1;
    }
  } else {
    std::cerr << "WARNING: no ode-params-file given" << std::endl
              << "Reading all parameters from command line" << std::endl
              << "Hopefully, no parameters are missing ..." << std::endl
              << std::endl;
  }
  
  po::notify(vm); 

  return 0;
  
}

//------------------------------------------------------------

bool parse_model_args_file(po::options_description& model_options,
                           po::variables_map& vm)
{
  
  if (vm.count("model-params-file")) {
    std::ifstream ifs(vm["model-params-file"].as<std::string>().c_str());
    try {
      po::store(po::parse_config_file(ifs, model_options), vm);
    }
    catch (std::exception& e) {
      std::cerr << "Error parsing params file: " << e.what()
                << std::endl;
      return 1;
    }
  } else {
    std::cerr << "WARNING: no model-params-file given" << std::endl
              << "Reading all parameters from command line" << std::endl
              << "Hopefully, no parameters are missing ..." << std::endl
              << std::endl;
  }
    
  po::notify(vm); 

  return 0;
  
}

//------------------------------------------------------------

bool parse_graph_args_file(po::options_description& graph_options,
                           po::variables_map& vm)
{

  if (vm.count("graph-params-file")) {
    std::ifstream ifs(vm["graph-params-file"].as<std::string>().c_str());
    try {
      po::store(po::parse_config_file(ifs, graph_options), vm);
    }
    catch (std::exception& e) {
      std::cerr << "Error parsing params file: " << e.what()
                << std::endl;
      return 1;
    }
  } else {
    std::cerr << "WARNING: no graph-params-file given" << std::endl
              << "Reading all parameters from command line" << std::endl
              << "Hopefully, no parameters are missing ..." << std::endl
              << std::endl;
  }
    
  po::notify(vm); 
  
  return 0;
  
}

//------------------------------------------------------------

template <class Params, class Eqs>
int init_ode_params(po::variables_map& vm,
                    ode::OdeSolver<Params, Eqs>& x)
{
  if (vm.count("file-id")) {
    // ode file id name 
    x.set_file_id(vm["file-id"].as<std::string>());      
  } else  {
    std::cerr << "ERROR: no file-id given" << std::endl;
    return 1;
  }
  if (vm.count("ic-file")) {
    // ode ic file name 
    x.set_ic_file_name(vm["ic-file"].as<std::string>().c_str());      
  } else  {
    std::cerr << "ERROR: no ic-file given" << std::endl;
    return 1;
  }   
  if (vm.count("tmax")) {
    x.set_tmax(vm["tmax"].as<double>());
  } else {
    std::cerr << "ERROR: no tmax given" << std::endl;
    return 1;
  }
  if (vm.count("dt")) {
    x.set_dt(vm["dt"].as<double>());
  } else {
    std::cerr << "ERROR: no dt given" << std::endl;
    return 1;
  }
  if (vm.count("nsave")) {
    x.set_nsave(vm["nsave"].as<unsigned int>());
  } else {
    std::cerr << "ERROR: no nsave given" << std::endl;
    return 1;
  }
  if (vm.count("step-algo")) {
    x.set_step_algo(vm["step-algo"].as<std::string>().c_str());
  } else {
    std::cerr << "ERROR: no step-algo given" << std::endl;
    return 1;
  }
  if (vm.count("abs-tol")) {
    x.set_abs_tol(vm["abs-tol"].as<double>());
  } else {
    std::cerr << "ERROR: no abs-tol given" << std::endl;
    return 1;
  }
  if (vm.count("rel-tol")) {
    x.set_rel_tol(vm["rel-tol"].as<double>());
  } else {
    std::cerr << "ERROR: no rel-tol given" << std::endl;
    return 1;
  }
  if (vm.count("verbose")) {
    x.set_verbose(true);
  }
  if (vm.count("check-convergence")) {
    x.set_convergence_check(true);
  }
  
  return 0;
}

//------------------------------------------------------------
  
int write_gp_script(std::string fname, double Qd, double Qi,
                    double tmax, double N, bool verbose)
{
  // file name
  std::ofstream gpFile;
  std::string gpFileName = fname + ".gp";
  if (verbose) std::cout << "... writing gnuplot output file" << std::endl;

  // open file
  try {
    gpFile.open(gpFileName.c_str(), std::ios::out);
  }
  
  catch (std::exception &e) {
    std::cerr << "... unable to open gnuplot output file " 
              << gpFileName << std::endl;
    std::cerr << "... Standard exception: " << e.what() << std::endl;      
    return 1; 
  }

  // writing to file
  gpFile << "### model parameters generated by ode_solver" << std::endl;
  gpFile << "N=" << N << std::endl;
  gpFile << "Qd=" << Qd << std::endl;
  gpFile << "Qi=" << Qi << std::endl;
  gpFile << "Tmax=" << tmax << std::endl;
  gpFile << "### end of model parameters" << std::endl << std::endl;
  
  // closing file
  gpFile.close();
  
  if (verbose) std::cout << " ... done\n";
  return 0;   
}

//------------------------------------------------------------

#endif // ODE_IO_UTILS_H
