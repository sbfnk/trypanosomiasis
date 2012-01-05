#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "ode_solver.hh"
#include "tryps.hh"
#include "ode_io_utils.hh"

//------------------------------------------------------------

template <class T>
int solve_ode_system(int argc, char* argv[]);

//------------------------------------------------------------

int main(int argc, char* argv[])
{
  int status;
  
  status = solve_ode_system<ode::OdeSolver<tryps::Params, tryps::Eqs> >(argc, argv);
  
  return status;
  
}

//------------------------------------------------------------

template <class T>
int solve_ode_system(int argc, char* argv[])
{
  // ode system variable
  T ode_sys;
  
  // generate main options
  po::options_description* main_options = generate_main_options();
  
  // generate ode options
  po::options_description* ode_options = generate_ode_options();
  
  // generate specific model options
  po::options_description* model_options = generate_model_options(ode_sys);
  
  // gathering options
  po::options_description all_options, all_ode_options;
  
  all_ode_options.add(*main_options).add(*ode_options);
  all_options.add(*main_options).add(*ode_options).add(*model_options);
  
  po::options_description help_options, long_help_options;
  help_options.add(*main_options);
  long_help_options.add(*main_options).add(*ode_options).add(*model_options);
  
  //------------------------------------------------------------
  // parsing command line and files
  
  // variables map for command line args
  po::variables_map vm;
  
  // parse command line arguments
  parse_comm_line_args(argc, argv, all_options, vm);
  
  // parse ode-params file
  parse_ode_args_file(all_ode_options, vm);
  
  // parse model-params file
  parse_model_args_file(*model_options, vm);
  
  //------------------------------------------------------------
  // help messages
  
  // print help messages
  if (vm.count("help"))
    {
      std::cout << help_options << std::endl;
      return 1;
    }
  if (vm.count("longhelp"))
    {
      std::cout << long_help_options << std::endl;
      return 1;
    }

  //------------------------------------------------------------
  
  // initializing ode parameters
  init_ode_params(vm, ode_sys);

  std::ofstream cmdLineFile;
  std::string cmdLineFileName = ode_sys.get_file_id() + ".cmd_line";
  cmdLineFile.open(cmdLineFileName.c_str(), std::ios::out);
  cmdLineFile << argv[0] << " ";
  for (int j = 1; j < argc; j++) {
    cmdLineFile << argv[j] << " ";
  }
  cmdLineFile << std::endl;
  cmdLineFile.close();

  // initializing model parameters
  init_model_params(vm, ode_sys);
  
  // print ode_system
  std::cout << ode_sys << std::endl;
  
  // init rhs
  ode_sys.init_rhs();
  
  // solve
  ode_sys.solve();
  
  //------------------------------------------------------------
  // write parameters to Gnuplot script
  
  // cleaning
  delete main_options;
  delete ode_options;
  delete model_options;

  return 0;
}

//------------------------------------------------------------
