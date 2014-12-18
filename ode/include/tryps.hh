#ifndef INFO_SIRS_MF_H
#define INFO_SIRS_MF_H

//------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <gsl/gsl_errno.h>

#include <boost/program_options.hpp>

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

namespace tryps
{

  //------------------------------------------------------------
  // Params structure
  
  struct Params
  {
    Params() : nvars(13) {};
    
    // No. of equations
    unsigned int nvars;
    
    // parameters
    double b[12], f[12], mu[12], gamma[12], n[12];
    double bv, tau, muv;
    double scaling;
    std::string hname[12];
    std::string vname;
    
    // overloading operator<<
    friend std::ostream& operator <<
      (std::ostream& os, const Params& x);
    
  }; // Params
  
  //------------------------------------------------------------
  // overloading operator<<
  
  std::ostream& operator<< (std::ostream& os, const Params& x)
  {
    os << std::endl
       << "Model parameters:\n"
       << "=================\n";
    for (unsigned int i = 0; i < 12; ++i) {
      os << "  " << x.hname[i] << std::endl
         << "    b     = " << x.b[i] << std::endl
         << "    f     = " << x.f[i] << std::endl
         << "    mu    = " << x.mu[i] << std::endl
         << "    gamma = " << x.gamma[i] << std::endl
         << "    n     = " << x.n[i] << std::endl;
    }

    os << "  " << x.vname << std::endl
       << "    bv    = " << x.bv << std::endl
       << "    tau   = " << x.tau << std::endl
       << "    mu    = " << x.muv << std::endl;
    
    return os;
    
  } // operator<<
  
  //------------------------------------------------------------
  // Equations structure
  
  struct Eqs
  {
    // rhs function
    static int rhs_eval (double t, const double y[], double rhs[], void* params)
    {         
      Params p = *(static_cast<Params*>(params));
      
      // local readable short variables

      for (unsigned int i = 0; i < 12; ++i) {
        double lambda = p.b[i] * p.tau * p.f[i] / p.n[i] * y[12];
        rhs[i] = p.scaling * lambda * (1 - y[i]) - (p.mu[i] + p.gamma[i]) * y[i];
      }
      double vlambda = 0;
      for (unsigned int i = 0; i < 12; ++i) {
        vlambda += p.f[i] * y[i];
      }
      vlambda *= p.bv * p.tau;
      rhs[12] = vlambda * (1 - y[12]) - (p.muv) * y[12];
       
      return GSL_SUCCESS;         
    }
    
  }; // Eqs

  //------------------------------------------------------------
  // generating model options
  
  po::options_description* generate_model_options(ode::OdeSolver<Params, Eqs>& dummy)
  {
    // dummy is for overloading
    
    po::options_description* opt =
      new po::options_description("Model parameters\n");

    for (unsigned int i = 0; i < 12; ++i) {
      std::stringstream s;
      s << (i+1);
      opt->add_options()
        ((("hname")+s.str()).c_str(), po::value<std::string>(),
         "host species name")
        (("b"+s.str()).c_str(), po::value<double>(),
         "host-specific susceptibility")
        (("f"+s.str()).c_str(), po::value<double>(),
         "biting preference")
        (("mu"+s.str()).c_str(), po::value<double>(),
         "mortality/birth rate")
        (("gamma"+s.str()).c_str(), po::value<double>(),
         "recovery rate")
        (("n"+s.str()).c_str(), po::value<double>(),
         "abundance")
      ;
    }
    
    opt->add_options()
      ("vname", po::value<std::string>(),
       "vector species name")
      ("bv", po::value<double>(),
       "vector-specific susceptibility")
      ("tau", po::value<double>(),
       "biting rate")
      ("muv", po::value<double>(),
       "vector mortality/birth rate")
      ("scaling", po::value<double>(),
       "scaling ")
      ;
    return opt;
    
  }
  
  //------------------------------------------------------------  
  //  initial model parameters from vm
  
  int init_model_params(po::variables_map& vm,
                        ode::OdeSolver<Params, Eqs>& x)
  {
    Params* model_params = x.get_model_params(); 
    
    for (unsigned int i = 0; i < 12; ++i) {
      std::stringstream s;
      s << (i+1);
      if (vm.count(std::string("hname")+s.str())) {
        model_params->hname[i]=vm[("hname"+s.str()).c_str()].as<std::string>();
      } else {
        std::cerr << "WARNING: no name given to species " << i << std::endl;
        std::cerr << "setting to species" << i  << std::endl;
        model_params->hname[i]=(("species"+s.str()).c_str());
      }
      if (vm.count(("b"+s.str()).c_str())) {
        model_params->b[i]=vm[("b"+s.str()).c_str()].as<double>();
      } else {
        std::cerr << "WARNING: no b given for species " << i << std::endl;
        std::cerr << "setting to 0" << std::endl;
        model_params->b[i]=0;
      }
      if (vm.count(("f"+s.str()).c_str())) {
        model_params->f[i]=vm[("f"+s.str()).c_str()].as<double>();
      } else {
        std::cerr << "WARNING: no f given for species " << i << std::endl;
        std::cerr << "setting to 0" << std::endl;
        model_params->f[i]=0;
      }
      if (vm.count(("mu"+s.str()).c_str())) {
        model_params->mu[i]=vm[("mu"+s.str()).c_str()].as<double>();
      } else {
        std::cerr << "WARNING: no mu given for species " << i << std::endl;
        std::cerr << "setting to 0" << std::endl;
        model_params->mu[i]=0;
      }
      if (vm.count(("gamma"+s.str()).c_str())) {
        model_params->gamma[i]=vm[("gamma"+s.str()).c_str()].as<double>();
      } else {
        std::cerr << "WARNING: no gamma given for species " << i << std::endl;
        std::cerr << "setting to 0" << std::endl;
        model_params->gamma[i]=0;
      }
      if (vm.count(("n"+s.str()).c_str())) {
        model_params->n[i]=vm[("n"+s.str()).c_str()].as<double>();
      } else {
        std::cerr << "WARNING: no n given for species " << i << std::endl;
        std::cerr << "setting to 0" << std::endl;
        model_params->n[i]=0;
      }
    }
    
    if (vm.count("vname")) {
      model_params->vname=vm["vname"].as<std::string>();
    } else {
      std::cerr << "WARNING: no name given to vector species" << std::endl;
      std::cerr << "setting to vspecies" << std::endl;
      model_params->vname="vspecies";
    }
    if (vm.count("bv")) {
      model_params->bv=vm["bv"].as<double>();
    } else {
      std::cerr << "WARNING: no bv given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->bv=0;
    }
    if (vm.count("muv")) {
      model_params->muv=vm["muv"].as<double>();
    } else {
      std::cerr << "WARNING: no muv given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->muv=0;
    }
    if (vm.count("tau")) {
      model_params->tau=vm["tau"].as<double>();
    } else {
      std::cerr << "WARNING: no tau given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->tau=0;
    }
    if (vm.count("scaling")) {
      model_params->scaling=vm["scaling"].as<double>();
    } else {
      std::cerr << "WARNING: no scaling given" << std::endl;
      std::cerr << "setting to 1" << std::endl;
      model_params->scaling=1;
    }

    double R0 = 0;
    for (size_t i = 0; i < 12; ++i) {
      R0 +=
        (model_params->scaling * model_params->bv * model_params->b[i] *
         model_params->tau * model_params->tau * model_params->f[i] *
         model_params->f[i]) /
        (model_params->muv * (model_params->mu[i] + model_params->gamma[i]) *
         model_params->n[i]);
    }

    R0 = sqrt(R0);
    std::cout << R0 << std::endl;
    
    return 0;
  }
  
} // namespace InfoSIRSmf

//------------------------------------------------------------

#endif // INFO_SIRS_MF_H
