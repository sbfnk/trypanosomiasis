#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

//------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <exception>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "convergence.hh"

//------------------------------------------------------------

namespace cv = convergence;

//------------------------------------------------------------

namespace ode
{
  //------------------------------------------------------------
  // declarations (needed for template friend function
  // see C++ FAQ lite [35.16]:
  // http://www.parashift.com/c++-faq-lite/templates.html
  
  template <class Params, class Eqs>
  class OdeSolver;

  template <class Params, class Eqs>
  std::ostream& operator<<
    (std::ostream& os, const OdeSolver<Params, Eqs>& x);
  
  //------------------------------------------------------------
  
  template <class Params, class Eqs>
  class OdeSolver
  {
  public:
    
    // constructors and destructors 
    OdeSolver();
    ~OdeSolver();
    
    // mutators 
    void set_nsave(const size_t nsave) { OdeSolver::nsave = nsave; }
    void set_step_algo(std::string step_algo)
    { OdeSolver::step_algo = step_algo; }
    void set_abs_tol(const double abs_tol) { OdeSolver::abs_tol = abs_tol; }
    void set_rel_tol(const double rel_tol) { OdeSolver::rel_tol = rel_tol; } 
    void set_tmax(const double tmax) { OdeSolver::tmax = tmax; }
    void set_dt(const double dt) { OdeSolver::dt = dt; }
    void set_convergence_check(const bool b)
    { OdeSolver::check_if_converged = b; }
    void set_file_id(std::string file_id) { OdeSolver::file_id = file_id; }
    void set_ic_file_name(std::string ic_file) { OdeSolver::ic_file = ic_file; }
    void set_rhs(double *rhs_ic) { OdeSolver::rhs = rhs_ic; }
    void set_verbose(const bool verbose)
    { OdeSolver::verbose = verbose; }
    
    // accessors 
    size_t get_nsave() const { return nsave; }
    std::string get_step_algo() const { return step_algo; }
    double get_abs_tol() const { return abs_tol; }
    double get_rel_tol() const { return rel_tol; }
    double get_tmax() const { return tmax; }
    double get_dt() const { return dt; }
    bool get_convergence_check() const { return check_if_converged; }
    double* get_rhs() const { return rhs; }
    std::string get_file_id() const {return file_id; }
    std::string get_ic_file_name() const {return ic_file; }
    bool get_verbose() const { return verbose; }
    Params* get_model_params() const { return model_params; }
    Eqs get_model_eqs() const { return model_eqs; }
    
    // gsl/odeiv solve stuff 
    void init_rhs();    
    void solve();
    
    // overloding operators
    friend std::ostream& operator<< <Params, Eqs>
    (std::ostream& os, const OdeSolver<Params, Eqs>& x);
    
  private:
    double *rhs;                 // rhs[navrs] for the variables
    
    // gsl/odeiv solve stuff
    const gsl_odeiv_step_type* set_step_type();    
    //    void init_Q();
    void write_rhs(std::ofstream& output_file, double t);
    void write_last_line();
         
    // ode parameters 
    size_t nsave;          // save solution every nsave time steps
    std::string step_algo; // name of stepping algorithm type
    double abs_tol, rel_tol;     // absolute/relative tolerances
    double tmax, dt;             // ...
    bool verbose;                // print verbose output
    bool check_if_converged;     // convergence check
    
    // output file
    std::string file_id; // file identifier
    std::string ic_file; // file identifier
    
    // Model object
    Params* model_params;
    Eqs model_eqs;
         
  }; // class OdeSolver 

  //------------------------------------------------------------
   
  // constructors and destructors 
   
  template <class Params, class Eqs>
  OdeSolver<Params, Eqs>::OdeSolver() : nsave(1), abs_tol(1e-6),
                                        rel_tol(1e-6), tmax(10), dt(1e-6),
                                        verbose(0),
                                        check_if_converged(false)
  {
    step_algo = "rkf45";
    file_id = "no_file_id";
    ic_file = "";
    rhs=NULL;
    model_params = new Params;
    
  } // OdeSolver 
   
  //------------------------------------------------------------
   
  template <class Params, class Eqs>
  OdeSolver<Params, Eqs>::~OdeSolver()
  {
    delete model_params;
    delete [] rhs;
    
  } // ~OdeSolver 
  
  //------------------------------------------------------------
  
   template <class Params, class Eqs>
   void OdeSolver<Params, Eqs>::init_rhs()
   {
     double *rhs_ic;
     unsigned int nv = model_params->nvars;     
     
     // opening ic_file
     std::string fname;
     if (ic_file == "") {
       fname = file_id + ".init";
     } else {
       fname = ic_file;
     }
     std::ifstream ic;
     try {
       ic.open(fname.c_str(), std::ios::in);
     }   
     catch (std::exception &e) {
       std::cerr << "... unable to open init file " << std::endl
                 << "... Standard exception: " << e.what() << std::endl;
       std::exit(1); 
     }
     
     // allocating rhs 
     if (verbose) std::cout << "... allocating rhs memory";      
     try {      
       rhs_ic = new double[nv]; // will be deleted in destructor
     }      
     catch (std::exception &e) {
       std::cerr << "... unable to alloc rhs\n"
                 << "... Standard exception: " << e.what() << std::endl;      
       std::exit(1); 
     }      
     if (verbose) std::cout << " ... done\n";
     
     // init rhs from ic_file 
     if (verbose) std::cout << "... reading content of ic-file";
     
     for(unsigned int i = 0; i < nv; i++)
       ic >> rhs_ic[i];
     
     set_rhs(rhs_ic);
     
     if (verbose) std::cout << " ... done\n";
     
     // closing ic_file 
     ic.close();
     
   } // init_rhs
  
  //------------------------------------------------------------
  
  template <class Params, class Eqs>
  void OdeSolver<Params, Eqs>::solve()
  {
    // some variables
    double t;
    int status;
    int nv = model_params->nvars;
    
    // convergence check object
    cv::ConvergenceCheck conv(nv);
    
    // opening output file
    std::string fname = file_id + ".dat";
    std::ofstream output_file;
    
    if (nsave > 0) {
      try {
        output_file.open(fname.c_str(), std::ios::out);
      }   
      catch (std::exception &e) {
        std::cerr << "... unable to open output file " << std::endl
                  << "... Standard exception: " << e.what() << std::endl;
        std::exit(1); 
      }
    }
    
    // setting step type 
    const gsl_odeiv_step_type *step_type = set_step_type();
    
    // init additional paramteres if needed - Qi, Qd, Qb      
    //     model_params->init_Q(rhs, verbose);
    
    // printing
    //if (verbose)
    //  std::cout << *this ;
    
    // allocating stepping function 
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type, nv);
    
    // allocating control function 
    gsl_odeiv_control *control = gsl_odeiv_control_y_new(abs_tol,rel_tol);
//     gsl_odeiv_control *control = gsl_odeiv_control_yp_new(abs_tol,rel_tol);
//     gsl_odeiv_control *control = gsl_odeiv_control_standard_new(abs_tol,rel_tol,10,1000);
    
    // allocating evolution function 
    gsl_odeiv_evolve *evolve  = gsl_odeiv_evolve_alloc(nv);
    
    // defining the system 
    gsl_odeiv_system sys =
      {&(model_eqs.rhs_eval), NULL, nv, static_cast<void*>(model_params) };
    
    // write t=0 rhs
    t=0;
    if (nsave > 0)
      write_rhs(output_file, t);
    
    size_t o_count = 1;
    
    if (verbose) 
      std::cout << std::endl << "Initial conditions at t = 0:\n"
                <<              "============================\n"
                << "\033[00;32mS-\033[0m" << " = " << rhs[0] << std::endl
                << "\033[00;31mI-\033[0m" << " = " << rhs[1] << std::endl
                << "\033[00;34mR-\033[0m" << " = " << rhs[2] << std::endl
                << "\033[01;32mS+\033[0m" << " = " << rhs[3] << std::endl
                << "\033[01;31mI+\033[0m" << " = " << rhs[4] << std::endl
                << "\033[01;34mR+\033[0m" << " = " << rhs[5] << std::endl
                << std::endl;
    
    
    //--- main loop ---//
    if (verbose) std::cout << "... doing main loop !!!\n";   
    
    double t_prev = 0;
    bool not_converged = true;
    while (t < tmax)
      {
        // stepping solution
        status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tmax,
                                         &dt, rhs);      
        
        if (status != GSL_SUCCESS)
          exit(1);
        
        // writing RHS to ofile
        if (nsave > 0)
          if (!(o_count % nsave)) 
            write_rhs(output_file, t);
               
        // check convergence
        if (check_if_converged)
          if ((t - t_prev) > conv.get_samples_interval()) {
            if (conv.check(rhs)) {
              std::cout << "... Ode_Solver converged at t = "
                        << t << " !!!" << std::endl;
              not_converged = false;
              break;              
            }
            
            t_prev = t;
          }
        
        
        // advance counter
        ++o_count;        
        
      }

    if (check_if_converged)
      if (not_converged)
        std::cout << "... Ode_Solver did not converge at t = tmax = "
                  << tmax << std::endl;
    
    if (verbose) std::cout << "... writing last line ";
    write_last_line();
    if (verbose) std::cout << "... done\n";
    
    if (verbose) std::cout << "... done ode\n";
    
    if (verbose) 
      std::cout << std::endl << "Summary for t = " << t << std::endl
                <<              "=======================\n"
                << "\033[00;32mS-\033[0m" << " = " << rhs[0] << std::endl
                << "\033[00;31mI-\033[0m" << " = " << rhs[1] << std::endl
                << "\033[00;34mR-\033[0m" << " = " << rhs[2] << std::endl
                << "\033[01;32mS+\033[0m" << " = " << rhs[3] << std::endl
                << "\033[01;31mI+\033[0m" << " = " << rhs[4] << std::endl
                << "\033[01;34mR+\033[0m" << " = " << rhs[5] << std::endl
                << std::endl;
    
    // closing output file
    if (nsave > 0)
      output_file.close();
    
    // free allocated memory 
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
  } // solve 
  
  //------------------------------------------------------------
   
  template <class Params, class Eqs>
  const gsl_odeiv_step_type* OdeSolver<Params, Eqs>::set_step_type()
  {

    if (verbose)
      std::cout << "... setting stepping algorithm to "
                << step_algo << std::endl;
      
    if (step_algo == "rk2") return gsl_odeiv_step_rk2; 
    if (step_algo == "rk4") return gsl_odeiv_step_rk2; 
    if (step_algo == "rkf45") return gsl_odeiv_step_rkf45; 
    if (step_algo == "rkck") return gsl_odeiv_step_rkck; 
    if (step_algo == "rk8pd") return gsl_odeiv_step_rk8pd; 
    if (step_algo == "rk2imp") return gsl_odeiv_step_rk2imp; 
    if (step_algo == "rk4imp")return gsl_odeiv_step_rk4imp; 
    if (step_algo == "bsimp") return gsl_odeiv_step_bsimp; 
    if (step_algo == "gear1") return gsl_odeiv_step_gear1; 
    if (step_algo == "gear2") return gsl_odeiv_step_gear2; 
    
    if (verbose)
      std::cout << "... wrong step_algo value, set to default rkf45\n";     
    return gsl_odeiv_step_rkf45;
    
  } // get_step_type

  //------------------------------------------------------------
   
  template <class Params, class Eqs>
  void OdeSolver<Params, Eqs>::write_rhs(std::ofstream& ofile, double t)
  {
    int nv = model_params->nvars;
    
    ofile << t; 
    for(int i = 0; i < nv; i++)
      ofile << '\t' << rhs[i];
    ofile << std::endl;
      
  } // write_rhs 
   
  //------------------------------------------------------------

  template <class Params, class Eqs>
  void OdeSolver<Params, Eqs>::write_last_line()
  {    
    int nv = model_params->nvars;
    std::string lastLineFileName(file_id + ".final");
    
    // open file_id.final
    std::ofstream lastLine(lastLineFileName.c_str(), std::ios::out);
    
    // writing to stringstream
    for(int i = 0; i < nv; i++)
      lastLine << rhs[i] << '\t';
    lastLine << std::endl;
    
    lastLine.close();           
    
  } // write_last_line
  
  //------------------------------------------------------------
   
  // overloding ode::operator<<
  template <class Params, class Eqs>
  std::ostream& operator<< (std::ostream& output,
                            const OdeSolver<Params, Eqs>& x)   
  {
    output << std::endl
           << "Ode solver parameters:" << std::endl
           << "======================" << std::endl
           << "file id ............ " << x.get_file_id() << std::endl
           << "tmax ............... " << x.get_tmax() << std::endl
           << "nsave .............. " << x.get_nsave() << std::endl
           << "verbose mode ....... " << x.get_verbose() << std::endl
           << "convergence check .. " << x.get_convergence_check() << std::endl
           << "output file ........ " << x.get_file_id()+".dat" << std::endl
           << "ic file ............ " << x.get_ic_file_name() << std::endl
           << "dt ................. " << x.get_dt() << std::endl
           << "step type .......... " << x.get_step_algo() << std::endl
           << "abs tol ............ " << x.get_abs_tol() << std::endl
           << "rel tol ............ " << x.get_rel_tol() << std::endl
           << *(x.model_params);
            
    return output;
      
  } /* operator<< */
   
  //------------------------------------------------------------

}

//------------------------------------------------------------

#endif // ODE_SOLVER_H

