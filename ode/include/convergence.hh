#include <vector>

//------------------------------------------------------------

namespace convergence
{
  
  class ConvergenceCheck
  {

  public:
    
    ConvergenceCheck(unsigned int nv = 0, unsigned int smp_int = 1,
                     unsigned int nsmpls = 3, double e = 1e-5) :
      nvars(nv), samples_interval(smp_int), nsamples(nsmpls), eps(e) {}
    
    ~ConvergenceCheck() {}
    
    // mutators
    void set_nvars(const unsigned int nv) { nvars = nv; }
    void set_samples_interval(const unsigned int i) { samples_interval = i; }
    void set_nsamples(const unsigned int i) { nsamples = i; }
    void set_eps(const double e) { eps = e; }
    
    
    // accessors
    unsigned int get_nvars() const { return nvars; }
    unsigned int get_samples_interval() const { return samples_interval; }
    unsigned int get_nsamples() const { return nsamples; }
    double get_eps() const { return eps; }
    std::vector<std::vector<double> > get_rhs() const { return rhs; } 
    
    // main function - check convergence
    bool check(double* rhs_in);
    
  private:
    unsigned int nvars; // number of equations
    unsigned int samples_interval; // interval between rhs samples
    unsigned int nsamples; // how many rhs to account for
    double eps; // upper boundary
    
    // rhs container
    std::vector<std::vector<double> > rhs; // array of nsamples rhs
    
  }; // class ConvergenceCheck
  
  //------------------------------------------------------------
  
  bool ConvergenceCheck::check(double* rhs_in)
  {
    // convert new rhs to std::vector
    std::vector<double> tmp_vec(nvars);
    for (unsigned int i = 0; i < tmp_vec.size(); i++)
      tmp_vec[i] = rhs_in[i];
    
    // and insert it in front
    rhs.insert(rhs.begin(), tmp_vec);
    
    // cut rhs tail
    if (rhs.size() > nsamples)
      rhs.pop_back();
    
    // perform check
    bool converge = false;
    
    if (rhs.size() == nsamples) {
      double norm[rhs.size()-1];
      converge = true;
      
      // assign norm
      for (unsigned int i = 0; i < rhs.size()-1; i++)
        for (unsigned int j = 0; j < rhs[i].size(); j++) 
          norm[i] += (rhs[i][j] - rhs[i+1][j])*(rhs[i][j] - rhs[i+1][j]);
      
      // compare to eps
      for (unsigned int i = 0; i < rhs.size()-1; i++)
        if (norm[i] > eps)
          converge = false;
      
//       for (unsigned int i = 0; i < rhs.size()-1; i++)
//         std::cout << norm[i] << '\t';
//       std::cout << std::endl << std::endl;
      
    }
    
    return converge;
    
  }

  
} // namespace convergence


  
