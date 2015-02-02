// -*- compile-command: "cd .. ; R CMD INSTALL .; cd -"; -*-
/*! \file rcpp_chronic.cpp
  \brief Implementation of the wrapper for running the chronic model form R
*/

#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>

#include "chronic_model.hpp"

// [[Rcpp::export]]
SEXP chronic(Rcpp::IntegerVector init,
             Rcpp::NumericVector params,
             Rcpp::IntegerVector times,
             Rcpp::IntegerVector stage1_passive = Rcpp::IntegerVector::create(),
             Rcpp::IntegerVector stage2_passive = Rcpp::IntegerVector::create(),
             Rcpp::NumericVector seed = R_NilValue,
             Rcpp::IntegerVector verbose = 0)
{
    std::map<std::string, double> cpp_param;
    std::map<std::string, int> cpp_init;

    // Get R variables
    Rcpp::CharacterVector initVariables = init.names();

    for (size_t i = 0; i < initVariables.length(); ++i) {
        cpp_init.emplace(std::string(initVariables[i]), init[i]);
    }

    Rcpp::CharacterVector paramVariables = params.names();

    for (size_t i = 0; i < paramVariables.length(); ++i) {
        cpp_param.emplace(std::string(paramVariables[i]), params[i]);
    }

    ChronicModel model(cpp_param, cpp_init,
                       Rcpp::as<std::vector<unsigned int> >(stage1_passive),
                       Rcpp::as<std::vector<unsigned int> >(stage2_passive));

    std::vector<int> cpp_times =
        Rcpp::as<std::vector<int> >(times);

    std::map<std::string, std::vector<double> > traj =
        model.Simulate(cpp_times);

    Rcpp::List Rtraj(traj.size());
    Rcpp::CharacterVector nameVector(traj.size());

    unsigned int count = 0;
    for (std::map<std::string, std::vector<double> >::iterator it =
             traj.begin(); it != traj.end(); ++it)
    {
        nameVector[count] = it->first;
        Rtraj[count] = it->second;
        count++;
    }

    Rtraj.attr("names") = nameVector;

    return(Rcpp::DataFrame(Rtraj));
}
