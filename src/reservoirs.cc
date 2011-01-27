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
      host newHost(lineVector, headings);
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
      vector newVector(lineVector, headings);
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

  betafunc_params p(

  
  int minIndex = 0;

  for (unsigned int i = 0; i < years.size(); ++i) {
    if (highYear == years[i]) {
      minIndex = i;
    }
  }

  for (std::vector<std::vector<std::string> >::iterator it = data.begin();
       it != data.end(); it++) {

    std::string currentCountry = (*it)[0];
    if (country.length() == 0 || country == currentCountry) {
    
      std::vector<double> xValues;
      std::vector<double> yValues;
      std::vector<double> ySqrtValues;

      std::stringstream line;
      if (name.length() > 0) {
        line << name << ":";
      } else {
        line << currentCountry << ":";
      }

      for (int i = minIndex + maxYears - 1; i > (minIndex - 1); --i) {
        int j = 0;
        std::istringstream s((*it)[i+1]); // index 0 is the country
        s >> j;
        if (j > 0) {
          xValues.push_back(years[i]);
          yValues.push_back(j);
          ySqrtValues.push_back(sqrt(j));
        }
      }

      if (verbose >= 2) {
        std::cout << "VALUES(x) (" << xValues.size() << ") ";
        for (unsigned int i = 0; i < xValues.size(); ++i) {
          std::cout << xValues[i] << " ";
        }
        std::cout << std::endl;
        
        std::cout << "VALUES(y) (" << yValues.size() << ") ";
        for (unsigned int i = 0; i < yValues.size(); ++i) {
          std::cout << yValues[i] << " ";
        }
        std::cout << std::endl;

        if (verbose >= 2) {
          std::cout << "VALUES(sqrt(y)) (" << ySqrtValues.size() << ") ";
          for (unsigned int i = 0; i < ySqrtValues.size(); ++i) {
            std::cout << ySqrtValues[i] << " ";
          }
          std::cout << std::endl;
        }
      }

      if (yValues.size() > 1) {

        int country_maxBreakpoints;
        
        country_maxBreakpoints =
          std::min(static_cast<int>(xValues.size() / minYears) - 1,
                   maxBreakpoints);

        bool emerged = false, extinct = false;
        if (xValues[0] - minYears >= years[minIndex + maxYears - 1]) {
          line << " emerged in " << xValues[0] << ",";
          emerged = true;
        }
        if (xValues[xValues.size() - 1] + minYears <= years[minIndex]) {
          line << " extinct since " << (xValues[xValues.size() - 1] + 1) << ",";
          extinct = true;
        }

        if (!(emerged && extinct)) {
          bool emerging = false;
          bool receding = false;
          if (regression) {
            
            double null_rSquared = -1;
            double alt_rSquared = -1;
            std::vector<double > p_null_rSquared (permutations - 1, -1);

            int low_bp = 0;	// minimum: 0 breakpoints
            int high_bp = country_maxBreakpoints > 0 ? country_maxBreakpoints : -1; // maximum

            bool do_low = true;

            std::vector<std::vector<double> > permuted_yValues;
            std::vector<double> null_bp, null_residuals, null_params;
            std::vector<double> alt_bp, alt_residuals, alt_params;

            while (low_bp < high_bp || high_bp == -1) {
              if (do_low) {
                null_rSquared = -1;
                std::vector<double> null_starting_bp;
                create_bp_vector(low_bp, xValues, null_starting_bp, minYears);

                grid_search(null_starting_bp, yValues, null_rSquared,
                            null_bp, null_residuals, null_params, minYears,
                            verbose);

                if (verbose >= 2) {
                  std::cout << print_results(null_rSquared, null_bp, null_params,
                                             null_residuals, xValues, &yValues);
                  // std::vector<double> y(xValues.size());
                  // std::vector<double> best_bp(1, null_bp.size());
                  // for (std::vector<double>::iterator it = null_bp.begin();
                  //      it != null_bp.end(); it++) {
                  //   best_bp.push_back(*it);
                  // }
                  // for (std::vector<double>::iterator it = xValues.begin();
                  //      it != xValues.end(); it++) {
                  //   best_bp.push_back(*it);
                  // }
                
                  // segmented(&null_params[0], &y[0], null_params.size(),
                  //           xValues.size(), &best_bp[0]);

                  // for (unsigned int i = 0; i < xValues.size(); ++i) {
                  //   std::cout << xValues[i] << " " << y[i] << std::endl;
                  // }
                }
              
                if (high_bp > 0) {
                  // don't need to permutations if I
                  // have nothing to compare to anyway

                  std::vector<double> scaled_residuals;
                  if (verbose >= 2) {
                    std::cout << "SCALED RESIDUALS (" << null_residuals.size()
                              << ") ";
                  }
                  for (unsigned int i = 0; i < null_residuals.size(); ++i) {
                    scaled_residuals.push_back(null_residuals[i] / ySqrtValues[i]);
                    if (verbose >= 2) {
                      std::cout << scaled_residuals[scaled_residuals.size() - 1] << " ";
                    }
                  }
                  if (verbose >= 2) {
                    std::cout << std::endl;
                  }
                
                  permuted_yValues.clear();
                  for (unsigned int i = 0; i < (permutations - 1); ++i) {
                    permuted_yValues.push_back(scaled_residuals);
                    std::random_shuffle(permuted_yValues[i].begin(), permuted_yValues[i].end());
                    if (verbose >= 2) {
                      std::cout << "PERMUTED (" << permuted_yValues[i].size() << ") ";
                    }
                    for (unsigned int j = 0; j < permuted_yValues[i].size(); ++j) {
                      permuted_yValues[i][j] *= ySqrtValues[j];
                      permuted_yValues[i][j] += yValues[j];
                      if (verbose >= 2) {
                        std::cout << permuted_yValues[i][j] << " ";
                      }
                    }
                    if (verbose >= 2) {
                      std::cout << std::endl;
                    }
                  }

                  for (std::vector<double>::iterator it = p_null_rSquared.begin();
                       it != p_null_rSquared.end(); it++) {
                    (*it) = -1;
                  }
                  std::vector<std::vector<double> > p_null_bp (permutations - 1);
                  std::vector<std::vector<double> > p_null_residuals (permutations - 1);
                  std::vector<std::vector<double> > p_null_params (permutations - 1);

                  create_bp_vector(low_bp, xValues, null_starting_bp, minYears);

                  grid_search(null_starting_bp, permuted_yValues, p_null_rSquared,
                              p_null_bp, p_null_residuals, p_null_params,
                              minYears, verbose);
              
                  if (verbose >= 2) {
                    std::cout << print_results(p_null_rSquared, p_null_bp,
                                               p_null_params, p_null_residuals,
                                               xValues);
                  }
                }
              }

              // now test alternative model (if we do so)
              if (high_bp > 0) {

                alt_rSquared = -1;
            
                std::vector<double> alt_starting_bp;
                create_bp_vector(high_bp, xValues, alt_starting_bp, minYears);

                grid_search(alt_starting_bp, yValues, alt_rSquared,
                            alt_bp, alt_residuals, alt_params,
                            minYears, verbose);
              
                if (verbose >= 2) {
                  std::cout << print_results(alt_rSquared, alt_bp,
                                             alt_params, alt_residuals,
                                             xValues, &yValues);
                  // std::vector<double> y(xValues.size());
                  // std::vector<double> best_bp(1, alt_bp.size());
                  // for (std::vector<double>::iterator it = alt_bp.begin();
                  //      it != alt_bp.end(); it++) {
                  //   best_bp.push_back(*it);
                  // }
                  // for (std::vector<double>::iterator it = xValues.begin();
                  //      it != xValues.end(); it++) {
                  //   best_bp.push_back(*it);
                  // }
                
                  // segmented(&alt_params[0], &y[0], null_params.size(),
                  //           xValues.size(), &best_bp[0]);

                  // for (unsigned int i = 0; i < xValues.size(); ++i) {
                  //   std::cout << xValues[i] << " " << y[i] << std::endl;
                  // }
                }
                std::vector<std::vector<double> > p_alt_bp (permutations - 1);
                std::vector<std::vector<double> > p_alt_residuals (permutations - 1);
                std::vector<std::vector<double> > p_alt_params (permutations - 1);
                std::vector<double> p_alt_rSquared (permutations - 1, -1);
              
                create_bp_vector(high_bp, xValues, alt_starting_bp, minYears);
                grid_search(alt_starting_bp, permuted_yValues, p_alt_rSquared,
                            p_alt_bp, p_alt_residuals, p_alt_params,
                            minYears, verbose);
              
                if (verbose >= 2) {
                  std::cout << print_results(p_alt_rSquared, p_alt_bp,
                                             p_alt_params, p_alt_residuals,
                                             xValues);
                }
            
                unsigned int p_count = 1;
                double compare = null_rSquared / alt_rSquared;
                for (unsigned int i = 0; i < permutations - 1; ++i) {
                  if (verbose >= 2) {
                    std::cout << "COMPARE " << p_null_rSquared[i] / p_alt_rSquared[i]
                              << " " << compare;
                  }
                  if (p_null_rSquared[i] / p_alt_rSquared[i] >= compare) {
                    if (verbose >= 2) {
                      std::cout << " (+)";
                    }
                    ++p_count;
                  }
                  if (verbose >= 2) {
                    std::cout << std::endl;
                  }
                }
              
                double pValue = p_count / static_cast<double>(permutations);
                if (verbose >= 1) {
                  std::cout << low_bp << " vs " << high_bp << ": p=" << pValue
                            << std::endl;
                }

                if (pValue > (alpha / high_bp)) {
                  // accepted
                  do_low = false;
                  --high_bp;
                } else {
                  // rejected
                  do_low = true;
                  ++low_bp;
                }
              } else {
                high_bp = -2;
              }
            }

            std::vector<double> *params;
            std::vector<double> *bp;

            if (do_low && high_bp > 0) {
              bp = &alt_bp;
              params = &alt_params;
            } else {
              bp = &null_bp;
              params = &null_params;
            }
            if (verbose >= 1) {
              std::cout << "Number of breakpoints: " << bp->size();
              if (bp->size() > 0) {
                std::cout << " ("
                          << xValues[static_cast<unsigned int>((*bp)[0])];
                for (std::vector<double>::iterator it = bp->begin() + 1;
                     it != bp->end(); it++) {
                  std::cout << " " << xValues[static_cast<unsigned int>(*it)];
                }
                std::cout << ")";
              }
              std::cout << std::endl;
            
              std::cout << "Parameters: ";
              for (std::vector<double>::iterator it = params->begin();
                   it != params->end(); it++) {
                std::cout << *it << " ";
              }
              std::cout << std::endl;
            }
          
            // calculate slope
            double paramSum = 0;
            for (std::vector<double>::iterator it = params->begin() + 1;
                 it != params->end(); it++) {
              paramSum += (*it);
            }

            // calculate average over that period
            unsigned int periodSum = 0;
            unsigned int lastBp = 0;
            if (bp->size() > 0) {
              lastBp = static_cast<unsigned int>((*bp)[bp->size() - 1]);
            }
            for (unsigned int i = lastBp + 1; i < xValues.size(); ++i) {
              periodSum += yValues[i];
            }

            double caseAvg = periodSum /
              static_cast<double>(xValues.size() - lastBp - 1);
          
            if (caseAvg < 100) caseAvg = 100;
          
            if (verbose) {
              std::cout << "Slope " << paramSum
                        << ", average " << caseAvg << std::endl;
            }
            
            if (paramSum / caseAvg > threshold) {
              emerging = true;
            }
            if (paramSum / caseAvg < (-threshold)) {
              receding = true;
            }

            if (emerging) {
              line << " emerging ";
              if (bp->size() > 0) {
                line << "since " << xValues[static_cast<int>((*bp)[bp->size()-1])]
                     << " ";
              }              
            } else if (!extinct) {
              if (receding) {
                line << " receding ";
                if (bp->size() > 0) {
                  line << "since " << xValues[static_cast<int>((*bp)[bp->size()-1])]
                       << " ";
                }
              } else {
                line <<  " stable ";
                if (bp->size() > 0) {
                  line << "since " << xValues[static_cast<int>((*bp)[bp->size()-1])]
                       << " ";
                }
              }
            }
          }
        }
        std::string outStr = line.str();
        std::cout << outStr.substr(0, outStr.size()-1) << std::endl;
      }
    }
  }
}
   
//------------------------------------------------------------
