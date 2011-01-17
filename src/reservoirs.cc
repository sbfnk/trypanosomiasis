#include <iostream>

#include <boost/multi_array.hpp>

#include "reservoirs.hh"

//------------------------------------------------------------

int main(int argc, char* argv[])
{

  po::options_description main_options
    ("Usage: reservoirs [options]... \n\nOptions");

  main_options.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("country,c", po::value<std::string>(),
     "choose country")
    ;

  po::positional_options_description file_option;
  file_option.add("input-file", -1);

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
  
  
  // boost::multi_array<unsigned int, 2> mixingStructure
  //   (boost::extents[38][38]);
  // boost::multi_array<double, 2> mixingMatrix
  //   (boost::extents[38][38]);
  // std::vector<double> params;

  // // human, domestic, wildlife -- n*(n-1)=6 interactions

  // // human-human
  // mixingStructure[0][0] = 0;

  // // human-domestic
  // for (unsigned int i = 1; i < 5; ++i) {
  //   mixingStructure[0][i] = 1;
  // }

  // // human-wildlife
  // for (unsigned int i = 5; i < 38; ++i) {
  //   mixingStructure[0][i] = 2;
  // }

  // // domestic-same domestic
  // for (unsigned int i = 1; i < 5; ++i) {
  //   mixingStructure[i][i] = 3;
  // }
  
  // // domestic-different domestic
  // for (unsigned int i = 1; i < 5; ++i) {
  //   for (unsigned int j = i+1; j < 5; ++j) {
  //     mixingStructure[i][j] = 4;
  //   }
  // }
  
  // // domestic-wildlife
  // for (unsigned int i = 1; i < 5; ++i) {
  //   for (unsigned int j = 5; j < 38; ++j) {
  //     mixingStructure[i][j] = 5;
  //   }
  // }

  // // wildlife-same wildlife
  // for (unsigned int i = 5; i < 38; ++i) {
  //   mixingStructure[i][i] = 6;
  // }

  // // wildlife-different wildlife
  // for (unsigned int i = 5; i < 38; ++i) {
  //   for (unsigned int j = i+1; j < 38; ++j) {
  //     mixingStructure[i][j] = 7;
  //   }
  // }

  // // assume symmetry (to be relaxed)
  // for (unsigned int j = 0; j < 38; ++j) {
  //   for (unsigned int i = j + 1; i < 38; ++i) {
  //     mixingStructure[i][j] = mixingStructure[j][i];
  //   }
  // }

  // // find maximum
  // unsigned int max = 0;
  // for (boost::multi_array<unsigned int, 2>::iterator cIt =
  //        mixingStructure.begin();
  //      cIt != mixingStructure.end(); cIt++) {
  //   for (boost::multi_array<unsigned int, 2>::subarray<1>::type::iterator
  //          rIt = cIt->begin(); rIt != cIt->end(); rIt++) {
  //     if (*rIt > max) max = *rIt;
  //   }
  // }
  
  // params.resize(max, 0);


  
  return 0;
}
   
//------------------------------------------------------------
