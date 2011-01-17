#include <iostream>

#include <boost/multi_array.hpp>

#include "tryps.hh"

//------------------------------------------------------------

double calculate_likelihoods(std::vector<double>& params,
                             boost::multi_array<unsigned int, 2> mixingStructure,
                             boost::multi_array<double, 2>& data)
{
  std::vector<double> FoI (params.size(), .0);
  for (unsigned int i = 0; i < mixingStructure.size(); ++i) {
    for (unsigned int j = 0; j < mixingStructure.size() ++j) {
      FoI[i] += params[mixingStructure[i][j] 
}

void loop_params(std::vector<double>& params, unsigned int idx, double& ml,
                 boost::multi_array<unsigned int, 2> mixingStructure,
                 boost::multi_array<double, 2>& data)
{
  if ((idx + 1) < params) {
    for (double p = .0; p < 10; p += 0.1) {
      params[idx] = p;
      loop_params(params, idx+1, ml, mixingStructure, data);
    }
  } else {
    for (double p = .0; p < 10; p += 0.1) {
      params[idx] = p;
      double lh = calculate_likelihoods(params, mixingStructure, data);
      if (lh > ml) ml = lh;
    }
  }
}

int main(int argc, char* argv[])
{

  boost::multi_array<unsigned int, 2> mixingStructure
    (boost::extents[38][38]);
  boost::multi_array<double, 2> mixingMatrix
    (boost::extents[38][38]);
  std::vector<double> params;

  // human, domestic, wildlife -- n*(n-1)=6 interactions

  // human-human
  mixingStructure[0][0] = 0;

  // human-domestic
  for (unsigned int i = 1; i < 5; ++i) {
    mixingStructure[0][i] = 1;
  }

  // human-wildlife
  for (unsigned int i = 5; i < 38; ++i) {
    mixingStructure[0][i] = 2;
  }

  // domestic-same domestic
  for (unsigned int i = 1; i < 5; ++i) {
    mixingStructure[i][i] = 3;
  }
  
  // domestic-different domestic
  for (unsigned int i = 1; i < 5; ++i) {
    for (unsigned int j = i+1; j < 5; ++j) {
      mixingStructure[i][j] = 4;
    }
  }
  
  // domestic-wildlife
  for (unsigned int i = 1; i < 5; ++i) {
    for (unsigned int j = 5; j < 38; ++j) {
      mixingStructure[i][j] = 5;
    }
  }

  // wildlife-same wildlife
  for (unsigned int i = 5; i < 38; ++i) {
    mixingStructure[i][i] = 6;
  }

  // wildlife-different wildlife
  for (unsigned int i = 5; i < 38; ++i) {
    for (unsigned int j = i+1; j < 38; ++j) {
      mixingStructure[i][j] = 7;
    }
  }

  // assume symmetry (to be relaxed)
  for (unsigned int j = 0; j < 38; ++j) {
    for (unsigned int i = j + 1; i < 38; ++i) {
      mixingStructure[i][j] = mixingStructure[j][i];
    }
  }

  // find maximum
  unsigned int max = 0;
  for (boost::multi_array<unsigned int, 2>::iterator cIt =
         mixingStructure.begin();
       cIt != mixingStructure.end(); cIt++) {
    for (boost::multi_array<unsigned int, 2>::subarray<1>::type::iterator
           rIt = cIt->begin(); rIt != cIt->end(); rIt++) {
      if (*rIt > max) max = *rIt;
    }
  }
  
  params.resize(max, 0);


  
  return 0;
}
   
//------------------------------------------------------------
