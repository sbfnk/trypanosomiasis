/*! \file chronic_model.hpp
  \brief Header file for the chronic model
*/

#ifndef CHRONIC_MODEL_HPP
#define CHRONIC_MODEL_HPP

#include <map>
#include <vector>

#include "event.hpp"

class ChronicModel {

public:

    ChronicModel(std::map<std::string, double> params,
                 std::map<std::string, int> init,
                 std::vector<unsigned int> p1s =
                 std::vector<unsigned int>(),
                 std::vector<unsigned int> p2s =
                 std::vector<unsigned int>(),
                 std::vector<unsigned int> verbose =
                 std::vector<unsigned int>());

    std::map<std::string, std::vector<double> >
    Simulate(std::vector<int> times);

    std::map<std::string, std::vector<double> >
    Simulate(std::vector<int> times, unsigned int seed);

    double generateEventList(std::vector<double> & result);

private:

    std::map<std::string, int> states;
    std::map<std::string, double> rates;

    std::vector<Event> eventList;

    unsigned int verbose;
    
    std::vector<unsigned int> passive_stage1;
    std::vector<unsigned int> passive_stage2;

};

#endif
