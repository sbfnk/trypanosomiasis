#ifndef EVENT_HPP
#define EVENT_HPP

#include <vector>

struct Event
{
    Event(double r, std::vector<std::string> m =
          std::vector<std::string>()) :
        rate(r), multipliers(m) {;}

    void addMultiplier(std::string m) { multipliers.push_back(m); }
    void addAction(std::string state, int change)
    {
        actions.push_back(std::make_pair(state, change));
    }

    double rate;
    std::vector<std::string> multipliers;
    std::vector<std::pair<std::string, int> > actions;
};

#endif
