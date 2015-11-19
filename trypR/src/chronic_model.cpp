#include <sys/time.h>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "chronic_model.hpp"

ChronicModel::ChronicModel(std::map<std::string, double> params,
                           std::map<std::string, int> init,
                           std::vector<unsigned int> p1s,
                           std::vector<unsigned int> p2s,
                           std::vector<unsigned int> verbose) :
    states(init), rates(params), passive_stage1(p1s), passive_stage2(p2s),
    verbose(verbose.size() > 0 ? verbose[0] : 0)
{
    Event chronic_infection(params["lambda"] * params["pc"],
                            std::vector<std::string>(1, "S"));
    chronic_infection.addAction("S", -1);
    chronic_infection.addAction("Ic", +1);
    eventList.push_back(chronic_infection);

    Event chronic_infection_from_Ic(params["lambda"] * params["pc"] *
                                    params["delta"],
                                    std::vector<std::string>(1, "S"));
    chronic_infection_from_Ic.addAction("S", -1);
    chronic_infection_from_Ic.addAction("Ic", +1);
    chronic_infection_from_Ic.addMultiplier("Ic");
    eventList.push_back(chronic_infection_from_Ic);

    Event chronic_infection_from_I1(params["lambda"] * params["pc"],
                                    std::vector<std::string>(1, "S"));
    chronic_infection_from_I1.addAction("S", -1);
    chronic_infection_from_I1.addAction("Ic", +1);
    chronic_infection_from_I1.addMultiplier("I1");
    eventList.push_back(chronic_infection_from_I1);

    Event stage1_infection(params["lambda"] * params["pc"],
                            std::vector<std::string>(1, "S"));
    stage1_infection.addAction("S", -1);
    stage1_infection.addAction("I1", +1);
    eventList.push_back(stage1_infection);

    Event stage1_infection_from_Ic(params["beta"] * (1 - params["pc"]) *
                                   params["delta"],
                                   std::vector<std::string>(1, "S"));
    stage1_infection_from_Ic.addAction("S", -1);
    stage1_infection_from_Ic.addAction("I1", +1);
    stage1_infection_from_Ic.addMultiplier("Ic");
        eventList.push_back(stage1_infection_from_Ic);

    Event stage1_infection_from_I1(params["beta"] * (1 - params["pc"]),
                                   std::vector<std::string>(1, "S"));
    stage1_infection_from_I1.addAction("S", -1);
    stage1_infection_from_I1.addAction("I1", +1);
    stage1_infection_from_I1.addMultiplier("I1");
    eventList.push_back(stage1_infection_from_I1);

    Event chronic_recovery(params["rc"],
                           std::vector<std::string>(1, "Ic"));
    chronic_recovery.addAction("S", +1);
    chronic_recovery.addAction("Ic", -1);
    eventList.push_back(chronic_recovery);

    Event stage1_progression(params["r1"],
                             std::vector<std::string>(1, "I1"));
    stage1_progression.addAction("I2", +1);
    stage1_progression.addAction("I1", -1);
    eventList.push_back(stage1_progression);

    Event stage2_removal(params["r2"],
                         std::vector<std::string>(1, "I2"));
    stage2_removal.addAction("S", +1);
    stage2_removal.addAction("I2", -1);
    eventList.push_back(stage2_removal);

    if (params.find("p1") != params.end()) {
        Event stage1_passive_detection(params["p1"],
                                       std::vector<std::string>(1, "I1"));
        stage1_passive_detection.addAction("I1", -1);
        stage1_passive_detection.addAction("S", +1);
        stage1_passive_detection.addAction("Z1pass", +1);
        eventList.push_back(stage1_passive_detection);
        states["Z1pass"] = 0;
    }

    if (params.find("p2") != params.end()) {
        Event stage2_passive_detection(params["p2"],
                                       std::vector<std::string>(1, "I2"));
        stage2_passive_detection.addAction("I2", -1);
        stage2_passive_detection.addAction("S", +1);
        stage2_passive_detection.addAction("Z2pass", +1);
        eventList.push_back(stage2_passive_detection);
        states["Z2pass"] = 0;
    }

}

std::map<std::string, std::vector<double> >
ChronicModel::Simulate(std::vector<int> times)
{
    unsigned int seed;
    struct timeval tv;
    gettimeofday(&tv, 0);
    seed = tv.tv_sec + tv.tv_usec;

    return (Simulate(times, seed));
}

std::map<std::string, std::vector<double> >
ChronicModel::Simulate(std::vector<int> times, unsigned int seed)
{
    boost::mt19937 randGen(seed);
    boost::uniform_01<boost::mt19937> gen(randGen);

    int next_time = -1;
    int last_time = -1;

    if (times.size() > 0) {
        next_time = 0;
        last_time = times.back();
    }

    std::map<std::string, std::vector<double> > traj;
    traj.emplace("time", std::vector<double>());
    traj.emplace("S", std::vector<double>());
    traj.emplace("Ic", std::vector<double>());
    traj.emplace("I1", std::vector<double>());
    traj.emplace("I2", std::vector<double>());
    if (rates.find("p1") != rates.end()) traj.emplace("Z1pass", std::vector<double>());
    if (rates.find("p2") != rates.end()) traj.emplace("Z2pass", std::vector<double>());

    double time = 0;
    while (time < last_time)
    {
        while (times[next_time] <= time)
        {
            if (static_cast<int>(passive_stage1.size()) > next_time)
            {
                states["I1"] -= passive_stage1[next_time];
            }
            if (static_cast<int>(passive_stage2.size()) > next_time)
            {
                states["I2"] -= passive_stage2[next_time];
            }

            traj["time"].push_back(times[next_time]);
            for (std::map<std::string, int>::iterator it =
                     states.begin(); it != states.end(); ++it)
            {
                traj[it->first].push_back(it->second);
            }

            if (states["I1"] < 0) states["I1"] = 0;
            if (states["I2"] < 0) states["I2"] = 0;
            ++next_time;
        }

        std::vector<double> eventRates;

        // event time
        double eventSum = generateEventList(eventRates);

        if (verbose) {
            std::cout << std::endl << "Time " << time << ": ";
            for (std::vector<double>::iterator it = eventRates.begin();
                 it != eventRates.end(); ++it)
            {
                std::cout << *it << " ";
            }
            std::cout << "-- " << eventSum << std::endl;
        }
        double dt = -log(gen()) / eventSum;
        time += dt;

        // choose event
        double randEvent = gen() * eventSum;

        if (verbose) {
            std::cout << "Random event " << randEvent << ": ";
        }

        int chosenEvent = 0;
        double rateAcc = 0;

        while ((chosenEvent < static_cast<int>(eventList.size())) &
               (rateAcc + eventRates[chosenEvent] < randEvent))
        {
            rateAcc += eventRates[chosenEvent];
            chosenEvent++;
        }

        if (verbose) {
            std::cout << chosenEvent << std::endl;
        }

        if (chosenEvent >= 0)
        {
            for (std::vector<std::pair<std::string, int> >::iterator
                     it = eventList[chosenEvent].actions.begin();
                 it != eventList[chosenEvent].actions.end(); ++it)
            {
                states[it->first] += it->second;
            }
        }
    }
    while ((next_time < static_cast<int>(times.size())) &
           (times[next_time] <= time))
    {
        traj["time"].push_back(times[next_time]);
        for (std::map<std::string, int>::iterator it = states.begin();
             it != states.end(); ++it)
        {
            traj[it->first].push_back(it->second);
        }
        ++next_time;
    }

    return(traj);
}

double ChronicModel::generateEventList(std::vector<double>& result)
{
    double sum(0);

    for (std::vector<Event>::iterator it = eventList.begin();
         it != eventList.end(); ++it)
    {
        double rate = it->rate;
        for (std::vector<std::string>::iterator it2 = it->multipliers.begin();
             it2 != it->multipliers.end(); ++it2)
        {
            rate *= states[*it2];
        }
        result.push_back(rate);
        sum += rate;
    }

    return(sum);
}
