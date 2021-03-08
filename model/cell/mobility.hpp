#ifndef CISE_PANDEMIC_MOBILITY_HPP
#define CISE_PANDEMIC_MOBILITY_HPP

#include <numeric>
#include <cadmium/celldevs/utils/grid_utils.hpp>
#include <cadmium/json/json.hpp>
#include "state.hpp"

using namespace cadmium::celldevs;


struct MobilityReduction {
    virtual ~MobilityReduction() = default;
    virtual std::vector<double> mobility_reduction(double simulation_clock, cell_position const &cell, sird const &cell_state) const {
        return std::vector<double>(cell_state.susceptible.size(), 1);
    }
};

class PeriodicMobilityReduction: public MobilityReduction {
    mutable int first_day;
public:
    const int time_scale;
    const std::vector<std::vector<double>> lockdown_rates;
    const std::vector<int> phase_durations;
    int time_sum;

    PeriodicMobilityReduction(int time_scale, std::vector<std::vector<double>> &lr, std::vector<int> &pd, bool scramble):
        time_scale(time_scale), lockdown_rates(lr), phase_durations(pd) {
        time_sum = 0;
        for(int phase_duration: phase_durations) {
            time_sum += phase_duration;
        }
        first_day = (scramble) ? -1 : 0;
    }

    std::vector<double> mobility_reduction(double simulation_clock, cell_position const &cell_id, sird const &cell_state) const override {
        if (first_day < 0) {
            first_day = 0;
            int first_phase = accumulate(cell_id.begin(),cell_id.end(), 0) % (int)phase_durations.size();
            for (int i = 0; i < first_phase; ++i) {
                first_day += phase_durations.at(i);
            }
        }
        int phase_day = (int)(first_day + simulation_clock * time_scale) % time_sum;
        int phase_group = 0;
        while (phase_day >= phase_durations.at(phase_group)) {
            phase_day -= phase_durations.at(phase_group++);
        }
        phase_group = (phase_group + accumulate(cell_id.begin(),cell_id.end(),0)) % (int)phase_durations.size();
        return lockdown_rates.at(phase_group);
    }
};


struct InfectedMobilityReduction: MobilityReduction {
    const std::vector<std::vector<double>> infected_effect;

    InfectedMobilityReduction(std::vector<std::vector<double>> &ie): infected_effect(ie) {}

    std::vector<double> mobility_reduction(double simulation_clock, cell_position const &cell_id, sird const &cell_state) const override {
        std::vector<double> lockdown_factors = {};
        for(int i = 0; i < cell_state.infected.size(); i++) {
            lockdown_factors.push_back(std::max(1 - cell_state.infected_ratio() * infected_effect.at(0).at(i), 0.0));
        }
        return lockdown_factors;
    }
};


class InfectedPhaseMobilityReduction: public MobilityReduction {
    mutable int phase;
    const std::vector<std::vector<double>> lockdown_rates;
    const std::vector<double> phase_thresholds;
    const std::vector<double> threshold_buffers;

    bool shouldGoToNextPhase(sird const &cell_state) const {
        return (phase < phase_thresholds.size() - 1 && cell_state.infected_ratio() >= phase_thresholds[phase + 1]);
    }

    bool shouldGoToPreviousPhase(sird const &cell_state) const {
        return (phase > 0 && (cell_state.infected_ratio() + threshold_buffers[phase]) < phase_thresholds[phase]);
    }

    void updatePhase(sird const &cell_state) const {
        if(shouldGoToNextPhase(cell_state)) {
            ++phase;
        } else if (shouldGoToPreviousPhase(cell_state)) {
            --phase;
        }
    }

public:
    InfectedPhaseMobilityReduction(std::vector<std::vector<double>> &lr, std::vector<double> &pt, std::vector<double> &tb):
            lockdown_rates(lr), phase_thresholds(pt), threshold_buffers(tb), phase(0) {}

    std::vector<double> mobility_reduction(double simulation_clock, cell_position const &cell_id, sird const &cell_state) const override {
        updatePhase(cell_state);
        return lockdown_rates.at(phase);
    }
};

#endif //CISE_PANDEMIC_MOBILITY_HPP
