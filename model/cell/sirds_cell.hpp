/**
 * Copyright (c) 2020, Román Cárdenas Rodríguez
 * ARSLab - Carleton University
 * GreenLSI - Polytechnic University of Madrid
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CISE_PANDEMIC_SIRDS_CELL_HPP
#define CISE_PANDEMIC_SIRDS_CELL_HPP

#include <cmath>
#include <random>
#include <cadmium/json/json.hpp>
#include <cadmium/celldevs/cell/grid_cell.hpp>

#include "state.hpp"
#include "vicinity.hpp"

using cadmium::json;
using namespace cadmium::celldevs;

struct sirds_config {

    std::vector<double> susceptibility;
    std::vector<double> virulence;
    std::vector<double> recovery;
    std::vector<double> mortality;
    std::vector<double> immunity;
    int precision;
    int time_scaler;

    sirds_config(): susceptibility({1.0}), virulence({0.0}), recovery({0.0}), mortality({0.0}), immunity({1.0}), precision(100), time_scaler(1) {}

    sirds_config(std::vector<double> &s, std::vector<double> &v, std::vector<double> &r, std::vector<double> &m, std::vector<double> &i, int p, int t):
            susceptibility(s), virulence(v), recovery(r), mortality(m), immunity(i), precision(p), time_scaler(t) {}
};

void from_json(const cadmium::json& j, sirds_config &c) {
    j.at("susceptibility").get_to(c.susceptibility);
    j.at("virulence").get_to(c.virulence);
    j.at("recovery").get_to(c.recovery);
    j.at("mortality").get_to(c.mortality);
    j.at("immunity").get_to(c.immunity);
    uint n_decimals = (j.contains("n_decimals")) ? j["n_decimals"].get<uint>() : 2;
    c.precision = (int) pow(10, n_decimals);
    int time_scaler = (j.contains("time_scaler"))? j["time_scaler"].get<int>() : 1;
    c.time_scaler = (time_scaler > 0) ? time_scaler : 1;
}

template <typename T>
class sirds_cell : public grid_cell<T, sird, mc> {
public:

    using grid_cell<T, sird, mc>::cell_id;
    using grid_cell<T, sird, mc>::simulation_clock;
    using grid_cell<T, sird, mc>::state;
    using grid_cell<T, sird, mc>::map;
    using grid_cell<T, sird, mc>::neighbors;

	std::vector<double> susceptibility;
	std::vector<double> virulence;
	std::vector<double> recovery;
	std::vector<double> mortality;
	std::vector<double> immunity;
	int precision = 100;
	int time_scaler = 1;

	sirds_cell() : grid_cell<T, sird, mc>()  {}

	sirds_cell(cell_position const &cell_id, cell_unordered<mc> const &neighborhood, sird &initial_state,
               cell_map<sird, mc> const &map_in, std::string const &delay_id, sirds_config &config) :
			    grid_cell<T, sird, mc>(cell_id, neighborhood, initial_state, map_in, delay_id) {
		susceptibility = config.susceptibility;
		virulence = config.virulence;
		recovery = config.recovery;
		mortality = config.mortality;
        immunity = config.immunity;
		precision = config.precision;
		time_scaler = config.time_scaler;
		auto s = state.current_state;
	}

	[[nodiscard]] unsigned int inline n_age_segments() const {
		return state.current_state.population.size();
	}

	// user must define this function. It returns the next cell state and its corresponding timeout
	[[nodiscard]] sird local_computation() const override {
		auto res = state.current_state;
		auto new_i = new_infections(res);
		auto new_r = new_recoveries(res);
		auto new_d = new_deaths(res);
		auto new_s = new_susceptible(res);

		for (int i = 0; i < n_age_segments(); i++) {
			res.recovered[i] = std::round((res.recovered[i] + new_r[i] - new_s[i]) * precision) / precision;
			res.deceased[i] = std::round((res.deceased[i] + new_d[i]) * precision) / precision;
			res.infected[i] = std::round((res.infected[i] + new_i[i] - (new_r[i] + new_d[i])) * precision) / precision;
			res.susceptible[i] = 1 - (res.recovered[i] + res.infected[i] + res.deceased[i]);
			// We avoid any possible negative value due to rounding
			if (res.susceptible[i] < 0) {
                res.infected[i] += res.susceptible[i];
                res.susceptible[i] = 0;
                if (res.infected[i] < 0) {
                    res.recovered[i] += res.infected[i];
                    res.infected[i] = 0;
                    if (res.recovered[i] < 0) {
                        res.deceased[i] += res.recovered[i];
                        res.recovered[i] = 0;
                        if (res.deceased[i] < 0) {
                            throw std::bad_exception();
                        }
                    }
                }
			}
		}
		return res;
	}

	// It returns the delay to communicate cell's new state.
	T output_delay(sird const &cell_state) const override { return ((float)1) / time_scaler; }

	[[nodiscard]] std::vector<double> new_infections(sird const &last_state) const {
        double n_effect = 0;
        for (auto neighbor: neighbors) {
            sird n_state = state.neighbors_state.at(neighbor);
            mc n_vicinity = state.neighbors_vicinity.at(neighbor);
            for (int k = 0; k < n_age_segments(); k++) {
                n_effect += n_state.population[k] * n_state.infected[k] * n_vicinity.connectivity * n_vicinity.mobility[k] * virulence[k];
            }
        }

        auto new_inf = std::vector<double>();
        for (int n = 0; n < n_age_segments(); n++) {
            double ratio = std::min(1.0, n_effect / last_state.population[n]);
            new_inf.push_back(last_state.susceptible[n] * susceptibility[n] * ratio);
        }
        return new_inf;
	}

	[[nodiscard]] std::vector<double> new_recoveries(sird const &last_state) const {
		auto new_r = std::vector<double>();
		for (int i = 0; i < n_age_segments(); i++) {
			new_r.push_back(last_state.infected[i] * recovery[i]);
		}
		return new_r;
	}

	[[nodiscard]] std::vector<double> new_deaths(sird const &last_state) const {
		auto new_d = std::vector<double>();
		for (int i = 0; i < n_age_segments(); i++) {
			new_d.push_back(last_state.infected[i] * mortality[i]);
		}
		return new_d;
	}

    [[nodiscard]] std::vector<double> new_susceptible(sird const &last_state) const {
        auto new_s = std::vector<double>();
        for (int i = 0; i < n_age_segments(); i++) {
            new_s.push_back(last_state.recovered[i] * (1 - immunity[i]));
        }
        return new_s;
    }
};

#endif //CISE_PANDEMIC_SIRDS_CELL_HPP
