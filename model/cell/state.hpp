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

#ifndef CISE_PANDEMIC_STATE_HPP
#define CISE_PANDEMIC_STATE_HPP

#include <iostream>
#include <cadmium/json/json.hpp>

struct sird {
    std::vector<int> population;
    std::vector<double> susceptible;
    std::vector<double> infected;
    std::vector<double> recovered;
    std::vector<double> deceased;
    std::vector<double> mobility;

    sird() : population({0}), susceptible({1}), infected({0}), recovered({0}), deceased({0}), mobility({1}){}
    sird(std::vector<int> &pop, std::vector<double> &s, std::vector<double> &i, std::vector<double> &r, std::vector<double> &d, std::vector<double> &m) :
            population(pop), susceptible(s), infected(i), recovered(r), deceased(d), mobility(m){}

    [[nodiscard]] int total_population() const {
        int res = 0;
        for (auto const &i: population)
            res += i;
        return res;
    }

    [[nodiscard]] double susceptible_ratio(int i) const {
        return susceptible[i] * population[i] / total_population();
    }

    [[nodiscard]] double infected_ratio(int i) const {
        return infected[i] * population[i] / total_population();
    }

    [[nodiscard]] double recovered_ratio(int i) const {
        return recovered[i] * population[i] / total_population();
    }

    [[nodiscard]] double deceased_ratio(int i) const {
        return deceased[i] * population[i] / total_population();
    }

    [[nodiscard]] double susceptible_ratio() const {
        double res = 0;
        for (int i = 0; i < population.size(); i++) {
            res += susceptible_ratio(i);
        }
        return res;
    }

    [[nodiscard]] double infected_ratio() const {
        double res = 0;
        for (int i = 0; i < population.size(); i++) {
            res += infected_ratio(i);
        }
        return res;
    }

    [[nodiscard]] double recovered_ratio() const {
        double res = 0;
        for (int i = 0; i < population.size(); i++) {
            res += recovered_ratio(i);
        }
        return res;
    }

    [[nodiscard]] double deceased_ratio() const {
        double res = 0;
        for (int i = 0; i < population.size(); i++) {
            res += deceased_ratio(i);
        }
        return res;
    }
};

// Required for comparing states and detect any change
inline bool operator != (const sird &x, const sird &y) {
    return x.population != y.population || x.susceptible != y.susceptible || x.infected != y.infected ||
           x.recovered != y.recovered || x.deceased != y.deceased || x.mobility != y.mobility;
}

// Required for printing the state of the cell
std::ostream &operator << (std::ostream &os, const sird &x) {
    os << "<";

    for (int i = 0; i < x.population.size(); i++) {
        os << x.population[i] << "," << x.susceptible[i] << "," << x.infected[i] << "," << x.recovered[i] << "," << x.deceased[i] << ",";
    }

    os << x.total_population() << "," << x.susceptible_ratio() << "," << x.infected_ratio() << "," << x.recovered_ratio() << "," << x.deceased_ratio() << ">";
    return os;
}

// Required for creating SIR objects from JSON file
[[maybe_unused]] void from_json(const cadmium::json& j, sird &s) {
    j.at("population").get_to(s.population);
    j.at("susceptible").get_to(s.susceptible);
    j.at("infected").get_to(s.infected);
    j.at("recovered").get_to(s.recovered);
    j.at("deceased").get_to(s.deceased);

    int length = s.population.size();
    if (length != s.susceptible.size() || length != s.infected.size() || length != s.recovered.size() || length != s.deceased.size()) {
        throw "Vector sizes of JSON do not coincide.";
    }
    s.mobility = std::vector<double>(length, 1);  // By default, mobility is set to 1 for all age segments.
}

#endif //CISE_PANDEMIC_STATE_HPP
