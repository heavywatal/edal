// -*- mode: c++; coding: utf-8 -*-
/*! @file patch.cpp
    @brief Implementation of Patch class
*/
#include "patch.h"

#include <cmath>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "cxxwtils/debug.hpp"
#include "cxxwtils/iostr.hpp"
#include "cxxwtils/prandom.hpp"


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Patch::append(Individual&& ind) {
    if (prandom().bernoulli(0.5)) {
        females_.push_back(std::move(ind));
    } else {
        males_.push_back(std::move(ind));
    }
}

std::vector<Individual> Patch::mate_and_reproduce() const {
    std::vector<Individual> offsprings;
    if (males_.empty()) {return offsprings;}
    offsprings.reserve(females_.size() * 8);
    for (const auto& mother: females_) {
        std::vector<double> prefs;
        prefs.reserve(males_.size());
        for (const auto& male: males_) {
            prefs.push_back(mother.mating_probability(male));
        }
        std::vector<double> upper_bounds(males_.size());
        std::partial_sum(prefs.begin(), prefs.end(), upper_bounds.begin());
        const double dart = prandom().uniform(upper_bounds.back());
        size_t father_i = 0;
        while (upper_bounds[father_i] < dart) {++father_i;}
        const Individual& father = males_[father_i];
        const size_t num_children = mother.poisson_offsprings();
        for (size_t i=0; i<num_children; ++i) {
            offsprings.push_back(Individual(mother.gametogenesis(), father.gametogenesis()));
        }
    }
    return offsprings;
}

double Patch::effective_num_competitors(const Individual& focal) const {
    double n = 0;
    auto impl = [&] (const std::vector<Individual>& members) {
        for (const auto& ind: members) {
//            n += focal.resource_overlap(ind); // too slow
            n += focal.preference_overlap(ind) * focal.morphology_overlap(ind);
        }
    };
    impl(females_);
    impl(males_);
    return n;
}

void Patch::viability_selection() {
    auto choose = [this] (const std::vector<Individual>& members) {
        std::vector<size_t> indices;
        indices.reserve(members.size());  // TODO
        size_t i = 0;
        for (auto& ind: members) {
            if (prandom().bernoulli(ind.survival_probability(effective_num_competitors(ind)))) {
                indices.push_back(i);
            }
            ++i;
        }
        return indices;
    };
    const auto indices_female = choose(females_);
    const auto indices_male = choose(males_);
    std::vector<Individual> tmp;
    tmp.reserve(std::max(females_.capacity(), males_.capacity()));
    // to prevent re-allocation, not indices.size()
    auto extract = [&tmp](
        const std::vector<size_t>& indices, std::vector<Individual>* members) {
        for (auto i: indices) {
            tmp.push_back(std::move(members->operator[](i)));
        }
        members->swap(tmp);
        tmp.clear();
    };
    extract(indices_female, &females_);
    extract(indices_male, &males_);
}

std::map<Individual, size_t> Patch::summarize() const {
    std::map<Individual, size_t> genotypes;
    for (const auto& ind: females_) {
        ++genotypes[ind];
    }
    for (const auto& ind: males_) {
        ++genotypes[ind];
    }
    return genotypes;
}

std::ostream& operator<< (std::ostream& ost, const Patch& patch) {
    return ost << patch.summarize();
    // operator<< for std::map is defined in "iostr.hpp"
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Patch::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Patch patch(20, Individual({15,0,15,0}));
    std::cerr << patch.size();
    for (size_t i=0; i<10; ++i) {
        for (auto& child: patch.mate_and_reproduce()) {
            patch.append(std::move(child));
        }
        std::cerr << " b " << patch.size();
        patch.viability_selection();
        std::cerr << " d " << patch.size();
    }
    std::cerr << std::endl;
    std::cerr << patch << std::endl;
}