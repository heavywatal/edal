// -*- mode: c++; coding: utf-8 -*-
/*! @file patch.cpp
    @brief Implementation of Patch class
*/
#include "patch.h"

#include <cmath>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/prandom.hpp"


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Patch::append(const Individual& ind) {
    if (prandom().bernoulli(0.5)) {
        females_.push_back(ind);
    } else {
        males_.push_back(ind);
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
            n += focal.habitat_overlap_roughgarden(ind) * focal.morphology_overlap_roughgarden(ind);
        }
    };
    impl(females_);
    impl(males_);
    return n;
}

void Patch::viability_selection() {
    auto impl = [this] (std::vector<Individual>* members) {
        std::vector<Individual> survivors;
        for (const auto& ind: *members) {
            if (prandom().bernoulli(ind.survival_probability(effective_num_competitors(ind)))) {
                survivors.push_back(ind);
            }
        }
        members->swap(survivors);
    };
    impl(&females_);
    impl(&males_);
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
    Patch patch(20);
    std::cerr << patch.size();
    for (size_t i=0; i<10; ++i) {
        for (const auto& child: patch.mate_and_reproduce()) {
            patch.append(child);
        }
        std::cerr << " b " << patch.size();
        patch.viability_selection();
        std::cerr << " d " << patch.size();
    }
    std::cerr << std::endl;
    std::cerr << patch << std::endl;
}