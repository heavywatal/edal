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

void Patch::change_sex_half(size_t n) {
    n /= 2;
    for (auto it=members_.begin()+n; it!=members_.end(); ++it) {
        it->change_sex();
    }
}

std::vector<Individual> Patch::mate_and_reproduce() const {
    std::poisson_distribution<size_t> poisson(Individual::AVG_NUM_OFFSPINRGS());
    std::vector<Individual> offsprings;
    std::vector<size_t> male_indices;
    male_indices.reserve(members_.size());
    for (size_t i=0; i<members_.size(); ++i) {
        if (members_[i].is_male()) {
            male_indices.push_back(i);
        }
    }
    if (male_indices.empty()) {return offsprings;}
    offsprings.reserve(members_.size() * Individual::AVG_NUM_OFFSPINRGS());
    for (const auto& mother: members_) {
        if (mother.is_male()) continue;
        std::vector<double> prefs;
        prefs.reserve(male_indices.size());
        for (const size_t i: male_indices) {
            prefs.push_back(mother.mating_probability(members_[i]));
        }
        std::vector<double> upper_bounds(male_indices.size());
        std::partial_sum(prefs.begin(), prefs.end(), upper_bounds.begin());
        std::uniform_real_distribution<> uniform(0.0, upper_bounds.back());
        const double dart = uniform(rng_);
        size_t father_i = 0;
        while (upper_bounds[father_i] < dart) {++father_i;}
        const Individual& father = members_[male_indices[father_i]];
        const size_t num_children = poisson(rng_);
        for (size_t i=0; i<num_children; ++i) {
            offsprings.emplace_back(
                mother.gametogenesis(rng_), father.gametogenesis(rng_),
                std::bernoulli_distribution(0.5)(rng_));
        }
    }
    return offsprings;
}

double Patch::effective_num_competitors(const Individual& focal) const {
    double n = 0;
    auto impl = [&] (const std::vector<Individual>& members) {
        for (const auto& ind: members) {
            n += focal.resource_overlap(ind);
        }
    };
    impl(members_);
    return n;
}

void Patch::viability_selection() {
    std::vector<size_t> indices;
    indices.reserve(members_.size());
    size_t i = 0;
    for (auto& ind: members_) {
        const double p = ind.survival_probability(effective_num_competitors(ind));
        if (std::bernoulli_distribution(p)(rng_)) {
            indices.push_back(i);
        }
        ++i;
    }
    std::vector<Individual> tmp;
    tmp.reserve(members_.size());
    for (auto i: indices) {
        tmp.push_back(std::move(members_[i]));
    }
    members_.swap(tmp);
}

std::pair<size_t, size_t> Patch::choose_patch(size_t row, size_t col) const {
    if (!std::bernoulli_distribution(Individual::MIGRATION_RATE())(rng_))
        {return {row, col};}
    switch (std::uniform_int_distribution<size_t>(0, 7)(rng_)) {
      case 0:        ++col; break;
      case 1: ++row; ++col; break;
      case 2: ++row;        break;
      case 3: ++row; --col; break;
      case 4:        --col; break;
      case 5: --row; --col; break;
      case 6: --row;        break;
      case 7: --row; ++col; break;
    }
    return {row, col};
}

std::map<Individual, size_t> Patch::summarize() const {
    std::map<Individual, size_t> genotypes;
    for (const auto& ind: members_) {
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