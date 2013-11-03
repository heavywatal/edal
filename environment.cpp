// -*- mode: c++; coding: utf-8 -*-
#include "environment.h"

#include <cmath>

#include "cxxwtils/prandom.hpp"


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {
} // namespace
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
    offsprings.reserve(females_.size() * Individual::AVG_NUM_OFFSPINRGS_ + 8);
    for (const auto& mother: females_) {
        std::vector<double> freqs;
        freqs.reserve(males_.size());
        for (const auto& male: males_) {
            freqs.push_back(mother.mating_frequencies(male));
        }
        std::vector<double> upper_bounds(males_.size());
        std::partial_sum(freqs.begin(), freqs.end(), upper_bounds.begin());
        const double dart = prandom().uniform(upper_bounds.back());
        size_t father_i = 0;
        while (upper_bounds[father_i] < dart) {++father_i;}
        const Individual& father = males_[father_i];
        const size_t num_children = prandom().poisson(Individual::AVG_NUM_OFFSPINRGS_);
        for (size_t i=0; i<num_children; ++i) {
            offsprings.push_back(Individual(mother.gametogenesis(), father.gametogenesis()));
        }
    }
    return offsprings;
}

double Patch::effective_population_size() const {
    double n = 0;
    for (const auto& ind: females_) {
        n += ind.effective_population_size_i();
    }
    for (const auto& ind: males_) {
        n += ind.effective_population_size_i();
    }
    return n;
}

void Patch::viability_selection() {
    const double N_e = effective_population_size();
    std::vector<Individual> survivors;
    for (const auto& ind: females_) {
        if (ind.survive(N_e)) {
            survivors.push_back(ind);
        }
    }
    females_.swap(survivors);
    survivors.clear();
    for (const auto& ind: males_) {
        if (ind.survive(N_e)) {
            survivors.push_back(ind);
        }
    }
    males_.swap(survivors);
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void patch_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Patch patch(20);
    std::cerr << patch.size() << std::endl;
    for (size_t i=0; i<10; ++i) {
        for (const auto& child: patch.mate_and_reproduce()) {
            patch.append(child);
        }
        std::cerr << patch.size() << std::endl;
        patch.viability_selection();
        std::cerr << patch.size() << std::endl;
    }
}