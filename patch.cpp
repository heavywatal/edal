// -*- mode: c++; coding: utf-8 -*-
#include "patch.h"

#include <cmath>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

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
    offsprings.reserve(females_.size() * Individual::AVG_NUM_OFFSPINRGS_ + 8);
    for (const auto& mother: females_) {
        std::vector<double> prefs;
        prefs.reserve(males_.size());
        for (const auto& male: males_) {
            prefs.push_back(mother.mating_preference(male));
        }
        std::vector<double> upper_bounds(males_.size());
        std::partial_sum(prefs.begin(), prefs.end(), upper_bounds.begin());
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

double Patch::effective_num_competitors_f(const size_t index) const {
    double n = 0;
    for (size_t i=0; i<females_.size(); ++i) {
        if (i==index) continue;
        n += females_[index].effective_num_competitors(females_[i]);
    }
    for (const auto& ind: males_) {
        n += females_[index].effective_num_competitors(ind);
    }
    return n;
}

double Patch::effective_num_competitors_m(const size_t index) const {
    double n = 0;
    for (const auto& ind: females_) {
        n += males_[index].effective_num_competitors(ind);
    }
    for (size_t i=0; i<males_.size(); ++i) {
        if (i==index) continue;
        n += males_[index].effective_num_competitors(males_[i]);
    }
    return n;
}

void Patch::viability_selection() {
    std::vector<Individual> survivors;
    for (size_t i=0; i<females_.size(); ++i) {
        if (females_[i].survive(effective_num_competitors_f(i))) {
            survivors.push_back(females_[i]);
        }
    }
    females_.swap(survivors);
    survivors.clear();
    for (size_t i=0; i<males_.size(); ++i) {
        if (males_[i].survive(effective_num_competitors_m(i))) {
            survivors.push_back(males_[i]);
        }
    }
    males_.swap(survivors);
}

std::string Patch::str() const {
    std::ostringstream ost;
    for (const auto& ind: females_) {
        ost << ind.str() << "\n";
    }
    for (const auto& ind: males_) {
        ost << ind.str() << "\n";
    }
    return ost.str();
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