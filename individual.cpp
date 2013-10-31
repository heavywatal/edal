// -*- mode: c++; coding: utf-8 -*-
#include "individual.h"

#include <cmath>
#include <iostream>

#include "cxxwtils/prandom.hpp"

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;
constexpr double Individual::STRENGTH_OF_MATING_PREFERENCE_;
constexpr size_t Individual::NUM_FITTEST_OFSPINRGS_;
constexpr size_t Individual::MAX_CARRYING_CAPACITY_;

double Individual::mate(const Individual& male) {
    if (phenotype(trait::choosiness) == 0.5) {
        return 1.0;
    }
    double exponent = phenotype(trait::female_trait);
    if (phenotype(trait::choosiness) > 0.5) {
        exponent -= male.phenotype(trait::male_trait);
    } else {
        exponent -= (1.0 - male.phenotype(trait::male_trait));
    }
    exponent *= exponent;
    double true_choosiness = (2 * phenotype(trait::choosiness) - 1);
    true_choosiness *= true_choosiness;
    exponent *= true_choosiness;
    exponent /= STRENGTH_OF_MATING_PREFERENCE_;
    exponent /= STRENGTH_OF_MATING_PREFERENCE_;
    exponent /= -2.0;
    return std::exp(exponent);
}

Individual::Loci Individual::recombination(const Loci& lhs, const Loci& rhs) {
    const Loci filter(prandom().randrange(pow<NUM_LOCI_>(2)));
    return (lhs & filter) | (rhs & ~filter);
}

std::vector<Individual::Loci> Individual::gametogenesis() const {
    std::vector<Loci> gamete;
    gamete.reserve(trait::size);
    for (size_t i=0; i<trait::size; ++i) {
        gamete.push_back(recombination(genotype_.first[i], genotype_.second[i]));
    }
    return gamete;
}

void individual_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Individual ind;
    std::cerr << ind.gametogenesis().at(0) << std::endl;
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
}
