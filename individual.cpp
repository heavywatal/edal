// -*- mode: c++; coding: utf-8 -*-
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "cxxwtils/prandom.hpp"

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;
constexpr double Individual::STRENGTH_OF_MATING_PREFERENCE_;
constexpr size_t Individual::NUM_FITTEST_OFSPINRGS_;
constexpr size_t Individual::MAX_CARRYING_CAPACITY_;
constexpr double Individual::MU_LOCUS_;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {
inline double preference_impl(const double habitat_preference,
                              const double preference_strength,
                              const double env_characterstic) {
    double exponent = habitat_preference;
    exponent -= env_characterstic;
    exponent *= preference_strength;
    //exponent /= sigma_p;
    exponent *= exponent;
    exponent /= -2.0;
    return std::exp(exponent);
}

inline double abundance(const double height, const double diameter) {
    constexpr double mean_height = 0.5;
    constexpr double sd_height = 0.5;
    constexpr double c0 = 0.5;
    constexpr double c1 = c0 - 0.01;
    constexpr double k = 1.0;
    double h_exponent = height;
    h_exponent -= mean_height;
    h_exponent /= sd_height;
    h_exponent *= h_exponent;
    h_exponent *= -2;
    double theta_u = 1.0;
    theta_u /= (c0 - c1 * height);
    double result = k;
    result *= std::exp(h_exponent);
    result *= std::exp(theta_u * diameter);
    result *= theta_u;
    return result;
}

} // namespace
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Individual::Individual(const std::vector<size_t>& values): genotype_{{}, {}} {
    genotype_.first.reserve(trait::size);
    genotype_.second.reserve(trait::size);
    for (auto x: values) {
        if (x > NUM_LOCI_) {
            genotype_.first.push_back(FULL_BITS);
            genotype_.second.push_back(std::pow(2, x - NUM_LOCI_) - 1);
        } else {
            genotype_.first.push_back(std::pow(2, x) - 1);
            genotype_.second.push_back(0);
        }
    }
    for (size_t i=0; i<2; ++i) {
        genotype_.first.push_back(FULL_BITS);
        genotype_.second.push_back(FULL_BITS);
    }
    for (size_t i=0; i<3; ++i) {
        genotype_.first.push_back(HALF_BITS);
        genotype_.second.push_back(HALF_BITS);
    }
}

double Individual::preference(const double height, const double diameter) const {
    double result = 1.0;
    result *= preference_impl(phenotype(trait::height_preference),
                              phenotype(trait::height_tolerance),
                              height);
    result *= preference_impl(phenotype(trait::diameter_preference),
                              phenotype(trait::diameter_tolerance),
                              diameter);
    return result;
}



double Individual::denominator() const {
    return integral([this](const double x, const double y)->double {
        return preference(x, y) * abundance(x, y);
    });
}

double Individual::denominator_prime() const {
    return integral([this](const double x, const double y)->double {
        return abundance(x, y);
    });
}

double Individual::denominator_2() const {
    return integral([this](const double x, const double y)->double {
        double result = preference(x, y);
        result *= result;
        return result *= abundance(x, y);
    });
}

double Individual::denominator_prime_2() const {
    return integral([this](const double x, const double y)->double {
        return abundance(x, y);
    });
}

double Individual::fitness(const double height, const double diameter) const {
    constexpr double sigma_h = 1.0;
    constexpr double sigma_d = 1.0;
    double lhs = phenotype(trait::toepad_size);
    lhs -= height;
    lhs /= sigma_h;
    lhs *= lhs;
    double rhs = phenotype(trait::limb_length);
    rhs -= diameter;
    rhs /= sigma_d;
    rhs *= rhs;
    double exponent = lhs + rhs;
    exponent *= -0.5;
    return std::exp(exponent);
}

double Individual::carrying_capacity() const {
    return integral([this](const double height, const double diameter)->double {
        double result = fitness(height, diameter);
        result *= preference(height, diameter);
        result *= abundance(height, diameter);
        return result;
    }) / denominator();
}

double Individual::carrying_capacity_prime() const {
    return integral([this](const double height, const double diameter)->double {
        double result = fitness(height, diameter);
        result *= preference(height, diameter);
        result *= abundance(height, diameter);
        return result;
    }) / denominator_prime();
}

bool Individual::survive(const double effective_carrying_capacity,
                         const double effective_population_size) const {
    constexpr double alpha = 1.0;
    double denom = NUM_FITTEST_OFSPINRGS_;
    denom -= 1.0;
    denom *= std::pow(effective_population_size, alpha);
    denom /= effective_carrying_capacity;
    return prandom().bernoulli(1.0 / denom);
}

double Individual::mate(const Individual& male) {
    const double choosiness = phenotype(trait::choosiness);
    if (choosiness == 0.5) {
        // random mating
        return 1.0;
    }
    double exponent = phenotype(trait::female_trait);
    if (choosiness > 0.5) {
        // assortative mating
        exponent -= male.phenotype(trait::male_trait);
    } else {
        // disassortative mating
        exponent -= 1.0;
        exponent += male.phenotype(trait::male_trait);
    }
    exponent *= exponent;
    double true_choosiness = (2 * choosiness - 1);
    true_choosiness *= true_choosiness;
    exponent *= true_choosiness;
    exponent /= STRENGTH_OF_MATING_PREFERENCE_;
    exponent /= STRENGTH_OF_MATING_PREFERENCE_;
    exponent /= -2.0;
    return std::exp(exponent);
}

std::vector<Individual::Loci> Individual::mutate(std::vector<Individual::Loci> haplotype) {
    const size_t num_mutations = prandom().poisson(MU_LOCUS_ * NUM_LOCI_ * trait::size);
    std::vector<size_t> indices(NUM_LOCI_ * trait::size);
    std::iota(indices.begin(), indices.end(), 0);
    indices = prandom().sample(std::move(indices), num_mutations);
    for (const auto& i: indices) {
        haplotype[i / NUM_LOCI_].flip(i % NUM_LOCI_);
    }
    return haplotype;
}

Individual::Loci Individual::recombination(const Loci& lhs, const Loci& rhs) {
    const Loci filter(prandom().randrange(wtl::pow<NUM_LOCI_>(2)));
    return (lhs & filter) | (rhs & ~filter);
}

std::vector<Individual::Loci> Individual::gametogenesis() const {
    std::vector<Loci> gamete;
    gamete.reserve(trait::size);
    for (size_t i=0; i<trait::size; ++i) {
        gamete.push_back(recombination(genotype_.first[i], genotype_.second[i]));
    }
    return mutate(gamete);
}


std::string Individual::str() const {
    std::ostringstream ost;
    for (const auto& loci: genotype_.first) {
        ost << loci << " ";
    }
    ost << "\n";
    for (const auto& loci: genotype_.second) {
        ost << loci << " ";
    }
    ost << "\n";
    return ost.str();
//    return str_haplotype(genotype_.first) + "\n" + str_haplotype(genotype_.second) + "\n";
}

void individual_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Individual ind;
    std::cerr << ind.gametogenesis().at(0) << std::endl;
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
}
