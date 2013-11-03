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
constexpr size_t Individual::AVG_NUM_OFFSPINRGS_;
constexpr double Individual::MU_LOCUS_;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

template <class Func> inline
double integral(Func func) {
    constexpr int precision = 10;
    constexpr double delta = 1.0 / precision;
    double result = 0.0;
    for (size_t i=0; i<=precision; ++i) {
        for (size_t j=0; j<=precision; ++j) {
            result += func(i * delta, j * delta);
        }
    }
    result *= delta;
    result *= delta;
    return result;
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

inline double preference_impl(const double habitat_preference,
                              const double preference_strength,
                              const double env_characterstic) {
    constexpr double sigma_p = 1.0;

    double exponent = habitat_preference;
    exponent -= env_characterstic;
    exponent *= preference_strength;
    exponent /= sigma_p;
    exponent *= exponent;
    exponent /= -2.0;
    return std::exp(exponent);
}

inline double carrying_capacity_coef() {
    constexpr size_t max_carrying_capacity = 30;

    std::vector<size_t> trait_values(trait::size * 2 - 1);
    std::iota(trait_values.begin(), trait_values.end(), 0);
    double kii = 0;
    double kii_prime = 0;
    for (const auto x0: trait_values) {
     for (const auto x1: trait_values) {
      for (const auto x2: trait_values) {
       for (const auto x3: trait_values) {
        Individual ind(std::vector<size_t>{x0, x1, x2, x3});
        //std::cout << ind.carrying_capacity() << " ";
        //std::cout << ind.carrying_capacity_prime() << std::endl;
        kii += ind.carrying_capacity();
        kii_prime += ind.carrying_capacity_prime();
    }}}}
    std::cout << "K(I, I):  " << kii << std::endl;
    std::cout << "K(I, I'): " << kii_prime << std::endl;
    return max_carrying_capacity * kii / kii_prime;
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

double Individual::habitat_preference(const double height, const double diameter) const {
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
    return integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
}

double Individual::denominator_2() const {
    return integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        result *= result;
        return result *= abundance(height, diameter);
    });
}

double Individual::denominator_prime() const {
    static const double dprime = integral([this](const double height, const double diameter)->double {
        return abundance(height, diameter);
    });
    return dprime;
}

double Individual::denominator_prime_2() const {
    return denominator_prime();
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
        result *= habitat_preference(height, diameter);
        result *= abundance(height, diameter);
        return result;
    }) / denominator();
}

double Individual::carrying_capacity_prime() const {
    return integral([this](const double height, const double diameter)->double {
        double result = fitness(height, diameter);
        result *= habitat_preference(height, diameter);
        result *= abundance(height, diameter);
        return result;
    }) / denominator_prime();
}

double Individual::effective_population_size_i() const {
    double n = 1.0;
    n /= std::sqrt(denominator_2());
    n /= std::sqrt(denominator_prime_2());
    n *= integral([this](const double height, const double diameter) {
        double result = habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
    return n;
}

bool Individual::survive(const double effective_population_size) const {
    static const double ccc = carrying_capacity_coef();
    constexpr double alpha = 1.0;
    double denom = AVG_NUM_OFFSPINRGS_;
    denom -= 1.0;
    denom *= std::pow(effective_population_size, alpha);
    denom /= ccc;
    denom /= carrying_capacity_prime();
    return prandom().bernoulli(1.0 / denom);
}

double Individual::mating_preference(const Individual& male) const {
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

double Individual::mating_frequencies(const Individual& male) const {
    double probability = mating_preference(male);
    probability /= std::sqrt(denominator_2());
    probability /= std::sqrt(denominator_prime_2());
    probability *= integral(
        [this, male](const double height, const double diameter) {
        double result = habitat_preference(height, diameter);
        return result *= male.habitat_preference(height, diameter);
    });
    return probability;
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
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void individual_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Individual ind;
    std::cerr << ind.gametogenesis().at(0) << std::endl;
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
}
