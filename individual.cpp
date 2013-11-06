// -*- mode: c++; coding: utf-8 -*-
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/prandom.hpp"

size_t Individual::AVG_NUM_OFFSPINRGS_ = 4;
double Individual::STRENGTH_OF_MATING_PREFERENCE_ = 0.05;
double Individual::MU_LOCUS_ = 1e-5;
double Individual::MU_NEUTRAL_ = 1e-5;

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;

boost::program_options::options_description& Individual::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Individual"};
    desc.add_options()
        ("mu_locus,u", po::value<double>(&MU_LOCUS_)->default_value(MU_LOCUS_))
        ("mu_neutral,n", po::value<double>(&MU_NEUTRAL_)->default_value(MU_NEUTRAL_))
    ;
    return desc;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

template <class Func> inline
double integral(Func func) {
    constexpr int precision = 20;
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
    constexpr double c0 = 2.0;
    constexpr double c1 = 1.0;
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
                              const double env_characterstic) {
    constexpr double sigma_p = 0.2;

    double exponent = habitat_preference;
    exponent -= env_characterstic;
    exponent /= sigma_p;
    exponent *= exponent;
    exponent /= -2.0;
    return std::exp(exponent);
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
    for (size_t i=0; i<4; ++i) {
        genotype_.first.push_back(HALF_BITS);
        genotype_.second.push_back(HALF_BITS);
    }
    phenotype_ = init_phenotype();
    denominator_ = denom_();
    sqrt_denominator_2_ = sqrt_denom_2_();
    effective_carrying_capacity_ = effective_carrying_capacity();
}

double Individual::habitat_preference(const double height, const double diameter) const {
    double result = preference_impl(phenotype_[trait::height_preference], height);
    result *= preference_impl(phenotype_[trait::diameter_preference], diameter);
    return result;
}


double Individual::denom_() const {
    return integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
}

double Individual::sqrt_denom_2_() const {
    return std::sqrt(integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        result *= result;
        return result *= abundance(height, diameter);
    }));
}

double Individual::fitness(const double height, const double diameter) const {
    constexpr double sigma_h = 1.0;
    constexpr double sigma_d = 1.0;
    double lhs = phenotype_[trait::toepad_size];
    lhs -= height;
    lhs /= sigma_h;
    lhs *= lhs;
    double rhs = phenotype_[trait::limb_length];
    rhs -= diameter;
    rhs /= sigma_d;
    rhs *= rhs;
    double exponent = lhs + rhs;
    exponent *= -0.5;
    return std::exp(exponent);
}

double Individual::effective_carrying_capacity() const {
    constexpr size_t max_carrying_capacity = 30;
    double result = max_carrying_capacity;
    result /= denominator_;
    result *= integral([this](const double height, const double diameter)->double {
        double result = fitness(height, diameter);
        result *= habitat_preference(height, diameter);
        result *= abundance(height, diameter);
        return result;
    });
    return result;
}

double Individual::effective_num_competitors(const Individual& opponent) const {
    double n = integral([this, &opponent](const double height, const double diameter) {
        double result = habitat_preference(height, diameter);
        result *= opponent.habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
    n /= sqrt_denominator_2_;
    n /= opponent.sqrt_denominator_2_;
    return n;
}

bool Individual::survive(const double effective_population_size) const {
    double denom = AVG_NUM_OFFSPINRGS_;
    denom -= 1.0;
    denom *= effective_population_size;
    denom /= effective_carrying_capacity_;
    return prandom().bernoulli(1.0 / denom);
}

double Individual::mating_preference(const Individual& male) const {
    const double choosiness = phenotype_[trait::choosiness];
    if (choosiness == 0.5) {
        // random mating
        return 1.0;
    }
    double exponent = phenotype_[trait::female_trait];
    if (choosiness > 0.5) {
        // assortative mating
        exponent -= male.phenotype_[trait::male_trait];
    } else {
        // disassortative mating
        exponent -= 1.0;
        exponent += male.phenotype_[trait::male_trait];
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
    std::string sep{","};
    ost << wtl::str_join(genotype_.first, sep) << sep;
    ost << wtl::str_join(genotype_.second, sep) << sep;
    ost << wtl::str_join(phenotype_, sep);
    return ost.str();
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void individual_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Individual ind;
    std::cerr << ind.gametogenesis().at(0) << std::endl;
    std::cerr << ind.effective_carrying_capacity() << std::endl;
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
}
