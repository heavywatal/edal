// -*- mode: c++; coding: utf-8 -*-
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/prandom.hpp"

size_t Individual::CARRYING_CAPACITY = 60;
size_t Individual::AVG_NUM_OFFSPINRGS_ = 4;
double Individual::HABITAT_SIGMA_ = 0.2;
double Individual::MATING_SIGMA_ = 0.2;
double Individual::ADAPTIVE_T_SIGMA_ = 1.0;
double Individual::ADAPTIVE_L_SIGMA_ = 1.0;
double Individual::MU_LOCUS_ = 1e-4;
double Individual::MU_NEUTRAL_ = 1e-4;

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;

boost::program_options::options_description& Individual::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Individual"};
    desc.add_options()
        ("birth_rate,b", po::value<size_t>(&AVG_NUM_OFFSPINRGS_)->default_value(AVG_NUM_OFFSPINRGS_))
        ("carrying_capacity,K", po::value<size_t>(&CARRYING_CAPACITY)->default_value(CARRYING_CAPACITY))
        ("habitat_s,p", po::value<double>(&HABITAT_SIGMA_)->default_value(HABITAT_SIGMA_))
        ("mating_s,s", po::value<double>(&MATING_SIGMA_)->default_value(MATING_SIGMA_))
        ("adaptive_t_s,t", po::value<double>(&ADAPTIVE_T_SIGMA_)->default_value(ADAPTIVE_T_SIGMA_))
        ("adaptive_l_s,l", po::value<double>(&ADAPTIVE_L_SIGMA_)->default_value(ADAPTIVE_L_SIGMA_))
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
    h_exponent *= -0.5;
    double theta_u = 1.0;
    theta_u /= (c0 - c1 * height);
    double d_exponent = -diameter;
    d_exponent *= theta_u;
    double result = k;
    result *= std::exp(h_exponent + d_exponent);
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
    auto impl = [](const double ind_preference,
                   const double env_characterstics) {
        double exponent = ind_preference;
        exponent -= env_characterstics;
        exponent /= HABITAT_SIGMA_;
        exponent *= exponent;
        return exponent *= -0.5;
    };
    double exponent = impl(phenotype_[trait::height_preference], height);
    exponent += impl(phenotype_[trait::diameter_preference], diameter);
    return std::exp(exponent);
}

double Individual::denom_() const {
    return integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
}

double Individual::denom_new() const {
    constexpr size_t alpha = 1.0;
    constexpr size_t h0 = 1.0;
    constexpr size_t h1 = 1.0;
    const double y0 = phenotype_[trait::toepad_size];
    const double y1 = phenotype_[trait::limb_length];
    const double minus1_2a = std::pow(-1, 2 * alpha);
    const double minus1_2a_1 = std::pow(-1, 2 * alpha + 1);
    const double alpha_2 =  wtl::pow<2>(alpha);
    const double alpha_3 =  wtl::pow<3>(alpha);
    double z = 0.0;
    z +=  8 * minus1_2a_1 * alpha_3 * h1 * y1;
    z += 60 * minus1_2a_1 * alpha_2 * h1 * wtl::pow<2>(y1);
    z += 60 * minus1_2a_1 * alpha_2 * h0 * wtl::pow<2>(y0);
    z += 24 * minus1_2a   * alpha_3 * h0 * wtl::pow<2>(y0);
    z += 12 * minus1_2a   * alpha   * h1 * wtl::pow<2>(y1);
    z += 12 * minus1_2a   * alpha   * h0 * wtl::pow<2>(y0);
    z += 24 * minus1_2a   * alpha_3 * h1 * wtl::pow<2>(y1);
    z += 12 * minus1_2a_1 * alpha;
    z += 24 * minus1_2a_1 * alpha_3;
    z += 12 * minus1_2a   *           h0;
    z +=  2 * minus1_2a   *           h1;
    z += 60 * minus1_2a   * alpha_2;
    z +=  8 * minus1_2a_1 *           h1 * y1;
    z += 12 * minus1_2a_1 * alpha_2 * h0;
    z +=  6 * minus1_2a_1 * alpha   * h0;
    z += 24 * minus1_2a   *           h1 * wtl::pow<2>(y1);
    z +=      minus1_2a   * alpha_3 * h1;
    z -=      minus1_2a   * alpha   * h1;
    z +=  6 * minus1_2a   * alpha_3 * h0;
    z += 24 * minus1_2a_1 *           h0 * y0;
    z +=  2 * minus1_2a_1 * alpha_2 * h1;
    z += 24 * minus1_2a   *           h0 * wtl::pow<2>(y0);
    z += 24 * minus1_2a_1;
    z += 24 * minus1_2a_1 * alpha_3 * h0 * y0;
    z += 60 * minus1_2a   * alpha_2 * h0 * y0;
    z += 20 * minus1_2a   * alpha_2 * h1 * y1;
    z +=  4 * minus1_2a_1 * alpha   * h1 * y1;
    z += 12 * minus1_2a_1 * alpha   * h0 * y0;
    z *= std::pow(4, -alpha);
    z *= std::sqrt(M_PI);
    z *= std::tgamma(alpha - 2);
    z /= std::tgamma(alpha + 3 / 2);
    return z /= -12;
}

double Individual::sqrt_denom_2_() const {
    return std::sqrt(integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        result *= result;
        return result *= abundance(height, diameter);
    }));
}

double Individual::fitness(const double height, const double diameter) const {
    auto impl = [](const double x, const double mu, const double sigma) {
        double exponent = x;
        exponent -= mu;
        exponent /= sigma;
        exponent *= exponent;
        return exponent;
    };
    double exponent = impl(phenotype_[trait::toepad_size], height, ADAPTIVE_T_SIGMA_);
    exponent += impl(phenotype_[trait::limb_length], diameter, ADAPTIVE_L_SIGMA_);
    exponent *= -0.5;
    return std::exp(exponent);
}

double Individual::effective_carrying_capacity() const {
    double result = CARRYING_CAPACITY;
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
    denom *= effective_population_size;  // ^ alpha for crowding strength
    denom /= effective_carrying_capacity_;
    denom += 1.0;
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
    double true_choosiness = choosiness;
    true_choosiness *= 2.0;
    true_choosiness -= 1.0;
    exponent *= true_choosiness;
    exponent /= MATING_SIGMA_;
    exponent *= exponent;
    exponent *= -0.5;
    return std::exp(exponent);
}

size_t Individual::poisson_offsprings() const {
    return prandom().poisson(AVG_NUM_OFFSPINRGS_);
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
