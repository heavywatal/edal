// -*- mode: c++; coding: utf-8 -*-
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/math/distributions.hpp>
#include <boost/program_options.hpp>

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/prandom.hpp"

size_t Individual::CARRYING_CAPACITY = 1000;
size_t Individual::AVG_NUM_OFFSPINRGS_ = 4;
double Individual::HEIGHT_PREFERENCE_ = 0.5;
double Individual::DIAMETER_PREFERENCE_ = 0.5;
double Individual::MATING_SIGMA_ = 0.2;
double Individual::TOEPAD_SELECTION_ = 0.5;
double Individual::LIMB_SELECTION_ = 0.5;
double Individual::HEIGHT_COMPETITION_ = 0.5;
double Individual::DIAMETER_COMPETITION_ = 0.5;
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
        ("height_pref", po::value<double>(&HEIGHT_PREFERENCE_)->default_value(HEIGHT_PREFERENCE_))
        ("diameter_pref", po::value<double>(&DIAMETER_PREFERENCE_)->default_value(DIAMETER_PREFERENCE_))
        ("mating_s,s", po::value<double>(&MATING_SIGMA_)->default_value(MATING_SIGMA_))
        ("toepad_s,t", po::value<double>(&TOEPAD_SELECTION_)->default_value(TOEPAD_SELECTION_))
        ("limb_s,l", po::value<double>(&LIMB_SELECTION_)->default_value(LIMB_SELECTION_))
        ("mu_locus,u", po::value<double>(&MU_LOCUS_)->default_value(MU_LOCUS_))
        ("mu_neutral,n", po::value<double>(&MU_NEUTRAL_)->default_value(MU_NEUTRAL_))
    ;
    return desc;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

template <class Func> inline
double integral(Func func) {
    constexpr int precision = 50;
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

inline double abundance_v2(const double height, const double diameter) {
    constexpr double mean_height = 0.5;
    constexpr double sd_height = 0.4;
    constexpr double c0 = 2.0;
    constexpr double c1 = c0 - 1.0;
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

constexpr double height_alpha = 2.0;

inline double pdf_beta(const double height, const double diameter) {
    static_cast<void>(diameter);
    double result = std::pow(height * (1 - height), height_alpha - 1);
    result *= std::pow(4, height_alpha - 1);
    return result;
}

inline double pdf_triangle(const double height, const double diameter) {
    if (height == 1.0 || 1.0 - height < diameter) {return 0.0;}
    double result = (1.0 - height - diameter);
    result /= (1.0 - height);
    return result;
}

inline double abundance_v3(const double height, const double diameter) {
    return pdf_beta(height, diameter) * pdf_triangle(height, diameter);
}

inline double abundance(const double height, const double diameter) {
    return abundance_v3(height, diameter);
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// slow...
//inline double pdf_beta_boost(const double height, const double diameter) {
//    static_cast<void>(diameter);
//    static boost::math::beta_distribution<> beta(height_alpha, height_alpha);
//    static const double normalizer = 1.0 / boost::math::pdf(beta, 0.5);
//    double p = boost::math::pdf(beta, height);
//    return p *= normalizer;
//}
//
//inline double pdf_triangle_boost(const double height, const double diameter) {
//    if (height == 1.0) {return 0.0;}
//    boost::math::triangular_distribution<> tri(0, 0, 1 - height);
//    const double normalizer = (1.0 - height) / 2.0;
//    double p = boost::math::pdf(tri, diameter);
//    return p *= normalizer;
//}
//
//inline double abundance_v3_boost(const double height, const double diameter) {
//    return pdf_beta_boost(height, diameter) * pdf_triangle_boost(height, diameter);
//}
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

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
//    sqrt_denominator_2_ = sqrt_denom_2_();
    effective_carrying_capacity_ = effective_carrying_capacity();
}

double Individual::denom_numerical() const {
    return integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
}

double Individual::denom_mathematica() const {
    const double a = height_alpha;
    const double h0 = HEIGHT_PREFERENCE_;
    const double h1 = DIAMETER_PREFERENCE_;
    const double y0 = phenotype_[trait::height_preference];
    const double y1 = phenotype_[trait::diameter_preference];
    double z = -12;
    z += 6 * h0;
    z += h1;
    z += a * (-24 + h1 + 6 * h0 * wtl::pow<2>(1 - 2 * y0)
              - 8 * h1 * y1 + 24 * h1 * wtl::pow<2>(y1));
    z += 4 * (3 * h0 * (-1 + y0) * y0 + h1 * y1 * (-1 + 3 * y1));
    z *= std::sqrt(M_PI);
    z *= std::pow(2, -2 * (a + 1));
    z /= -3;
    z /= std::tgamma(a + 3 / 2);
    //z *= std::tgamma(a);
    z *= std::tgamma(2 * a);
    z /= std::tgamma(a);
    return z;
}

double Individual::denom_maple() const {
    const double alpha = height_alpha;
    const double h0 = HEIGHT_PREFERENCE_;
    const double h1 = DIAMETER_PREFERENCE_;
    const double y0 = phenotype_[trait::height_preference];
    const double y1 = phenotype_[trait::diameter_preference];
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
    z *= std::tgamma(2 * alpha);
    z /= wtl::pow<2>(std::tgamma(alpha));
    return z /= 12;  // TODO: wrong sign?
}

double Individual::sqrt_denom_2_() const {
    return std::sqrt(integral([this](const double height, const double diameter)->double {
        double result = habitat_preference(height, diameter);
        result *= result;
        return result *= abundance(height, diameter);
    }));
}


double Individual::habitat_preference(const double height, const double diameter) const {
    auto impl = [](const double ind_preference,
                   const double env_characterstics,
                   const double habitat_pref_strength) {
        double exponent = ind_preference;
        exponent -= env_characterstics;
        exponent *= exponent;
        return exponent *= -habitat_pref_strength;
    };
    double exponent = impl(phenotype_[trait::height_preference], height, HEIGHT_PREFERENCE_);
    exponent += impl(phenotype_[trait::diameter_preference], diameter, DIAMETER_PREFERENCE_);
    return std::exp(exponent);
}

double Individual::habitat_overlap_v2(const Individual& other) const {
    double n = integral([this, &other](const double height, const double diameter) {
        double result = habitat_preference(height, diameter);
        result *= other.habitat_preference(height, diameter);
        return result *= abundance(height, diameter);
    });
//    n /= sqrt_denominator_2_;
//    n /= other.sqrt_denominator_2_;
    return n;
}

double Individual::habitat_overlap_v3(const Individual& other) const {
    auto impl = [](const double yi, const double yj, const double c) {
        double exponent = yi;
        exponent -= yj;
        exponent *= exponent;
        return exponent *= -c;
    };
    double exponent = impl(phenotype_[trait::height_preference],
                           other.phenotype_[trait::height_preference],
                           HEIGHT_COMPETITION_);
    exponent += impl(phenotype_[trait::diameter_preference],
                     other.phenotype_[trait::diameter_preference],
                     DIAMETER_COMPETITION_);
    return std::exp(exponent);
}

double Individual::fitness(const double height, const double diameter) const {
    auto impl = [](const double x, const double mu, const double selection) {
        double exponent = x;
        exponent -= mu;
        exponent *= exponent;
        return exponent *= -selection;
    };
    double exponent = impl(phenotype_[trait::toepad_size], height, TOEPAD_SELECTION_);
    exponent += impl(phenotype_[trait::limb_length], diameter, LIMB_SELECTION_);
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

bool Individual::survive(const double effective_num_competitors) const {
    double denom = AVG_NUM_OFFSPINRGS_;
    denom -= 1.0;
    denom *= effective_num_competitors;  // ^ alpha for crowding strength
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
    return std::exp(exponent) * habitat_overlap(male);
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

template <class Func> inline
std::string resource_abundance_test(Func func) {
    constexpr int precision = 25;
    constexpr double delta = 1.0 / precision;
    std::ostringstream ost;
    std::string sep(",");
    ost << "height" << sep << "diameter" << sep << "abundance\n";
    for (size_t i=0; i<=precision; ++i) {
        for (size_t j=0; j<=precision; ++j) {
            const double x = i * delta;
            const double y = j * delta;
            ost << x << sep << y << sep << func(x, y) << "\n";
        }
    }
    return ost.str();
}

void individual_unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Individual ind;
    std::cerr << ind.gametogenesis().at(0) << std::endl;
    std::cerr << ind.effective_carrying_capacity() << std::endl;
    std::cerr << ind.denom_numerical() << std::endl;
    std::cerr << ind.denom_maple() << std::endl;
    std::cerr << ind.denom_mathematica() << std::endl;
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
    wtl::Fout{"ignore/abundance_v2.csv"} << resource_abundance_test(abundance_v2);
    wtl::Fout{"ignore/abundance_beta.csv"} << resource_abundance_test(pdf_beta);
    wtl::Fout{"ignore/abundance_triangle.csv"} << resource_abundance_test(pdf_triangle);
    wtl::Fout{"ignore/abundance_v3.csv"} << resource_abundance_test(abundance_v3);
}
