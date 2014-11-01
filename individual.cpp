// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/gz.hpp"

double Individual::BETA_PARAM_ = 3.0;
double Individual::NORMAL_SIGMA_ = 0.1;
double Individual::C0_ = 0.5;
double Individual::C1_ = 0.5;
size_t Individual::CARRYING_CAPACITY_ = 160;
size_t Individual::AVG_NUM_OFFSPINRGS_ = 4;
double Individual::HEIGHT_PREFERENCE_ = 0.5;
double Individual::DIAMETER_PREFERENCE_ = 0.5;
double Individual::HEIGHT_COMPETITION_ = 0.5;
double Individual::DIAMETER_COMPETITION_ = 0.5;
double Individual::TOEPAD_SELECTION_ = 0.2;
double Individual::LIMB_SELECTION_ = 0.2;
double Individual::MATING_SIGMA_ = 0.2;
double Individual::MU_LOCUS_ = 1e-4;
double Individual::MU_NEUTRAL_ = 1e-4;
double Individual::MIGRATION_RATE_ = 1e-1;

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;

//! Symbols for the program options can be different from those in equations
/*! @ingroup biol_param
    @return Program options description

    Command line option      | Symbol           | Variable
    ------------------------ | ---------------- | ------------------------------
    `-a,--beta_param`        | \f$ \alpha \f$   | Individual::BETA_PARAM_
    `-K,--carrying_capacity` | \f$ K_0 \f$      | Individual::CARRYING_CAPACITY_
    `-b,--birth_rate`        | \f$ b \f$        | Individual::AVG_NUM_OFFSPINRGS_
    `-p,--height_pref`       | \f$ h_0 \f$      | Individual::HEIGHT_PREFERENCE_
    `-P,--diameter_pref`     | \f$ h_1 \f$      | Individual::DIAMETER_PREFERENCE_
    `-c,--height_compe`      | \f$ c_0 \f$      | Individual::HEIGHT_COMPETITION_
    `-C,--diameter_compe`    | \f$ c_1 \f$      | Individual::DIAMETER_COMPETITION_
    `-s,--toepad_select`     | \f$ s_0 \f$      | Individual::TOEPAD_SELECTION_
    `-S,--limb_select`       | \f$ s_1 \f$      | Individual::LIMB_SELECTION_
    `-f,--mating_sigma`      | \f$ \sigma_a \f$ | Individual::MATING_SIGMA_
    `-u,--mu_locus`          | -                | Individual::MU_LOCUS_
    `-U,--mu_neutral`        | -                | Individual::MU_NEUTRAL_
    `-m,--migration_rate`    | \f$ m \f$        | Individual::MIGRATION_RATE_
*/
boost::program_options::options_description& Individual::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Individual"};
    desc.add_options()
        ("beta_param,a", po::value<double>(&BETA_PARAM_)->default_value(BETA_PARAM_))
        ("carrying_capacity,K", po::value<size_t>(&CARRYING_CAPACITY_)->default_value(CARRYING_CAPACITY_))
        ("birth_rate,b", po::value<size_t>(&AVG_NUM_OFFSPINRGS_)->default_value(AVG_NUM_OFFSPINRGS_))
        ("height_pref,p", po::value<double>(&HEIGHT_PREFERENCE_)->default_value(HEIGHT_PREFERENCE_))
        ("diameter_pref,P", po::value<double>(&DIAMETER_PREFERENCE_)->default_value(DIAMETER_PREFERENCE_))
        ("height_compe,c", po::value<double>(&HEIGHT_COMPETITION_)->default_value(HEIGHT_COMPETITION_))
        ("diameter_compe,C", po::value<double>(&DIAMETER_COMPETITION_)->default_value(DIAMETER_COMPETITION_))
        ("toepad_select,s", po::value<double>(&TOEPAD_SELECTION_)->default_value(TOEPAD_SELECTION_))
        ("limb_select,S", po::value<double>(&LIMB_SELECTION_)->default_value(LIMB_SELECTION_))
        ("mating_sigma,f", po::value<double>(&MATING_SIGMA_)->default_value(MATING_SIGMA_))
        ("mu_locus,u", po::value<double>(&MU_LOCUS_)->default_value(MU_LOCUS_))
        ("mu_neutral,U", po::value<double>(&MU_NEUTRAL_)->default_value(MU_NEUTRAL_))
        ("migration_rate,m", po::value<double>(&MIGRATION_RATE_)->default_value(MIGRATION_RATE_))
    ;
    return desc;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

//! Beta distribution function for resource abundance
inline double pdf_beta(const double height, const double diameter) {
    static_cast<void>(diameter);
    static const double k =
        std::tgamma(2 * Individual::BETA_PARAM_)
        / wtl::pow<2>(std::tgamma(Individual::BETA_PARAM_));
    double result = std::pow(height * (1 - height), Individual::BETA_PARAM_ - 1);
    return result *= k;
}

//! Triangular distribution function for resource abundance
inline double pdf_triangle(const double height, const double diameter) {
    if (height == 1.0 || 1.0 - height < diameter) {return 0.0;}
    double result = 1.0 - height - diameter;
    result /= wtl::pow<2>(1.0 - height);
    return result *= 2.0;
}

//! Product Beta(u, v) Tri(u, v)
inline double abundance(const double height, const double diameter) {
    return pdf_beta(height, diameter) * pdf_triangle(height, diameter);
}

inline double pdf_normal(const double height, const double diameter) {
    return std::exp(- 0.5 * wtl::pow<2>(height - 0.5) / Individual::NORMAL_SIGMA_);
}

inline double pdf_exp(const double height, const double diameter) {
    double theta = 1.0;
    theta /= Individual::C0_ - Individual::C1_ * diameter;
    return std::exp(- theta * diameter) * theta;
}

inline double abundance_old(const double height, const double diameter) {
    return pdf_normal(height, diameter) * pdf_exp(height, diameter);
}

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
}

double Individual::calc_denom_numerical() const {
    return wtl::integrate([this](const double height) {
        return wtl::integrate([this, height](const double diameter) {
            double result = habitat_preference(height, diameter);
            return result *= abundance(height, diameter);
        }, 0.0, 1.0 - height, NUM_STEPS_);
    }, 0.0, 1.0, NUM_STEPS_);
}

double Individual::calc_denom_mathematica() const {
    const double a = BETA_PARAM_;
    const double h0 = HEIGHT_PREFERENCE_;
    const double h1 = DIAMETER_PREFERENCE_;
    const double y0 = phenotype_[trait::height_preference];
    const double y1 = phenotype_[trait::diameter_preference];
    double d = 12;
    d -= 6*h0*(1 + 2*(-1 + y0)*y0);
    d += h1*(-1 + 4*(1 - 3*y1)*y1);
    d -= a*(-24 + h1 + 6*h0*wtl::pow<2>(1 - 2*y0) - 8*h1*y1 + 24*h1*wtl::pow<2>(y1));
    d /= (12 + 24*a);
    return d;
}


double Individual::calc_denom_maple() const {
    const double alpha = BETA_PARAM_;
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


double Individual::habitat_preference_v2(const double height, const double diameter) const {
    auto impl = [](double u, const double y, const double h) {
        u -= y;
        u *= u;
        return u *= h;
    };
    double exponent = -impl(phenotype_[trait::height_preference], height, HEIGHT_PREFERENCE_);
    exponent -= impl(phenotype_[trait::diameter_preference], diameter, DIAMETER_PREFERENCE_);
    return std::exp(exponent);
}

double Individual::habitat_preference(const double height, const double diameter) const {
    auto impl = [](double u, const double y, const double h) {
        u -= y;
        u *= u;
        return u *= h;
    };
    double x = 1.0;
    x -= impl(phenotype_[trait::height_preference], height, HEIGHT_PREFERENCE_);
    x -= impl(phenotype_[trait::diameter_preference], diameter, DIAMETER_PREFERENCE_);
    return std::max(x, 0.0);
}

double Individual::calc_xi_normalizer() const {
    //-h0*y0**2/2 + h0*y0/3 - h0/12 - h1*y1**2/2 + h1*y1/3 - h1/12 + 1/2
    auto impl = [](double y) {
        double result = 1.0;
        result /= 12.0;
        result -= y / 3.0;
        y *= y;
        return result += y * 0.5;
    };
    double x = 0.5;
//    auto impl = [](double y) {
//        double result = 1.0;
//        result /= 3.0;
//        result -= y;
//        y *= y;
//        return result += y;
//    };
//    double x = 1.0;
    x -= HEIGHT_PREFERENCE_ * impl(phenotype_[trait::height_preference]);
    x -= DIAMETER_PREFERENCE_ * impl(phenotype_[trait::diameter_preference]);
    return x;
}

double Individual::calc_xi_normalizer_numerical() const {
    return wtl::integrate([this](const double height) {
        return wtl::integrate([this, height](const double diameter) {
            return habitat_preference(height, diameter);
        }, 0.0, 1.0 - height, NUM_STEPS_);
    }, 0.0, 1.0, NUM_STEPS_);
}

double Individual::habitat_overlap_v2(const Individual& other) const {
    double n = wtl::integrate([this, &other](const double height) {
        return wtl::integrate([this, &other, height](const double diameter) {
            double result = habitat_preference(height, diameter);
            result *= other.habitat_preference(height, diameter);
            return result *= abundance(height, diameter);
        }, 0.0, 1.0 - height, NUM_STEPS_);
    }, 0.0, 1.0, NUM_STEPS_);
    return n;
}

double Individual::habitat_overlap(const Individual& other) const {
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
    return effective_carrying_capacity_unnormalized() / calc_xi_normalizer();
}

double Individual::effective_carrying_capacity_v3() const {
    return effective_carrying_capacity_unnormalized() / calc_denom();
}

double Individual::effective_carrying_capacity_unnormalized() const {
    double result = CARRYING_CAPACITY_;
    result *= wtl::integrate([this](const double height) {
        return wtl::integrate([this, height](const double diameter) {
            double result = fitness(height, diameter);
            result *= habitat_preference(height, diameter);
            result *= abundance(height, diameter);
            return result;
        }, 0.0, 1.0 - height, NUM_STEPS_);
    }, 0.0, 1.0, NUM_STEPS_);
    return result;
}

double Individual::survival_probability(const double effective_num_competitors) const {
    double denom = AVG_NUM_OFFSPINRGS_;
    denom -= 1.0;
    denom *= effective_num_competitors;  // ^ alpha for crowding strength
    denom /= effective_carrying_capacity();
    denom += 1.0;
    return 1.0 / denom;
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
    const std::string sep{","};
    ost << wtl::str_join(genotype_.first, sep) << sep;
    ost << wtl::str_join(genotype_.second, sep) << sep;
    ost << wtl::str_join(phenotype_, sep);
    return ost.str();
}

std::string Individual::header() {
    std::vector<std::string> names{
    "toepad", "limb",
    "height_pref", "diameter_pref",
    "male", "female", "choosiness",
    "neutral"};
    std::ostringstream ost;
    const std::string sep(",");
    for (const auto& s: names) {
        ost << s << "_L" << sep;
    }
    for (const auto& s: names) {
        ost << s << "_R" << sep;
    }
    for (const auto& s: names) {
        ost << s << "_P";
        if (s == "neutral") {ost << "\n";}
        else {ost << sep;}
    }
    return ost.str();
}

std::string Individual::str_detail() const {
    std::ostringstream ost;
    ost << Individual::header();
    ost << *this << std::endl;
    ost << gametogenesis() << std::endl;
    ost << "Ke: " << effective_carrying_capacity() << std::endl;
    ost << "DI numer: " << calc_denom_numerical() << std::endl;
    ost << "DI mathe: " << calc_denom_mathematica() << std::endl;
    ost << "DI maple: " << calc_denom_maple() << std::endl;
    ost << "Xi denom ana: " << calc_xi_normalizer() << std::endl;
    ost << "Xi denom num: " << calc_xi_normalizer_numerical() << std::endl;
    return ost.str();
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

template <class Func> inline
std::string test_resource_abundance(Func func) {
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

void Individual::write_resource_abundance() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    wtl::Fout{"ignore/abundance_beta.csv"} << ::test_resource_abundance(pdf_beta);
    wtl::Fout{"ignore/abundance_triangle.csv"} << ::test_resource_abundance(pdf_triangle);
    wtl::Fout{"ignore/abundance_v3.csv"} << ::test_resource_abundance(abundance);
    wtl::Fout{"ignore/abundance_normal.csv"} << ::test_resource_abundance(pdf_normal);
    wtl::Fout{"ignore/abundance_exp.csv"} << ::test_resource_abundance(pdf_exp);
    wtl::Fout{"ignore/abundance_old.csv"} << ::test_resource_abundance(abundance_old);
}

std::string Individual::possible_ke() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    constexpr size_t max_trait = Individual::NUM_LOCI_ * 2;
    constexpr size_t half = max_trait / 2;
    std::ostringstream ost;
    std::string sep(",");
    ost << "toepad,limb,height_pref,diameter_pref,Ke_v3,Ke_v3u,Ke_v3a,Dxi,DI\n";
    for (size_t toe=0; toe<=max_trait; ++toe) {
        for (size_t limb=0; limb<=max_trait; ++limb) {
            for (size_t hpref=0; hpref<=max_trait; ++hpref) {
                for (size_t dpref=0; dpref<=max_trait; ++dpref) {
                    Individual ind(std::vector<size_t>{toe, limb, hpref, dpref, half, half, half, half});
                    ost << toe << sep << limb << sep << hpref << sep << dpref << sep
                        << ind.effective_carrying_capacity_v3() << sep
                        << ind.effective_carrying_capacity_unnormalized() << sep
                        << ind.effective_carrying_capacity() << sep
                        << ind.calc_xi_normalizer() << sep
                        << ind.calc_denom() << "\n";
                }
            }
        }
    }
    return ost.str();
}

std::string Individual::sojourn_time(const bool normalizing) const {
    constexpr size_t max_trait = Individual::NUM_LOCI_ * 2;
    constexpr double inv_max = 1.0 / max_trait;
    const std::string sep = ",";
    std::vector<std::string> lines;
    lines.reserve(max_trait * max_trait);
    std::vector<double> preference{phenotype_.begin()+2, phenotype_.begin()+4};
    for (size_t ih=0; ih<max_trait; ++ih) {
        const double height = ih * inv_max;
        for (size_t id=0; id<max_trait; ++id) {
            const double diameter = id * inv_max;
            double result = 1.0;
            if (normalizing) {
                result /= calc_denom();
            }
            result *= habitat_preference(height, diameter);
            result *= abundance(height, diameter);
            std::ostringstream ost;
            ost.precision(16);
            ost << preference[0] << sep << preference[1] << sep
                << height << sep << diameter << sep
                << std::boolalpha << normalizing << sep
                << std::scientific << result;
            lines.push_back(ost.str());
        }
    }
    return wtl::str_join(lines, "\n");
}

std::string Individual::test_sojourn_time() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    constexpr size_t max_trait = Individual::NUM_LOCI_ * 2;
    constexpr size_t half = max_trait / 2;
    std::ostringstream ost;
    std::string sep(",");
    ost << "height_pref,diameter_pref,height,diameter,normalized,time\n";
    for (size_t hpref=0; hpref<=max_trait; ++hpref) {
        for (size_t dpref=0; dpref<=max_trait; ++dpref) {
            Individual ind(std::vector<size_t>{half, half, hpref, dpref, half, half, half, half});
            ost << ind.sojourn_time(true) << "\n";
            ost << ind.sojourn_time(false) << "\n";
        }
    }
    return ost.str();
}

void Individual::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Individual ind;
    std::cerr << ind.str_detail();
    Individual offspring(ind.gametogenesis(), ind.gametogenesis());
    std::cerr << offspring << std::endl;
}
