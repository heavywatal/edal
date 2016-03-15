// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/gz.hpp>

double Individual::BETA_PARAM_ = 3.0;
double Individual::NORMAL_SIGMA_ = 0.3;
double Individual::C0_ = 1.0;
double Individual::C1_ = 0.5;
size_t Individual::CARRYING_CAPACITY_ = 1000;
double Individual::AVG_NUM_OFFSPINRGS_ = 4;
double Individual::HEIGHT_PREFERENCE_ = 2.0;
double Individual::DIAMETER_PREFERENCE_ = 2.0;
double Individual::TOEPAD_SELECTION_ = 2.0;
double Individual::LIMB_SELECTION_ = 2.0;
double Individual::PREF_COMPETITION_ = 2.0;
double Individual::MORPH_COMPETITION_ = 2.0;
double Individual::MATING_SIGMA_ = 0.05;
double Individual::MU_LOCUS_ = 1e-4;
unsigned long Individual::MUTATION_MASK_ = 0;
double Individual::MIGRATION_RATE_ = 0.005;

constexpr size_t Individual::NUM_LOCI_;
constexpr unsigned long Individual::FULL_BITS;
constexpr unsigned long Individual::HALF_BITS;
constexpr double Individual::INV_NUM_LOCI_;
std::map<std::vector<double>, double> Individual::KE_CACHE_;

//! Symbols for the program options can be different from those in equations
/*! @ingroup biol_param
    @return Program options description

    Command line option      | Symbol         | Variable
    ------------------------ | -------------- | --------------------------------
    `-a,--beta_param`        | \f$\alpha\f$   | Individual::BETA_PARAM_
    `-K,--carrying_capacity` | \f$K_{\max}\f$ | Individual::CARRYING_CAPACITY_
    `-b,--birth_rate`        | \f$b\f$        | Individual::AVG_NUM_OFFSPINRGS_
    `-p,--height_pref`       | \f$h_0\f$      | Individual::HEIGHT_PREFERENCE_
    `-P,--diameter_pref`     | \f$h_1\f$      | Individual::DIAMETER_PREFERENCE_
    `-s,--toepad_select`     | \f$s_0\f$      | Individual::TOEPAD_SELECTION_
    `-S,--limb_select`       | \f$s_1\f$      | Individual::LIMB_SELECTION_
    `-c,--pref_compe`        | \f$c_y\f$      | Individual::PREF_COMPETITION_
    `-C,--morph_compe`       | \f$c_x\f$      | Individual::MORPH_COMPETITION_
    `-f,--mating_sigma`      | \f$\sigma_a\f$ | Individual::MATING_SIGMA_
    `-u,--mu_locus`          | -              | Individual::MU_LOCUS_
    `-U,--mutation_mask`     | -              | Individual::MUTATION_MASK_
    `-m,--migration_rate`    | \f$m\f$        | Individual::MIGRATION_RATE_
*/
boost::program_options::options_description& Individual::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Individual"};
    desc.add_options()
        ("beta_param,a", po::value<double>(&BETA_PARAM_)->default_value(BETA_PARAM_))
        ("carrying_capacity,K", po::value<size_t>(&CARRYING_CAPACITY_)->default_value(CARRYING_CAPACITY_))
        ("birth_rate,b", po::value<double>(&AVG_NUM_OFFSPINRGS_)->default_value(AVG_NUM_OFFSPINRGS_))
        ("height_pref,p", po::value<double>(&HEIGHT_PREFERENCE_)->default_value(HEIGHT_PREFERENCE_))
        ("diameter_pref,P", po::value<double>(&DIAMETER_PREFERENCE_)->default_value(DIAMETER_PREFERENCE_))
        ("toepad_select,s", po::value<double>(&TOEPAD_SELECTION_)->default_value(TOEPAD_SELECTION_))
        ("limb_select,S", po::value<double>(&LIMB_SELECTION_)->default_value(LIMB_SELECTION_))
        ("pref_compe,c", po::value<double>(&PREF_COMPETITION_)->default_value(PREF_COMPETITION_))
        ("morph_compe,C", po::value<double>(&MORPH_COMPETITION_)->default_value(MORPH_COMPETITION_))
        ("mating_sigma,f", po::value<double>(&MATING_SIGMA_)->default_value(MATING_SIGMA_))
        ("mu_locus,u", po::value<double>(&MU_LOCUS_)->default_value(MU_LOCUS_))
        ("mutation_mask,U", po::value<unsigned long>(&MUTATION_MASK_)->default_value(MUTATION_MASK_))
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
        / wtl::pow(std::tgamma(Individual::BETA_PARAM_), 2);
    double result = std::pow(height * (1 - height), Individual::BETA_PARAM_ - 1);
    return result *= k;
}

//! Triangular distribution function for resource abundance
inline double pdf_triangle(const double height, const double diameter) {
    if (height == 1.0 || 1.0 - height < diameter) {return 0.0;}
    double result = 1.0;
    result -= height;
    result -= diameter;
    result /= wtl::pow(1.0 - height, 2);
    return result *= 2.0;
}

//! Product Beta(u, v) Tri(u, v)
inline double abundance(const double height, const double diameter) {
    return pdf_beta(height, diameter) * pdf_triangle(height, diameter);
}

//! Distribution of tree height given diameter
inline double pdf_normal(const double height, const double diameter) {
    return std::exp(- 0.5 * wtl::pow(height - 0.5, 2) / wtl::pow(Individual::NORMAL_SIGMA_, 2));
}

//! Distribution of twig diameter given height
inline double pdf_exp(const double height, const double diameter) {
    const double theta = Individual::C0_ - Individual::C1_ * height;
    if (theta <= 0) {return 0.0;}
    double lambda = 1.0;
    lambda /= theta;
    return lambda *= std::exp(- lambda * diameter);
}

//! Product of Normal(u, v) and Exponential(u, v)
inline double abundance_old(const double height, const double diameter) {
    return pdf_normal(height, diameter) * pdf_exp(height, diameter);
}

//! Integral over \f$0 \leq u \leq 1\f$ and \f$0 \leq v \leq 1 - u\f$
template <class Func> inline
double integrate_triangle(Func&& func) {
    return wtl::integrate([&func](const double height) {
        return wtl::integrate([&func, height](const double diameter) {
            return func(height, diameter);
        }, 0.0, 1.0 - height, Individual::NUM_STEPS_);
    }, 0.0, 1.0, Individual::NUM_STEPS_);
}

//! Integral over \f$0 \leq u,v \leq 1\f$
template <class Func> inline
double integrate_square(Func&& func) {
    return wtl::integrate([&func](const double height) {
        return wtl::integrate([&func, height](const double diameter) {
            return func(height, diameter);
        }, 0.0, 1.0, Individual::NUM_STEPS_);
    }, 0.0, 1.0, Individual::NUM_STEPS_);
}

//! Make exponent of Gaussian function
inline double gaussian_exponent(double x, const double mu, const double sigma) {
    x -= mu;
    x /= sigma;
    x *= x;
    return x *= -0.5;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Individual::Individual(const std::vector<unsigned long>& flags): genotype_{{}, {}} {
    genotype_.first.reserve(trait::size);
    for (auto x: flags) {
        genotype_.first.emplace_back(x);
    }
    for (size_t i=genotype_.first.size(); i<trait::size; ++i) {
        genotype_.first.push_back(HALF_BITS);
    }
    assert(genotype_.first.size() == trait::size);
    genotype_.second = genotype_.first;
    assert(genotype_.second.size() == trait::size);
    phenotype_ = calc_phenotype();
    ke_ = effective_carrying_capacity_cache();
}

double Individual::habitat_preference_exp(const double height, const double diameter) const {
    double exponent = gaussian_exponent(phenotype_[trait::height_preference],
                                        height, HEIGHT_PREFERENCE_);
    return exponent += gaussian_exponent(phenotype_[trait::diameter_preference],
                                         diameter, DIAMETER_PREFERENCE_);
}

double Individual::habitat_preference_quadratic(const double height, const double diameter) const {
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

double Individual::calc_DI_analytical() const {
    const double a = BETA_PARAM_;
    const double h0 = HEIGHT_PREFERENCE_;
    const double h1 = DIAMETER_PREFERENCE_;
    const double y0 = phenotype_[trait::height_preference];
    const double y1 = phenotype_[trait::diameter_preference];
    double d = 12;
    d -= 6*h0*(1 + 2*(-1 + y0)*y0);
    d += h1*(-1 + 4*(1 - 3*y1)*y1);
    d -= a*(-24 + h1 + 6*h0*wtl::pow(1 - 2*y0, 2) - 8*h1*y1 + 24*h1*wtl::pow(y1, 2));
    d /= (12 + 24*a);
    return d;
}

double Individual::calc_DI_numerical() const {
    return integrate_triangle([this](const double u, const double v) {
        double result = habitat_preference_quadratic(u, v);
        return result *= abundance(u, v);
    });
}

double Individual::calc_Dxi_analytical() const {
    auto impl = [](double y) {
        double result = 1.0;
        result /= 12.0;
        result -= y / 3.0;
        y *= y;
        y *= 0.5;
        return result += y;
    };
    double x = 0.5;
    x -= HEIGHT_PREFERENCE_ * impl(phenotype_[trait::height_preference]);
    x -= DIAMETER_PREFERENCE_ * impl(phenotype_[trait::diameter_preference]);
    return x;
}

double Individual::calc_Dxi_numerical() const {
    return integrate_triangle([this](const double u, const double v) {
        return habitat_preference_quadratic(u, v);
    });
}

double Individual::fitness(const double height, const double diameter) const {
    double exponent = gaussian_exponent(phenotype_[trait::toepad_size],
                                        height, TOEPAD_SELECTION_);
    return exponent += gaussian_exponent(phenotype_[trait::limb_length],
                                  diameter, LIMB_SELECTION_);
}

double Individual::effective_carrying_capacity_quad_unnormalized() const {
    double result = CARRYING_CAPACITY_;
    return result *= integrate_triangle([this](const double u, const double v) {
        double result = std::exp(fitness(u, v));
        result *= habitat_preference_quadratic(u, v);
        return result *= abundance(u, v);
    });
}

double Individual::effective_carrying_capacity_exp_unnormalized() const {
    double result = CARRYING_CAPACITY_;
    return result *= integrate_triangle([this](const double u, const double v) {
        double result = fitness(u, v);
        result += habitat_preference_exp(u, v);
        return std::exp(result) * abundance(u, v);
    });
}

double Individual::effective_carrying_capacity_old_exp_unnormalized() const {
    double result = CARRYING_CAPACITY_;
    return result *= integrate_square([this](const double u, const double v) {
        double result = fitness(u, v);
        result += habitat_preference_exp(u, v);
        return std::exp(result) * abundance_old(u, v);
    });
}

double Individual::effective_carrying_capacity_cache() const {
    std::vector<double> ecol_traits(phenotype_.begin(), phenotype_.begin()+4);
    if (KE_CACHE_.find(ecol_traits) == KE_CACHE_.end()) {
        const double ke = effective_carrying_capacity_exp_unnormalized();
        return KE_CACHE_[ecol_traits] = ke;
    } else {
        return KE_CACHE_[ecol_traits];
    }
}

double Individual::preference_overlap(const Individual& other) const {
    double exponent = gaussian_exponent(phenotype_[trait::height_preference],
                                  other.phenotype_[trait::height_preference],
                                  PREF_COMPETITION_);
    return exponent += gaussian_exponent(phenotype_[trait::diameter_preference],
                                   other.phenotype_[trait::diameter_preference],
                                   PREF_COMPETITION_);
}

double Individual::morphology_overlap(const Individual& other) const {
    double exponent = gaussian_exponent(phenotype_[trait::toepad_size],
                                  other.phenotype_[trait::toepad_size],
                                  MORPH_COMPETITION_);
    return exponent += gaussian_exponent(phenotype_[trait::limb_length],
                                   other.phenotype_[trait::limb_length],
                                   MORPH_COMPETITION_);
}

double Individual::survival_probability(const double effective_num_competitors) const {
    double denom = AVG_NUM_OFFSPINRGS_;
    denom -= 1.0;
    denom *= effective_num_competitors;  // ^ theta for crowding strength
    denom /= ke_;
    denom += 1.0;
    return 1.0 / denom;
}

double Individual::mating_probability(const Individual& male) const {
    double choosiness = phenotype_[trait::choosiness];
    if (choosiness == 0.5) {
        // no assortative mating
        return std::exp(preference_overlap(male));
    }
    choosiness *= 2.0;
    choosiness -= 1.0;
    double exponent = phenotype_[trait::female_trait];
    if (choosiness > 0.0) {
        // assortative mating
        exponent -= male.phenotype_[trait::male_trait];
    } else {
        // disassortative mating
        exponent -= 1.0;
        exponent += male.phenotype_[trait::male_trait];
    }
    exponent /= MATING_SIGMA_;
    exponent *= choosiness;
    exponent *= exponent;
    exponent *= -0.5;
    exponent += preference_overlap(male);
    return std::exp(exponent);
}

double Individual::mating_probability_debarre(const Individual& male) const {
    double choosiness = phenotype_[trait::choosiness];
    if (choosiness == 0.5) {
        // no assortative mating
        return std::exp(preference_overlap(male));
    }
    choosiness *= 2.0;
    choosiness -= 1.0;
    double exponent = phenotype_[trait::female_trait];
    exponent -= male.phenotype_[trait::male_trait];
    exponent /= MATING_SIGMA_;
    exponent *= exponent;
    exponent *= -0.5;
    double phi = choosiness;
    phi *= choosiness;
    if (choosiness > 0.0) {
        // assortative mating
        phi *= (1.0 - std::exp(exponent));
    } else {
        // disassortative mating
        phi *= std::exp(exponent);
    }
    phi -= 1.0;
    return -phi * std::exp(preference_overlap(male));
}

double Individual::mating_probability_TPG2013(const Individual& male) const {
    double choosiness = phenotype_[trait::choosiness];
    if (choosiness == 0.5) {
        // no assortative mating
        return std::exp(preference_overlap(male));
    }
    choosiness *= 2.0;
    choosiness -= 1.0;
    double exponent = phenotype_[trait::female_trait];
    exponent -= male.phenotype_[trait::male_trait];
    exponent /= MATING_SIGMA_;
    exponent *= choosiness;
    exponent *= exponent;
    exponent *= -0.5;
    if (choosiness > 0.0) {
        // assortative mating
        return std::exp(exponent += preference_overlap(male));
    } else {
        // disassortative mating
        double psi = 2.0;
        psi -= std::exp(exponent);
        return psi *= std::exp(preference_overlap(male));
    }
}

std::vector<Individual::Loci> Individual::gametogenesis(wtl::sfmt19937& rng) const {
    std::uniform_int_distribution<unsigned long> uniform_filter(0, FULL_BITS);
    std::bernoulli_distribution bernoulli(MU_LOCUS_ * NUM_LOCI_);
    std::uniform_int_distribution<size_t> uniform_pos(0, NUM_LOCI_ - 1);
    static const std::bitset<8> mask(MUTATION_MASK_);
    std::vector<Loci> gamete;
    gamete.reserve(trait::size);
    for (size_t i=0; i<trait::size; ++i) {
        if (mask.test(i)) {
            gamete.push_back(genotype_.first[i]);
            continue;
        }
        // recombination
        const Loci filter(uniform_filter(rng));
        gamete.push_back((genotype_.first[i] & filter) | (genotype_.second[i] & ~filter));
        // mutation
        if (bernoulli(rng)) {
            gamete.back().flip(uniform_pos(rng));
        }
    }
    return gamete;
}

std::ostream& operator<< (std::ostream& ost, const Individual& ind) {
    const char* sep = ",";
    return ost
        << wtl::str_join(ind.genotype_.first, sep) << sep
        << wtl::str_join(ind.genotype_.second, sep) << sep
        << wtl::str_join(ind.phenotype_, sep);
}

std::string Individual::header() {
    std::vector<std::string> names{
    "toepad", "limb",
    "height_pref", "diameter_pref",
    "male", "female", "choosiness",
    "neutral"};
    std::ostringstream ost;
    const char* sep = ",";
    for (const auto& s: names) {
        ost << s << "_L" << sep;
    }
    for (const auto& s: names) {
        ost << s << "_R" << sep;
    }
    for (const auto& s: names) {
        ost << s << "_P";
        if (s != "neutral") {ost << sep;}
    }
    return ost.str();
}

std::vector<double> Individual::intermediate_phenotypes() const {
    return {
        effective_carrying_capacity_quad_unnormalized(),
        effective_carrying_capacity_exp_unnormalized(),
        effective_carrying_capacity_old_exp_unnormalized(),
        calc_DI_analytical(),
        calc_DI_numerical(),
        calc_Dxi_analytical(),
        calc_Dxi_numerical(),
        integrate_triangle([this](const double u, const double v) {
            return habitat_preference_exp(u, v);
        }),
        integrate_square([this](const double u, const double v) {
            return habitat_preference_exp(u, v);
        })};
}

const std::vector<std::string> Individual::INTERMEDIATE_KEYS_ = {
    "Ke_v3u", "Ke_v3eu", "Ke_oeu", "DI_a", "DI",
    "Dxi_a", "Dxi", "Dxi_e", "Dxi_oe"
};

std::string Individual::str_detail() const {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::ostringstream ost;
    ost << Individual::header() << std::endl;
    ost << *this << std::endl;
    const auto values = Individual::intermediate_phenotypes();
    for (size_t i=0; i<values.size(); ++i) {
        ost << Individual::INTERMEDIATE_KEYS_[i] << ": "
            << values[i] << std::endl;
    }
    return ost.str();
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

//! CSV string of \f$F(u,v)\f$ in \f$0 \leq u,v \leq 1\f$
template <class Func> inline
std::string test_resource_abundance(Func func) {
    constexpr int precision = 25;
    constexpr double delta = 1.0 / precision;
    std::ostringstream ost;
    const char* sep = ",";
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
    wtl::Fout{"abundance_beta.csv"} << ::test_resource_abundance(pdf_beta);
    wtl::Fout{"abundance_triangle.csv"} << ::test_resource_abundance(pdf_triangle);
    wtl::Fout{"abundance_v3.csv"} << ::test_resource_abundance(abundance);
    wtl::Fout{"abundance_normal.csv"} << ::test_resource_abundance(pdf_normal);
    wtl::Fout{"abundance_exp.csv"} << ::test_resource_abundance(pdf_exp);
    wtl::Fout{"abundance_old.csv"} << ::test_resource_abundance(abundance_old);
}

std::string Individual::possible_phenotypes() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    constexpr size_t max_trait = NUM_LOCI_ * 2;
    constexpr size_t half = max_trait / 2;
    std::ostringstream ost;
    const char* sep = ",";
    ost << "toepad,limb,height_pref,diameter_pref,"
        << wtl::str_join(INTERMEDIATE_KEYS_, sep)
        << "\n";
    for (size_t toe=0; toe<=max_trait; ++toe, ++toe) {
        for (size_t limb=0; limb<=max_trait; ++limb, ++limb) {
            for (size_t hpref=0; hpref<=max_trait; ++hpref, ++hpref) {
                for (size_t dpref=0; dpref<=max_trait; ++dpref, ++dpref) {
                    Individual ind(std::vector<size_t>{toe, limb, hpref, dpref, half, half, half, half});
                    ost << toe << sep << limb << sep << hpref << sep << dpref << sep
                        << wtl::str_join(ind.intermediate_phenotypes(), sep) << "\n";
                }
            }
        }
    }
    return ost.str();
}

std::string Individual::possible_geographic() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    constexpr size_t max_trait = NUM_LOCI_ * 2;
    constexpr size_t half = max_trait / 2;
    constexpr double max_reciprocal = 1.0 / max_trait;
    std::ostringstream ost;
    const char* sep = ",";
    ost << "toepad,limb,height_pref,diameter_pref,"
        << "DI,Dxi,height,diameter,"
        << "resource,Xi_quad,Xi_exp,fitness\n";
    for (size_t toe=0; toe<=max_trait; ++toe, ++toe) {
        for (size_t limb=0; limb<=max_trait; ++limb, ++limb) {
            for (size_t hpref=0; hpref<=max_trait; ++hpref, ++hpref) {
                for (size_t dpref=0; dpref<=max_trait; ++dpref, ++dpref) {
                    Individual ind(std::vector<size_t>{toe, limb, hpref, dpref, half, half, half, half});
                    const double t_denom = ind.calc_DI_numerical();
                    const double xi_denom = ind.calc_Dxi_numerical();
                    for (size_t hi=0; hi<=max_trait; ++hi, ++hi) {
                        const double height = hi * max_reciprocal;
                        for (size_t di=0; di<=max_trait; ++di, ++di) {
                            const double diameter = di * max_reciprocal;
                            ost << toe << sep << limb << sep << hpref << sep << dpref << sep
                                << t_denom << sep << xi_denom << sep
                                << hi << sep << di << sep
                                << abundance(height, diameter) << sep
                                << ind.habitat_preference_quadratic(height, diameter) << sep
                                << ind.habitat_preference_exp(height, diameter) << sep
                                << ind.fitness(height, diameter)
                                << "\n";
                        }
                    }
                }
            }
        }
    }
    return ost.str();
}

void Individual::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Individual ind;
    std::cerr << ind.str_detail();
}
