// -*- mode: c++; coding: utf-8 -*-
#pragma once
#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include <iostream>
#include <vector>
#include <bitset>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace wtl {
namespace detail {

template <unsigned int N, class T> struct Pow {
    constexpr T operator()(T x) const {
        return x * Pow<N - 1, T>()(x);
    }
};

template <class T> struct Pow<1, T> {
    constexpr T operator()(const T& x) const {return x;}
};

template <class T> struct Pow<0, T> {
    constexpr T operator()(const T& x) const {return x * 0 + 1;}
};

} // namespace detail

// interger powers at compile time
template <unsigned int N, class T> inline constexpr
T pow(const T& x) {
    return detail::Pow<N, T>()(x);
}
} // namespace wtl

namespace trait {
enum {
    // morphological characters
    toepad_size,
    limb_length,

    // habitat preference characters
    height_preference,
    diameter_preference,

    // mating compatibility characters
    male_trait,
    female_trait,
    choosiness,

    // to calculate genetic divergence
    neutral,

    // for computational convenience
    size
};
} // namespace trait

class Individual {
    // parameters
  private:
    static size_t CARRYING_CAPACITY;
    static size_t AVG_NUM_OFFSPINRGS_;
    static double HABITAT_SIGMA_;
    static double MATING_SIGMA_;
    static double ADAPTIVE_T_SIGMA_;
    static double ADAPTIVE_L_SIGMA_;
    static double MU_LOCUS_;
    static double MU_NEUTRAL_;

    constexpr static size_t NUM_LOCI_ = 8;  // per trait
    constexpr static unsigned long FULL_BITS = wtl::pow<NUM_LOCI_>(2) - 1;
    constexpr static unsigned long HALF_BITS = wtl::pow<NUM_LOCI_ / 2>(2) - 1;
    constexpr static double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;

    typedef std::bitset<NUM_LOCI_> Loci;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:

    Individual(): genotype_{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS}},
        phenotype_(init_phenotype()),
        denominator_{denom_()}, sqrt_denominator_2_{sqrt_denom_2_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm}, phenotype_(init_phenotype()),
        denominator_{denom_()}, sqrt_denominator_2_{sqrt_denom_2_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<size_t>&);

    double effective_carrying_capacity() const;
    double effective_num_competitors(const Individual&) const;

    bool survive(const double effective_population_size) const;

    double mating_preference(const Individual& male) const;
    size_t poisson_offsprings() const;
    std::vector<Loci> gametogenesis() const;

    std::string str() const;

    static boost::program_options::options_description& opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    std::vector<double> init_phenotype() const {
        const size_t n(trait::size);
        std::vector<double> output(n);
        for (size_t i=0; i<n; ++i) {
            output[i] = (genotype_.first[i].count() + genotype_.second[i].count()) * INV_NUM_LOCI_;
        }
        return output;
    };
    double habitat_preference(const double height, const double diameter) const;
    double denom_numerical() const;
    double denom_mathematica() const;
    double denom_maple() const;
    double denom_() const {return denom_mathematica();}
    double sqrt_denom_2_() const;
    double fitness(const double height, const double diameter) const;

    static Loci recombination(const Loci&, const Loci&);
    static std::vector<Loci> mutate(std::vector<Loci>);

    friend void individual_unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member (genotype)
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;
    std::vector<double> phenotype_;
    double denominator_;
    double sqrt_denominator_2_;
    double effective_carrying_capacity_;
};

inline std::ostream& operator<< (std::ostream& ost, const Individual& ind) {
    return ost << ind.str();
}
extern void individual_unit_test();

#endif /* INDIVIDUAL_H_ */
