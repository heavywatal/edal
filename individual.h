// -*- mode: c++; coding: utf-8 -*-
#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include <iostream>
#include <vector>
#include <bitset>
#include <string>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

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
    height_tolerance,
    diameter_tolerance,

    // mating compatibility characters
    male_trait,
    female_trait,
    choosiness,
    size
};
} // namespace trait

class Individual {
  private:
    // parameters
    static constexpr size_t NUM_LOCI_ = 8;  // per trait
    static constexpr unsigned long FULL_BITS = wtl::pow<NUM_LOCI_>(2) - 1;
    static constexpr unsigned long HALF_BITS = wtl::pow<NUM_LOCI_ / 2>(2) - 1;
    static constexpr double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;
    static constexpr double STRENGTH_OF_MATING_PREFERENCE_ = 0.05;
    static constexpr size_t NUM_FITTEST_OFSPINRGS_ = 4;
    static constexpr size_t MAX_CARRYING_CAPACITY_ = 30;
    static constexpr double MU_LOCUS_ = 1e-5;

    typedef std::bitset<NUM_LOCI_> Loci;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:

    Individual(): genotype_{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, FULL_BITS, FULL_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, FULL_BITS, FULL_BITS, HALF_BITS, HALF_BITS, HALF_BITS}} {}

    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm} {}

    Individual(const std::vector<size_t>&);

    double preference(const double height, const double diameter) const;
    double denominator() const;
    double denominator_2() const;
    double denominator_prime() const;
    double denominator_prime_2() const;
    double fitness(const double height, const double diameter) const;
    double carrying_capacity() const;
    double carrying_capacity_prime() const;

    bool survive(const double, const double) const;

    double mate(const Individual& male);

    // get phenotype
    double phenotype(size_t t) const {
        return (genotype_.first[t].count() + genotype_.second[t].count()) * INV_NUM_LOCI_;
    }
    std::string str() const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    static Loci recombination(const Loci&, const Loci&);
    static std::vector<Loci> mutate(std::vector<Loci>);
    std::vector<Loci> gametogenesis() const;

    template <class Func> inline
    double integral(Func func) const {
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

    friend void individual_unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member (genotype)
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;

};

extern void individual_unit_test();


#endif /* INDIVIDUAL_H_ */
