// -*- mode: c++; coding: utf-8 -*-
#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <vector>
#include <bitset>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

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
    static constexpr unsigned long FULL_BITS = pow<NUM_LOCI_>(2) - 1;
    static constexpr unsigned long HALF_BITS = pow<NUM_LOCI_ / 2>(2) - 1;
    static constexpr double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;
    static constexpr double STRENGTH_OF_MATING_PREFERENCE_ = 0.05;
    static constexpr size_t NUM_FITTEST_OFSPINRGS_ = 4;
    static constexpr size_t MAX_CARRYING_CAPACITY_ = 30;

    typedef std::bitset<NUM_LOCI_> Loci;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:

    Individual(): genotype_{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS}} {}

    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm} {}

    double mate(const Individual& male);

    // get phenotype
    double phenotype(size_t t) const {
        return (scale(genotype_.first[t]) + scale(genotype_.second[t])) * 0.5;
    }

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    static constexpr double scale(const Loci& loci) {
        // completely additive
        return loci.count() * INV_NUM_LOCI_;
    }
    static constexpr double scale(const std::pair<Loci, Loci>& p) {
        // completely additive
        return (p.first.count() + p.second.count()) * INV_NUM_LOCI_;
    }
    static Loci recombination(const Loci&, const Loci&);
    std::vector<Loci> gametogenesis() const;

    friend void individual_unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member (genotype)
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;

};

extern void individual_unit_test();


#endif /* INDIVIDUAL_H_ */
