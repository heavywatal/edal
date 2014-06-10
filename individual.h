// -*- mode: c++; coding: utf-8 -*-
#pragma once
#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include <iostream>
#include <vector>
#include <bitset>
#include <string>

#include "cxxwtils/algorithm.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace trait {
enum {
    //! morphological characters
    toepad_size,
    limb_length,

    //! habitat preference characters
    height_preference,
    diameter_preference,

    //! mating compatibility characters
    male_trait,
    female_trait,
    choosiness,

    //! to calculate genetic divergence
    neutral,

    //! for computational convenience
    size
};
} // namespace trait

class Individual {
    // parameters
  private:
    static size_t CARRYING_CAPACITY;
    static size_t AVG_NUM_OFFSPINRGS_;
    static double HEIGHT_PREFERENCE_;
    static double DIAMETER_PREFERENCE_;
    static double HEIGHT_COMPETITION_;
    static double DIAMETER_COMPETITION_;
    static double TOEPAD_SELECTION_;
    static double LIMB_SELECTION_;
    static double MATING_SIGMA_;
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
        denominator_{denom_()}, //sqrt_denominator_2_{sqrt_denom_2_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm},
        phenotype_(init_phenotype()),
        denominator_{denom_()}, //sqrt_denominator_2_{sqrt_denom_2_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<size_t>&);

    //! K_e(I)
    double effective_carrying_capacity() const;

    //! C(I, J)
    /** @param other individual to interact
        @return C(I, J)
    */
    double habitat_overlap(const Individual& other) const {
        return habitat_overlap_v3(other);
    }

    //! C(I, J) in anolis_v2
    double habitat_overlap_v2(const Individual&) const;

    //! C(I, J) in anolis_v3
    double habitat_overlap_v3(const Individual&) const;

    //! w(I)
    bool survive(const double effective_num_competitors) const;

    //! ψ(I, I') [Psi]
    double mating_preference(const Individual& male) const;

    //! generate poisson random number with lambda = AVG_NUM_OFFSPINRGS_
    /** return the number of offsprings
    */
    size_t poisson_offsprings() const;

    //! generate with recombination
    /** return a gamete
    */
    std::vector<Loci> gametogenesis() const;

    //! CSV formated string
    /** @return CSV formated string
    */
    std::string str() const;

    //! the header correspongs to str()
    static std::string header();

    static boost::program_options::options_description& opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! calculate phenotypic values from genotype
    /** @return phenotypic values
    */
    std::vector<double> init_phenotype() const {
        constexpr size_t n = trait::size;
        std::vector<double> output(n);
        for (size_t i=0; i<n; ++i) {
            output[i] = (genotype_.first[i].count() + genotype_.second[i].count()) * INV_NUM_LOCI_;
        }
        return output;
    };

    //! Ξ(I, u, v) [Xi] anolis_v2
    double habitat_preference_v2(const double height, const double diameter) const;

    //! Ξ(I, u, v) [Xi] anolis_v3
    double habitat_preference_v3(const double height, const double diameter) const;

    //! Ξ(I, u, v) [Xi]
    /** @param height habitat environment
        @param diameter habitat environment
        @return Ξ(I, u, v)
    */
    double habitat_preference(const double height, const double diameter) const {
        return habitat_preference_v3(height, diameter);
    }

    //! D_I numerical computation (slow)
    double denom_numerical() const;

    //! D_I analytical computation by Mathematica (fast)
    double denom_mathematica() const;

    //! D_I analytical computation by Maple
    //! @bug something wrong
    double denom_maple() const;

    //! D_I
    /** @return D_I
    */
    double denom_() const {return denom_mathematica();}

    double sqrt_denom_2_() const;

    //! W(I, u, v)
    /** @param height habitat environmant
        @param diameter habitat environment
    */
    double fitness(const double height, const double diameter) const;

    //! free recombination
    /** @param lhs left arm of a chromosome
        @param rhs right arm of a chromosome
        @return haplotype after recombination
    */
    static Loci recombination(const Loci&, const Loci&);
    
    //! poisson process with lambda = MU_LOCUS_ * NUM_LOCI_ * trait::size
    /** @param haplotype (call-by-value)
        @return mutated haplotype
    */
    static std::vector<Loci> mutate(std::vector<Loci>);

    friend void individual_unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member (genotype)
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;
    std::vector<double> phenotype_;
    double denominator_;
//    double sqrt_denominator_2_;
    double effective_carrying_capacity_;
};

inline std::ostream& operator<< (std::ostream& ost, const Individual& ind) {
    return ost << ind.str();
}
extern void individual_unit_test();

#endif /* INDIVIDUAL_H_ */
