// -*- mode: c++; coding: utf-8 -*-
/** @file individual.h
    @brief Interface of Individual class
*/
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
  private:
    //! \f$k\f$ in \f$F(u,v)\f$
    static size_t CARRYING_CAPACITY;

    //! \f$b\f$ in \f$w(I)\f$
    static size_t AVG_NUM_OFFSPINRGS_;

    //! \f$h_0\f$ in \f$\Xi(y_0,y_1|u,v)\f$
    static double HEIGHT_PREFERENCE_;

    //! \f$h_1\f$ in \f$\Xi(y_0,y_1|u,v)\f$
    static double DIAMETER_PREFERENCE_;

    //! \f$c_0\f$ in \f$C(I,J)\f$
    static double HEIGHT_COMPETITION_;

    //! \f$c_1\f$ in \f$C(I,J)\f$
    static double DIAMETER_COMPETITION_;

    //! \f$s_0\f$ in \f$W(x_0,x_1|u,v)\f$
    static double TOEPAD_SELECTION_;

    //! \f$s_1\f$ in \f$W(x_0,x_1|u,v)\f$
    static double LIMB_SELECTION_;

    //! \f$\sigma_a\f$ in \f$\Psi()\f$
    static double MATING_SIGMA_;

    //! Mutation rate per locus per generation
    /** Mutations occur at equal rates across all loci;
        the probabilities of forward and backward mutations are equal.
    */
    static double MU_LOCUS_;

    //! Unused yet
    static double MU_NEUTRAL_;

    //! The number of loci per trait
    constexpr static size_t NUM_LOCI_ = 8;
    constexpr static unsigned long FULL_BITS = wtl::pow<NUM_LOCI_>(2) - 1;
    constexpr static unsigned long HALF_BITS = wtl::pow<NUM_LOCI_ / 2>(2) - 1;
    constexpr static double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;

    //! typedef for diallelic loci of a trait
    typedef std::bitset<NUM_LOCI_> Loci;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:

    Individual(): genotype_{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS}},
        phenotype_(init_phenotype()),
        denominator_{denom_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm},
        phenotype_(init_phenotype()),
        denominator_{denom_()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    Individual(const std::vector<size_t>&);

    //! \f$K_e(I)\f$
    double effective_carrying_capacity() const;

    //! \f$C(I, J)\f$
    /** @param other individual to interact
        @return \f$C(I, J)\f$
    */
    double habitat_overlap(const Individual& other) const {
        return habitat_overlap_v3(other);
    }

    //! \f$C(I, J)\f$ in anolis_v2
    double habitat_overlap_v2(const Individual&) const;

    //! \f$C(I, J)\f$ in anolis_v3
    double habitat_overlap_v3(const Individual&) const;

    //! \f$w(I)\f$
    bool survive(const double effective_num_competitors) const;

    //! \f$\Psi(I, I')\f$
    double mating_preference(const Individual& male) const;

    //! generates poisson random number with lambda = AVG_NUM_OFFSPINRGS_
    /** return the number of offsprings
    */
    size_t poisson_offsprings() const;

    //! Gametogenesis with free recombination among loci
    /** return a gamete
    */
    std::vector<Loci> gametogenesis() const;

    //! CSV formated string
    /** @return CSV formated string
    */
    std::string str() const;

    //! The header correspongs to str()
    static std::string header();

    //! Program options for this class
    static boost::program_options::options_description& opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! calculates phenotypic values from genotype
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

    //! \f$\Xi(I, u, v)\f$ anolis_v2
    double habitat_preference_v2(const double height, const double diameter) const;

    //! \f$\Xi(I, u, v)\f$ anolis_v3
    double habitat_preference_v3(const double height, const double diameter) const;

    //! \f$\Xi(I, u, v)\f$
    /** @param height habitat environment
        @param diameter habitat environment
        @return \f$\Xi(I, u, v)\f$
    */
    double habitat_preference(const double height, const double diameter) const {
        return habitat_preference_v3(height, diameter);
    }

    //! \f$D_I\f$ numerical computation (slow)
    double denom_numerical() const;

    //! \f$D_I\f$ analytical computation by Mathematica (fast)
    double denom_mathematica() const;

    //! \f$D_I\f$ analytical computation by Maple
    //! @bug something wrong
    double denom_maple() const;

    //! \f$D_I\f$
    /** @return \f$D_I\f$
    */
    double denom_() const {return denom_mathematica();}

    //! \f$W(I, u, v)\f$
    /** @param height habitat environmant
        @param diameter habitat environment
    */
    double fitness(const double height, const double diameter) const;

    //! free recombination
    /** @param lhs left arm of a chromosome
        @param rhs right arm of a chromosome
        @return haplotype after recombination
    */
    static Loci recombination(const Loci& lhs, const Loci& rhs);

    //! poisson process with \f$\lambda\f$ = MU_LOCUS_ * NUM_LOCI_ * trait::size
    /** @param haplotype (call-by-value)
        @return mutated haplotype
    */
    static std::vector<Loci> mutate(std::vector<Loci> haplotype);

    //! unit test for the class
    friend void individual_unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member

    //! Unlinked diallelic loci with equal effects
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;

    //! All traits are scaled to be between 0 and 1
    std::vector<double> phenotype_;

    //! \f$D_I\f$ can be calculated at birth
    double denominator_;

    //! \f$K_e(I)\f$ can be calculated at birth
    double effective_carrying_capacity_;
};

inline std::ostream& operator<< (std::ostream& ost, const Individual& ind) {
    return ost << ind.str();
}
extern void individual_unit_test();

#endif /* INDIVIDUAL_H_ */
