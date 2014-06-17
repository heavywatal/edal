// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.h
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
#include "cxxwtils/prandom.hpp"

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

/*! @brief sexual, diploid, additive, unlinked, diallelic

    Individuals are sexual and diploid.
    Each individual has a number of additive quantitative characters:
    - two “morphological” characters
      \f$ x_0 \f$ (toepad size) and \f$ x_1 \f$ (limb length) that control viability;
    - two “habitat preference” characters:
      the most preferred height \f$ y_0 \f$ and the most preferred diameter \f$ y_1 \f$
    - three “mating compatibility” characters \f$ m \f$, \f$ f \f$, and \f$ c \f$.

    The male display trait \f$ m \f$ is expressed in males only,
    whereas female mating preference \f$ f \f$ and tolerance \f$ c \f$ are expressed in females only.
    Other traits are expressed in both sexes.
    All traits are scaled to be between 0 and 1,
    and are controlled by different unlinked diallelic loci with equal effects.
    Mutations occur at equal rates across all loci;
    the probabilities of forward and backward mutations are equal.
*/
class Individual {
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /*! @addtogroup biol_param
        @{*/

    //! \f$ \alpha \f$ of Beta distribution in \f$ F(u,v) \f$
    static double BETA_PARAM_;

    //! Coefficient for \f$ K_e(I)\f$ --- **NOT FOUND in anolis_v3.pdf**
    static size_t CARRYING_CAPACITY_;

    //! \f$ b \f$ in \f$ w(I) \f$
    static size_t AVG_NUM_OFFSPINRGS_;

    //! \f$ h_0 \f$ in \f$ \Xi(y_0,y_1|u,v) \f$
    static double HEIGHT_PREFERENCE_;

    //! \f$ h_1 \f$ in \f$ \Xi(y_0,y_1|u,v) \f$
    static double DIAMETER_PREFERENCE_;

    //! \f$ c_0 \f$ in \f$ C(I,J) \f$
    static double HEIGHT_COMPETITION_;

    //! \f$ c_1 \f$ in \f$ C(I,J) \f$
    static double DIAMETER_COMPETITION_;

    //! \f$ s_0 \f$ in \f$ W(x_0,x_1|u,v) \f$
    static double TOEPAD_SELECTION_;

    //! \f$ s_1\f$ in \f$ W(x_0,x_1|u,v) \f$
    static double LIMB_SELECTION_;

    //! \f$ \sigma_a \f$ in \f$ \Psi(f,c|m) \f$
    static double MATING_SIGMA_;

    //! Mutation rate per locus per generation
    /*! Mutations occur at equal rates across all loci;
        the probabilities of forward and backward mutations are equal.
    */
    static double MU_LOCUS_;

    //! Unused yet
    static double MU_NEUTRAL_;

    //! Migration rate \f$ m \f$ per generation
    static double MIGRATION_RATE_;

    //! The number of loci per trait
    constexpr static size_t NUM_LOCI_ = 8;

    //! Genotype that produce trait value = 1.0, i.e., `11111111`
    constexpr static unsigned long FULL_BITS = wtl::pow<NUM_LOCI_>(2) - 1;

    //! Genotype that produce trait value = 0.5, i.e., `00001111`
    constexpr static unsigned long HALF_BITS = wtl::pow<NUM_LOCI_ / 2>(2) - 1;

    /** @} endgroup */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! Compile-time constant value used in init_phenotype()
    constexpr static double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;

    //! Precision of numerical integration
    constexpr static size_t NUM_STEPS_ = 32;

    //! typedef for diallelic loci of a trait
    typedef std::bitset<NUM_LOCI_> Loci;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:

    //! Default constructor for original individuals
    Individual(): genotype_{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS}},
        phenotype_(init_phenotype()),
        denominator_{calc_denom()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    //! Constructor for sexual reproduction
    Individual(const std::vector<Loci>& egg, const std::vector<Loci>& sperm):
        genotype_{egg, sperm},
        phenotype_(init_phenotype()),
        denominator_{calc_denom()},
        effective_carrying_capacity_{effective_carrying_capacity()} {}

    //! Initialization by phenotypic values
    Individual(const std::vector<size_t>&);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /** @addtogroup biol_proc
        @{*/

    //! Competition coefficient \f$C(I, J)\f$ in anolis_v3
    /*! @ingroup habitat_pareference
        @param other individual to interact
        @return \f$C(I, J)\f$
        @retval 1 for individuals with identical preferences

        following Roughgarden and others
    */
    double habitat_overlap(const Individual& other) const;

    //! \f$C(I, J)\f$ in anolis_v2
    /*! @ingroup habitat_pareference
    */
    double habitat_overlap_v2(const Individual&) const;

    //! \f$K_e(I)\f$
    /*! @ingroup natural_selection
    */
    double effective_carrying_capacity() const;

    //! Probability of survival \f$w(I)\f$
    /*! @ingroup natural_selection
        @param effective_num_competitors \f$ K_e(I) \f$
        @retval 1.0 if \f$ N_e(I) \ll K_e(I) \f$
        @retval 1/b if \f$ N_e(I) = K_e(I) \f$
        @see AVG_NUM_OFFSPINRGS_

        The probability that an individual survives to the age of reproduction is
        \f[
            w(I) = \frac 1 {1 + (b-1) \frac {N_e(I)^\alpha} {K_e(I)}}
        \f]
        where the parameter \f$ \alpha \f$ (= 1 for now) controls the strength of crowding,
        and \f$ b > 0 \f$ (Individual::AVG_NUM_OFFSPINRGS_) is a parameter
        (average number of offspring per female; see below).
        This is the Beverton-Holt model which represents
        a discrete-time analog of the logistic model
    */
    double survival_probability(const double effective_num_competitors) const;

    //! \f$ P(I, I') = \Xi(I, I') C(I, I')\f$
    /*! @ingroup mating
    */
    double mating_probability(const Individual& male) const {
        return mating_preference(male) * habitat_overlap(male);
    };

    //! generates poisson random number with \f$\lambda\f$ = Individual::AVG_NUM_OFFSPINRGS_
    /*! @ingroup mating
        @return the number of offsprings
    */
    size_t poisson_offsprings() const;

    //! Gametogenesis with free recombination among loci
    /*! @ingroup mating
        @return a gamete
    */
    std::vector<Loci> gametogenesis() const;

    //! Bernoulli trial whether migrating or not
    /*! @ingroup life_cycle
    */
    bool is_migrating() const {
        return prandom().bernoulli(MIGRATION_RATE_);
    }

    /** @} biol_proc */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! CSV formated string
    /*! @return CSV formated string
    */
    std::string str() const;

    //! The header correspongs to str()
    static std::string header();

    static boost::program_options::options_description& opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /** @addtogroup biol_proc
        @{*/

    //! \f$\Xi(I, u, v)\f$ quadratic approximation in anolis_v3
    /*! @ingroup habitat_pareference
        @param height habitat environment
        @param diameter habitat environment
        @return \f$\Xi(I, u, v)\f$
    */
    double habitat_preference(const double height, const double diameter) const;

    //! \f$\Xi(I, u, v)\f$ in anolis_v2
    /*! @ingroup habitat_pareference
    */
    double habitat_preference_v2(const double height, const double diameter) const;

    //! \f$D_I\f$ numerical computation (slow)
    /*! @ingroup habitat_pareference
    */
    double calc_denom_numerical() const;

    //! \f$D_I\f$ analytical computation by Mathematica (fast)
    /*! @ingroup habitat_pareference
    */
    double calc_denom_mathematica() const;

    //! \f$D_I\f$ analytical computation by Maple
    /*! @ingroup habitat_pareference
        @bug something wrong (unused)
    */
    double calc_denom_maple() const;

    //! \f$D_I\f$
    /*! @ingroup habitat_pareference
        @return \f$D_I\f$
    */
    double calc_denom() const {return calc_denom_mathematica();}

    //! \f$\Psi(I, I')\f$
    /*! @ingroup mating
    */
    double mating_preference(const Individual& male) const;

    //! Measure adaptation to habitat \f$W(x_0, x_1| u, v)\f$
    /*! @ingroup natural_selection
        @param height habitat environmant
        @param diameter habitat environment

        The optimum values of traits \f$ x_0 \f$ and \f$ x_1 \f$
        under environmental conditions \f$(u, v)\f$ are
        \f$x_0 = u\f$ and \f$ x_1 = v \f$
    */
    double fitness(const double height, const double diameter) const;

    //! Free recombination among unlinked loci
    /*! @ingroup mating
        @param lhs left arm of a chromosome
        @param rhs right arm of a chromosome
        @return haplotype after recombination
    */
    static Loci recombination(const Loci& lhs, const Loci& rhs);

    //! poisson process with \f$\lambda\f$ = MU_LOCUS_ * NUM_LOCI_ * trait::size
    /*! @ingroup mating
        @param haplotype (call-by-value)
        @return mutated haplotype
    */
    static std::vector<Loci> mutate(std::vector<Loci> haplotype);

    /** @} biol_proc */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! calculates phenotypic values from genotype
    /*! @return phenotypic values
    */
    std::vector<double> init_phenotype() const {
        constexpr size_t n = trait::size;
        std::vector<double> output(n);
        for (size_t i=0; i<n; ++i) {
            output[i] = (genotype_.first[i].count() + genotype_.second[i].count()) * INV_NUM_LOCI_;
        }
        return output;
    };

    friend double pdf_beta(const double height, const double diameter);
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

//! Stream operator for Individual
inline std::ostream& operator<< (std::ostream& ost, const Individual& ind) {
    return ost << ind.str();
}

//! Unit test for Individual
extern void individual_unit_test();

#endif /* INDIVIDUAL_H_ */
