/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef INDIVIDUAL_HPP_
#define INDIVIDUAL_HPP_
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>

#include <wtl/math.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace wtl {
    class sfmt19937;
}

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

    //! \f$ K_{\max} \f$ in \f$ K_e(I)\f$
    static size_t CARRYING_CAPACITY_;

    //! \f$ b \f$ in \f$ w(I) \f$
    static double AVG_NUM_OFFSPINRGS_;

    //! \f$ h_0 \f$ in \f$ \Xi(y_0,y_1 \mid u,v) \f$
    static double HEIGHT_PREFERENCE_;

    //! \f$ h_1 \f$ in \f$ \Xi(y_0,y_1 \mid u,v) \f$
    static double DIAMETER_PREFERENCE_;

    //! \f$ s_0 \f$ in \f$ W(x_0,x_1 \mid u,v) \f$
    static double TOEPAD_SELECTION_;

    //! \f$ s_1\f$ in \f$ W(x_0,x_1 \mid u,v) \f$
    static double LIMB_SELECTION_;

    //! \f$c_y\f$ in \f$C_y(I,J)\f$
    static double PREF_COMPETITION_;

    //! \f$c_x\f$ in \f$C_x(I,J)\f$
    static double MORPH_COMPETITION_;

    //! \f$ \sigma_a \f$ in \f$ \psi(f,c \mid m) \f$
    static double MATING_SIGMA_;

    //! Mutation rate per locus per generation
    /*! Mutations occur at equal rates across all loci;
        the probabilities of forward and backward mutations are equal.
    */
    static double MU_LOCUS_;

    //! Flag set to protect specific traits from mutations
    /*! e.g., 10 ("00001010") blocks mutations on limb and diameter preference
    */
    static unsigned long MUTATION_MASK_;

    //! Migration rate \f$m\f$ (i.e., \f$Nm\f$ makes the expected # of migrants)
    static double MIGRATION_RATE_;

    //! The number of loci per trait
    constexpr static size_t NUM_LOCI_ = 8;

    /** @} endgroup */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! \f$ \sigma \f$ of Normal distribution in v2 \f$ F(u,v) \f$
    static double NORMAL_SIGMA_;

    //! \f$ c_0 \f$ of Normal distribution in v2 \f$ F(u,v) \f$
    static double C0_;

    //! \f$ c_1 \f$ of Normal distribution in v2 \f$ F(u,v) \f$
    static double C1_;

    //! Genotype that produce trait value = 1.0, i.e., `11111111`
    constexpr static unsigned long FULL_BITS = wtl::pow(2, NUM_LOCI_) - 1;

    //! Genotype that produce trait value = 0.5, i.e., `00001111`
    constexpr static unsigned long HALF_BITS = wtl::pow(2, NUM_LOCI_ / 2) - 1;

    //! Compile-time constant value used in init_phenotype()
    constexpr static double INV_NUM_LOCI_ = 0.5 / NUM_LOCI_;

    //! Precision of numerical integration
    constexpr static size_t NUM_STEPS_ = 32;

    //! Store \f$K_e(I)\f$ with ecological traits as keys
    static std::map<std::vector<double>, double> KE_CACHE_;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  public:
    //! Uniform Random Number Generator
    using URNG = wtl::sfmt19937;

    //! typedef for diallelic loci of a trait
    typedef std::bitset<NUM_LOCI_> Loci;

    //! Default constructor for original individuals
    Individual(): Individual{
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS},
        {HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS, HALF_BITS}} {}

    //! Constructor for sexual reproduction
    Individual(const std::vector<Loci>& egg,
               const std::vector<Loci>& sperm,
               const bool is_male=false):
        genotype_{egg, sperm},
        phenotype_(calc_phenotype()),
        sex_(is_male),
        ke_(effective_carrying_capacity_cache()) {}

    //! Homozygous initialization by bit values
    Individual(const std::vector<unsigned long>&);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /** @addtogroup biol_proc
        @{*/

    //! \f$K_e(I)\f$ with quadratic \f$\Xi\f$ before normalization (v3u)
    /*! @ingroup natural_selection
        \f[
            K_e(I) = K_0 \int_0^1 \int_0^{1-u}
                F(u,v) W(x_0,x_1 \mid u,v) \Xi(y_0,y_1 \mid u,v) dv du
        \f]
    */
    double effective_carrying_capacity_quad_unnormalized() const;

    //! \f$K_e(I)\f$ with exponential \f$\Xi\f$
    /*! @ingroup natural_selection
    */
    double effective_carrying_capacity_exp_unnormalized() const;

    //! \f$K_e(I)\f$ with old resource distribution and exponential \f$\Xi\f$ before normalization
    /*! @ingroup natural_selection
    */
    double effective_carrying_capacity_old_exp_unnormalized() const;

    //! Find \f$K_e(I)\f$ in KE_CACHE_ or calculate
    double effective_carrying_capacity_cache() const;

    //! Exponent of \f$C_y(I, J)\f$: competition on habitat preference
    /*! @ingroup habitat_pareference
        @param other individual to interact
        @return \f$C_y(I, J)\f$
        @retval 1 for individuals with identical preferences

        following Roughgarden and others
        \f[
            C_y(I,J) = \exp(-\frac{(y_{0,I} - y_{0,J})^2}{2c_y^2}
                            -\frac{(y_{1,I} - y_{1,J})^2}{2c_y^2})
        \f]
    */
    double preference_overlap(const Individual& other) const;

    //! Exponent of \f$C_x(I, J)\f$: competition on morphology
    /*! @ingroup habitat_pareference
        @param other individual to interact
        @return \f$C_x(I, J)\f$
        @retval 1 for individuals with identical preferences

        \f[
            C_x(I,J) = \exp(-\frac{(x_{0,I} - x_{0,J})^2}{2c_x^2}
                            -\frac{(x_{1,I} - x_{1,J})^2}{2c_x^2})
        \f]
    */
    double morphology_overlap(const Individual& other) const;

    //! \f$C(I, J) = C_x(I,J) C_y(I,J)\f$ for competition
    /*! @ingroup habitat_pareference
    */
    double resource_overlap(const Individual& other) const {
        return std::exp(preference_overlap(other) + morphology_overlap(other));
    }

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

    //! \f$P(I,I') = \psi(I,I') C_y(I,I')\f$
    /*! @ingroup mating
        \f[
            \psi(f,c\mid m) = \left\{
              \begin{array}{ll}
                \exp \left( -(2c-1)^2 \frac{(f-m)^2}{2\sigma_a^2}\right)
                  & \mbox{if}\ c > 0.5,\\
                1 & \mbox{if}\ c=0.5,\\
                \exp \left( -(2c-1)^2 \frac{(f-(1-m))^2}{2\sigma_a^2}\right)
                  & \mbox{if}\ c<0.5,
            \end{array} \right.
        \f]
    */
    double mating_probability(const Individual& male) const;

    //! \f$P(I,I') = \psi(I,I') C_y(I,I')\f$ with Debarre 2012
    /*! @ingroup mating
        \f[
            \psi(f,c\mid m) = \left\{
              \begin{array}{ll}
                1 - (2c-1)^2\Big[1 - \exp\Big(-\frac {(f-m)^2}{2\sigma_a}\Big)\Big]
                  & \mbox{if}\ c > 0.5,\\
                1 & \mbox{if}\ c=0.5,\\
                1 - (2c-1)^2\Big[    \exp\Big(-\frac {(f-m)^2}{2\sigma_a}\Big)\Big]
                  & \mbox{if}\ c<0.5,
            \end{array} \right.
        \f]
    */
    double mating_probability_debarre(const Individual& male) const;

    //! \f$P(I,I') = \psi(I,I') C_y(I,I')\f$ with Thibert-Plante and Gavrilets 2013
    /*! @ingroup mating
        \f[
            \psi(f,c\mid m) = \left\{
              \begin{array}{ll}
                \exp \left( -(2c-1)^2 \frac{(f-m)^2}{2\sigma_a^2}\right)
                  & \mbox{if}\ c > 0.5,\\
                1 & \mbox{if}\ c=0.5,\\
                2 - \exp \left( -(2c-1)^2 \frac{(f-m)^2}{2\sigma_a^2}\right)
                  & \mbox{if}\ c<0.5,
            \end{array} \right.
        \f]
    */
    double mating_probability_TPG2013(const Individual& male) const;

    //! Gametogenesis with free recombination and mutation
    /*! @ingroup mating
        @return a gamete
    */
    std::vector<Loci> gametogenesis(URNG&) const;

    //! Getter
    static double MIGRATION_RATE() {return MIGRATION_RATE_;}

    //! Getter
    static size_t AVG_NUM_OFFSPINRGS() {return AVG_NUM_OFFSPINRGS_;}

    /** @} biol_proc */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! Get sex_
    bool is_male() const {return sex_;}

    //! Change sex_
    void change_sex() {sex_ = !sex_;}

    //! Operator required to be std::map key
    bool operator<(const Individual& other) const {
        return genotype_ < other.genotype_;
    }

    //! The header correspongs to str()
    static std::string header();

    //! detailed str() for testing/debugging
    std::string str_detail() const;

    friend std::ostream& operator<< (std::ostream& ost, const Individual& ind);

    //! test function
    static void unit_test();
    //! test function
    static void write_resource_abundance();
    //! test function
    static std::string possible_phenotypes();
    //! test function
    static std::string possible_geographic();

    static boost::program_options::options_description opt_description();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /** @addtogroup biol_proc
        @{*/

    //! Exponent of \f$\Xi(I, u, v)\f$ in anolis_v2
    /*! @ingroup habitat_pareference
        @param height habitat environment
        @param diameter habitat environment
        @return \f$\Xi(I, u, v)\f$ before normalization
        \f[
            \Xi(I,u,v) \propto \exp(-\frac{(u - y_0)^2}{2h_0^2}
                                    -\frac{(v - y_1)^2}{2h_1^2})
        \f]
    */
    double habitat_preference_exp(const double height, const double diameter) const;

    //! \f$\Xi(I, u, v)\f$ quadratic approximation in anolis_v3
    /*! @ingroup habitat_pareference
        @param height habitat environment
        @param diameter habitat environment
        @return \f$\Xi(I, u, v)\f$ before normalization
        \f[
            \Xi(I,u,v) \propto 1 - h_0 (u - y_0)^2 - h_1 (v - y_1)^2
        \f]
    */
    double habitat_preference_quadratic(const double height, const double diameter) const;

    //! \f$D_I\f$ analytical computation by Mathematica (fast but limited)
    /*! @ingroup habitat_pareference
        @return \f$D_I\f$
        @bug not applicable in the parameter range where
        \f$\Xi\f$ goes negative with high \f$h_0\f$ and \f$h_1\f$.
        \f[
            D_I = \frac {12
                - 6 h_0 (1 + 2 (-1 + y_0) y_0)
                + h_1 (-1 + 4 (1 - 3 y_1) y_1)
                - \alpha (-24 + h_1 + 6 h_0 (1 - 2 y_0)^2 - 8 h_1 y_1 + 24 h_1 y_1^2)
            } {12 + 24 \alpha}
        \f]
    */
    double calc_DI_analytical() const;

    //! \f$D_I\f$ numerical computation (slow)
    /*! @ingroup habitat_pareference
        @return \f$D_I\f$
        \f[
            D_I = \int_0^1 \int_0^{1-u} \Xi(y_0,y_1 \mid u,v) F(u,v) dv du
        \f]
    */
    double calc_DI_numerical() const;

    //! Calculate quadratic \f$\Xi(I, u, v)\f$ normalizer with analytical solution
    /*! @ingroup habitat_pareference
        @bug not applicable in the parameter range where
        \f$\xi\f$ goes negative with high \f$h_0\f$ and \f$h_1\f$.
        \f[
            \frac{1} 2
            -h_0 (\frac{y_0^2} 2 + \frac{y_0} 3 - \frac{1} {12})
            -h_1 (\frac{y_1^2} 2 + \frac{y_1} 3 - \frac{1} {12})
        \f]
    */
    double calc_Dxi_analytical() const;

    //! Numerical integration of quadratic \f$\xi(I, u, v)\f$
    /*! @ingroup habitat_pareference
        \f[
            \int_0^1 \int_0^{1-u} \xi(y_0,y_1 \mid u,v) dv du =
            \int_0^1 \int_0^{1-u} \{1 - h_0 (u - y_0)^2 - h_1 (v - y_1)^2\} dv du
        \f]
    */
    double calc_Dxi_numerical() const;

    //! Exponent of \f$W(x_0, x_1 \mid u, v)\f$: measure adaptation to habitat
    /*! @ingroup natural_selection
        @param height habitat environmant
        @param diameter habitat environment

        \f[
            W(x_0,x_1 \mid u,v) = \exp(-\frac{(x_0 - u)^2}{2s_0^2}
                                  -\frac{(x_1 - v)^2}{2s_1^2})
        \f]
    */
    double fitness(const double height, const double diameter) const;

    /** @} biol_proc */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! calculates phenotypic values from genotype
    /*! @return phenotypic values
    */
    std::vector<double> calc_phenotype() const {
        constexpr size_t n = trait::size;
        std::vector<double> output(n);
        for (size_t i=0; i<n; ++i) {
            output[i] = (genotype_.first[i].count() + genotype_.second[i].count()) * INV_NUM_LOCI_;
        }
        return output;
    };

    //! aggregation of intermediate phenotype values
    std::vector<double> intermediate_phenotypes() const;

    //! labels for intermediate phenotypes
    static const std::vector<std::string> INTERMEDIATE_KEYS_;

    friend double pdf_beta(const double height, const double diameter);
    friend double pdf_normal(const double height, const double diameter);
    friend double pdf_exp(const double height, const double diameter);
    template <class Func> friend double integrate_triangle(Func&& func);
    template <class Func> friend double integrate_square(Func&& func);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member

    //! Unlinked diallelic loci with equal effects
    std::pair<std::vector<Loci>, std::vector<Loci> > genotype_;

    //! All traits are scaled to be between 0 and 1
    std::vector<double> phenotype_;

    //! female=0, male=1
    bool sex_ = false;

    //! \f$K_e(I)\f$ calculated in advance
    double ke_;
};

//! Overload: output 15 instead of 00001111
inline std::ostream& operator<< (std::ostream& ost, const Individual::Loci& bs) {
    return ost << bs.to_ulong();
}

namespace std {
//! Less operator for genotype comparison
inline bool operator< (const Individual::Loci& lhs, const Individual::Loci& rhs) {
    return lhs.to_ulong() < rhs.to_ulong();
}
}

#endif /* INDIVIDUAL_HPP_ */
