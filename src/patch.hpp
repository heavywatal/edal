/*! @file patch.hpp
    @brief Interface of Patch class
*/
#pragma once
#ifndef PATCH_HPP_
#define PATCH_HPP_

#include <vector>

#include "individual.hpp"

namespace wtl {
    class sfmt19937;
}

class Patch {
  public:
    using URNG = wtl::sfmt19937;

    //! Construct an empty patch
    Patch() = default;

    //! Construct a patch with the same number of females and males
    /*! @param n The number of inital individuals in this patch
    */
    Patch(const size_t n): members_(n) {
        change_sex_half(n);
    }

    //! Construct a patch with a non-default Individual
    /*! @param n The number of inital individuals in this patch
        @param founder The individual to be copied
    */
    Patch(const size_t n, const Individual& founder):
        members_(n, founder) {
        change_sex_half(n);
    }

    //! Copy constructor
    Patch(const Patch& obj):
        members_{obj.members_} {}

    //! Add an individual to this patch
    /*! @param ind New individual to add
    */
    void append(Individual&& ind) {members_.push_back(std::move(ind));}

    //! @return The number of individuals in this patch
    size_t size() const {return members_.size();}

    //! @return Whether this patch is empty or not
    bool empty() const {return members_.empty();}

    //! Count genotypes
    std::map<Individual, size_t> summarize() const;

    //! All females mate with someone according to mating probability
    /*! @ingroup mating
        @return offsprings of all the females in the patch

        Each mating results in a number of offspring
        drawn from a Poisson distribution with parameter b.
        We assume that all adult females mate.
        This assumption implies that any costs of mate choice,
        which can easily prevent divergence and speciation, are absent.
        This assumption also means that the effective population size is
        increased relative to the actual number of adults.
    */
    std::vector<Individual> mate_and_reproduce(URNG&) const;

    //! Some individuals die depending on Individual::survival_probability()
    /*! @ingroup natural_selection
    */
    void viability_selection(URNG&);

    //! Unit test for Patch
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! Calculate \f$N_e(I)\f$ from \f$C(I,J)\f$
    /*! @ingroup natural_selection
        @return \f$N_e\f$ of the focal individual
        \f[
            N_e(I) = \sum_J C(I,J)
        \f]
    */
    std::vector<double> effective_num_competitors() const;

    //! Change sex of second half members
    void change_sex_half(size_t n);

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member
    //! Individuals
    std::vector<Individual> members_;
};

//! Stream operator for Patch
extern std::ostream& operator<< (std::ostream& ost, const Patch& patch);

#endif /* PATCH_HPP_ */
