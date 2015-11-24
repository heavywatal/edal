// -*- mode: c++; coding: utf-8 -*-
/*! @file patch.h
    @brief Interface of Patch class
*/
#pragma once
#ifndef PATCH_H_
#define PATCH_H_

#include <vector>

#include <cxxwtils/prandom.hpp>

#include "individual.h"

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Patch {
  public:

    //! Construct an empty patch
    Patch(): rng_{wtl::sfmt()()} {};
    
    //! Construct a patch with the same number of females and males
    /*! @param n The number of inital individuals in this patch
    */
    Patch(const size_t n): members_(n), rng_{wtl::sfmt()()} {
        change_sex_half(n);
    }

    //! Construct a patch with a non-default Individual
    /*! @param n The number of inital individuals in this patch
        @param founder The individual to be copied
    */
    Patch(const size_t n, const Individual& founder):
        members_(n, founder), rng_{wtl::sfmt()()} {
        change_sex_half(n);
    }

    //! Copy constructor
    Patch(const Patch& obj):
        members_{obj.members_}, rng_{wtl::sfmt()()} {}

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
    std::vector<Individual> mate_and_reproduce() const;

    //! Some individuals die depending on Individual::survival_probability()
    /*! @ingroup natural_selection
    */
    void viability_selection();


    /*! @brief Change row/col with probability \f$m\f$ = Individual::MIGRATION_RATE_

        > With probability \f$ m > 0 \f$, each offspring becomes a "migrant."
        > Each migrant goes to one of the 8 neighboring patches.
        > For patches at the boundary,
        > the probability \f$ m \f$ is reduced according to the number of neighbors they have.
    */
    std::pair<size_t, size_t> choose_patch(size_t row, size_t col) const;

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

    //! Random number generator
    mutable wtl::sfmt19937 rng_;
};

//! Stream operator for Patch
extern std::ostream& operator<< (std::ostream& ost, const Patch& patch);

#endif /* PATCH_H_ */
