// -*- mode: c++; coding: utf-8 -*-
/*! @file patch.h
    @brief Interface of Patch class
*/
#pragma once
#ifndef PATCH_H_
#define PATCH_H_

#include <vector>

#include "individual.h"

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Patch {
  public:

    //! Construct an empty patch
    Patch() = default;
    
    //! Construct a patch with the same number of females and males
    /*! @param n The number of inital individuals in this patch
    */
    Patch(const size_t n): females_(n / 2), males_(n - females_.size()) {}

    //! Add an individual to this patch
    /*! @param ind New individual to add
    */
    void append(Individual&&);

    //! @return The number of individuals in this patch
    size_t size() const {return females_.size() + males_.size();}

    //! @return Whether this patch is empty or not
    bool empty() const {return females_.empty() && males_.empty();}

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

    //! Unit test for Patch
    static void unit_test();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! Calculate \f$N_e(I)\f$ from \f$C(I,J)\f$ and \f$C'(I,J)\f$
    /*! @ingroup natural_selection
        @param focal individual
        @return \f$N_e\f$ of the focal individual
        \f[
            N_e(I) = \sum_J C(I,J) C'(I,J)
        \f]
    */
    double effective_num_competitors(const Individual&) const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member
    //! female individuals
    std::vector<Individual> females_;

    //! male individuals
    std::vector<Individual> males_;
};

//! Stream operator for Patch
extern std::ostream& operator<< (std::ostream& ost, const Patch& patch);

#endif /* PATCH_H_ */
