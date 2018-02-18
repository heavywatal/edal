/*! @file patch.hpp
    @brief Interface of Patch class
*/
#pragma once
#ifndef PATCH_HPP_
#define PATCH_HPP_

#include <iosfwd>
#include <vector>
#include <map>
#include <memory>

namespace wtl {class sfmt19937_64;}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

class Individual;

class Patch {
  public:
    using URBG = wtl::sfmt19937_64;

    //! Construct an empty patch
    Patch(unsigned int seed);

    //! Construct a patch with the same number of females and males
    /*! @param n The number of inital individuals in this patch
    */
    Patch(size_t n, unsigned int seed);

    //! Construct a patch with a non-default Individual
    /*! @param n The number of inital individuals in this patch
        @param founder The individual to be copied
    */
    Patch(size_t n, const Individual& founder, unsigned int seed);

    //! Add an individual to this patch
    /*! @param ind New individual to add
    */
    void emplace_back(Individual&& ind) {members_.emplace_back(std::forward<Individual>(ind));}

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

    /*! @brief Change row/col with probability \f$m\f$ = Individual::MIGRATION_RATE_

        > With probability \f$ m > 0 \f$, each offspring becomes a "migrant."
        > Each migrant goes to one of the 8 neighboring patches.
        > For patches at the boundary,
        > the probability \f$ m \f$ is reduced according to the number of neighbors they have.
    */
    std::vector<std::pair<unsigned, unsigned int>>
    make_destinations(size_t n, size_t row, size_t col, size_t num_rows, size_t num_cols) const;

    //! Some individuals die depending on Individual::survival_probability()
    /*! @ingroup natural_selection
    */
    void viability_selection();

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

    std::unique_ptr<URBG> engine_;
};

//! Stream operator for Patch
extern std::ostream& operator<< (std::ostream& ost, const Patch& patch);

} // namespace edal

#endif /* PATCH_HPP_ */
