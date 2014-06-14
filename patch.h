// -*- mode: c++; coding: utf-8 -*-
/** @file patch.h
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

    //! construct an empty patch
    Patch() = default;
    
    //! construct a patch with the same number of females and males
    /** @param n number of inital individuals in this patch
    */
    Patch(const size_t n): females_(n / 2), males_(n - females_.size()) {}

    //! add an individual to this patch
    /** @param ind new individual to add
    */
    void append(const Individual&);

    //! @return the number of individuals in this patch
    size_t size() const {return females_.size() + males_.size();}

    //! @return whether this patch is empty or not
    bool empty() const {return females_.empty() && males_.empty();}

    //! @return string representation for debugging
    std::string str() const;

    //! all the females reproduce
    /** @return offsprings of the females in this patch
    */
    std::vector<Individual> mate_and_reproduce() const;

    //! Some individuals die depending on \f$N_e(I)\f$
    void viability_selection();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    //! calculate \f$N_e(I)\f$
    //! @param focal individual
    //! @return \f$N_e\f$ of the focal individual
    double effective_num_competitors(const Individual&) const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member
    //! female individuals
    std::vector<Individual> females_;

    //! male individuals
    std::vector<Individual> males_;
};

inline std::ostream& operator<< (std::ostream& ost, const Patch& patch) {
    return ost << patch.str();
}
extern void patch_unit_test();

#endif /* PATCH_H_ */
