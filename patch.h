// -*- mode: c++; coding: utf-8 -*-
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

    Patch() = default;
    Patch(const size_t n): females_(n / 2), males_(n - females_.size()) {}

    void append(const Individual&);
    size_t size() const {return females_.size() + males_.size();}
    std::string str() const;

    std::vector<Individual> mate_and_reproduce() const;
    void viability_selection();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    double effective_num_competitors_f(const size_t index) const;
    double effective_num_competitors_m(const size_t index) const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member
    std::vector<Individual> females_;
    std::vector<Individual> males_;
};

inline std::ostream& operator<< (std::ostream& ost, const Patch& patch) {
    return ost << patch.str();
}
extern void patch_unit_test();

#endif /* PATCH_H_ */
