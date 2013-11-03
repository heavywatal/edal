// -*- mode: c++; coding: utf-8 -*-
#ifndef PATCH_H_
#define PATCH_H_

#include <vector>

#include "individual.h"

class Patch {
  public:

    Patch() = default;
    Patch(const size_t n): females_(n / 2), males_(n / 2) {}

    void append(const Individual&);
    size_t size() const {return females_.size() + males_.size();}

    std::vector<Individual> mate_and_reproduce() const;
    void viability_selection();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:

    double effective_population_size() const;

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data member
    std::vector<Individual> females_;
    std::vector<Individual> males_;
};

extern void patch_unit_test();

#endif /* PATCH_H_ */
