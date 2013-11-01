// -*- mode: c++; coding: utf-8 -*-
#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <vector>
#include <numeric>

#include "individual.h"
#include "environment.h"

void constant() {
    std::vector<size_t> trait_values(trait::size * 2 - 1);
    std::iota(trait_values.begin(), trait_values.end(), 0);
    double kii = 0;
    double kii_prime = 0;
    for (const auto x0: trait_values) {
    for (const auto x1: trait_values) {
    for (const auto x2: trait_values) {
    for (const auto x3: trait_values) {
        Individual ind(std::vector<size_t>{x0, x1, x2, x3});
        std::cout << ind.carrying_capacity() << " ";
        std::cout << ind.carrying_capacity_prime() << std::endl;
        kii += ind.carrying_capacity();
        kii_prime += ind.carrying_capacity_prime();
    }
    }
    }
    }
    std::cout << kii << std::endl;
    std::cout << kii_prime << std::endl;
}

void run() {
    individual_unit_test();
    constant();
}

#endif /* MAIN_H_ */
