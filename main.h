// -*- mode: c++; coding: utf-8 -*-
#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <sstream>
#include <vector>

#include "cxxwtils/prandom.hpp"

#include "individual.h"
#include "patch.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// parameters

size_t num_rows = 8;
size_t num_cols = 8;
size_t initial_patch_size = 20;
double migration_rate = 1e-1;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// data members

std::vector<std::vector<Patch> > matrix;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions

void init() {
    matrix.assign(num_rows, std::vector<Patch>(num_cols));
    matrix[0][0] = Patch(initial_patch_size);
}

std::pair<size_t, size_t> migrate(const size_t row_orig, const size_t col_orig) {
    const size_t x = prandom().randrange(8);
    size_t row = row_orig;
    size_t col = col_orig;
    switch (x) {
      case 0:        ++col; break;
      case 1: ++row; ++col; break;
      case 2: ++row;        break;
      case 3: ++row; --col; break;
      case 4:        --col; break;
      case 5: --row; --col; break;
      case 6: --row;        break;
      case 7: --row; ++col; break;
    }
    if (row > 7 | col > 7) {row = row_orig; col = col_orig;}
    return {row, col};
}

void life_cycle() {
    std::vector<std::vector<Patch> > next_generation(matrix);
    for (size_t row=0; row<num_rows; ++row) {
        for (size_t col=0; col<num_cols; ++col) {
            auto offsprings = matrix[row][col].mate_and_reproduce();
            for (const auto& child: offsprings) {
                if (prandom().bernoulli(migration_rate)) {
                    auto new_coords = migrate(row, col);
                    next_generation[new_coords.first][new_coords.second].append(child);
                } else {
                    next_generation[row][col].append(child);
                }
            }
        }
    }
    for (auto& row: next_generation) {
        for (auto& patch: row) {
            patch.viability_selection();
        }
    }
    matrix.swap(next_generation);
}

template <class Vec2D, class Func> inline
void print_matrix(const Vec2D& m, Func func) {
    std::ostringstream ost;
    for (const auto& row: m) {
        for (const auto& cell: row) {
            ost << func(cell) << " ";
        }
        ost << "\n";
    }
    std::cout << ost.str();
}

void test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    individual_unit_test();
    patch_unit_test();
}

void run() {
    test();
    init();
    for (size_t i=0; i<30; ++i) {
        life_cycle();
        print_matrix(matrix, [](const Patch& p) {return p.size();});
        std::cout << std::endl;
    }
}

#endif /* MAIN_H_ */
