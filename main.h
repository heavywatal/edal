// -*- mode: c++; coding: utf-8 -*-
#pragma once
#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/getopt.hpp"
#include "cxxwtils/prandom.hpp"

#include "individual.h"
#include "patch.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// parameters

size_t num_rows = 8;
size_t num_cols = 8;
size_t initial_patch_size = 20;
double migration_rate = 1e-1;

unsigned int seed_ = std::random_device{}();

const char* HOME_ = std::getenv("HOME");
const fs::path HOME_DIR{HOME_};
const fs::path TMP_DIR{HOME_DIR / "tmp"};
fs::path WORK_DIR;
fs::path TOP_DIR;
fs::path LOAD_DIR;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// data members

std::vector<std::vector<Patch> > matrix;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions

//void check_flags(int argc, char* argv[]) {
//    namespace po = boost::program_options;
//    po::options_description description("Operation");
//    description.add_options()
//            ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
//            ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "verbose output")
//            ("test", po::value<bool>()->default_value(false)->implicit_value(true))
//            ("top_dir", po::value<std::string>()
//                ->default_value((TMP_DIR / wtl::strftime("out%Y%m%d")).string()))
//            ("row", po::value<size_t>(&num_rows)->default_value(num_rows))
//            ("col", po::value<size_t>(&num_cols)->default_value(num_cols))
//            ("migration_rate,m", po::value<double>(&migration_rate)->default_value(migration_rate))
//            ("seed", po::value<unsigned int>(&seed_)->default_value(seed_))
//            ;
//    description.add(Individual::opt_description());
//    description.add(Patch::opt_description());
//    po::variables_map vm;
//    po::store(po::parse_command_line(argc, argv, description), vm);
//    if (vm["help"].as<bool>()) {
//        description.print(std::cout);
//        exit(0);
//    }
//    TOP_DIR = fs::path(vm["top_dir"].as<std::string>());
//    prandom().seed(seed_); // TODO: want to read seed?
//    vm.notify();
//    if (vm["verbose"].as<bool>()) {
//        std::cout << description << std::endl;
//        std::cout << flags_into_string(description, vm) << std::endl;
//    }
//}

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
    const size_t n = 20;
    for (size_t i=0; i<n; ++i) {
        std::cout << "T = " << i << std::endl;
        print_matrix(matrix, [](const Patch& p) {return p.size();});
        std::cout << std::endl;
        life_cycle();
    }
    std::cout << "T = " << n << std::endl;
    print_matrix(matrix, [](const Patch& p) {return p.size();});
    std::cout << std::endl;
    for (const auto& row: matrix) {
        for (const auto& patch: row) {
            std::cout << patch;
        }
    }
}

#endif /* MAIN_H_ */
