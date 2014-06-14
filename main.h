// -*- mode: c++; coding: utf-8 -*-
/** @file main.h
    @brief Simulation functions
*/
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
#include "cxxwtils/os.hpp"
#include "cxxwtils/gz.hpp"
#include "cxxwtils/multiprocessing.hpp"

#include "individual.h"
#include "patch.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// parameters

//! The number of patches on Y axis
size_t NUM_ROWS = 8;

//! The number of patches on X axis
size_t NUM_COLS = 8;

//! The number of individuals in an original patch
size_t INITIAL_PATCH_SIZE = 40;

//! The number of generations to observe
size_t OBSERVATION_PERIOD = 1000;

//! Migration rate per generation
double MIGRATION_RATE = 1e-1;

//! Seed for random number generator
unsigned int SEED = std::random_device{}();
std::string LABEL;
std::string CONFIG_STRING;

const char* HOME_ = std::getenv("HOME");
const fs::path HOME_DIR{HOME_};
const fs::path TMP_DIR{HOME_DIR / "tmp"};
fs::path WORK_DIR;
fs::path TOP_DIR;
fs::path LOAD_DIR;

//! The number of CPU cores to use
size_t PPN = 4;
bool VERBOSE = false;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// data members

//! Two-dimensional matrix of Patch
std::vector<std::vector<Patch> > population;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions

inline boost::program_options::options_description opt_description() {HERE;
    namespace po = boost::program_options;
    po::options_description description("Operation");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
        ("verbose,v", po::value<bool>(&VERBOSE)
            ->default_value(VERBOSE)->implicit_value(true), "verbose output")
        ("test", po::value<bool>()->default_value(false)->implicit_value(true))
        ("ppn", po::value<size_t>(&PPN)->default_value(PPN))
        ("label", po::value<std::string>(&LABEL)->default_value("default"))
        ("top_dir", po::value<std::string>()
            ->default_value((TMP_DIR / wtl::strftime("out%Y%m%d")).string()))
        ("row", po::value<size_t>(&NUM_ROWS)->default_value(NUM_ROWS))
        ("col", po::value<size_t>(&NUM_COLS)->default_value(NUM_COLS))
        ("time,T", po::value<size_t>(&OBSERVATION_PERIOD)->default_value(OBSERVATION_PERIOD))
        ("migration_rate,m", po::value<double>(&MIGRATION_RATE)->default_value(MIGRATION_RATE))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED))
    ;
    return description;
}

inline void test() {HERE;
    individual_unit_test();
    patch_unit_test();
}

inline void check_flags(int argc, char* argv[]) {HERE;
    std::vector<std::string> arguments(argv, argv + argc);
    std::cout << wtl::str_join(arguments, " ") << std::endl;
    std::cout << wtl::iso8601datetime() << std::endl;

    namespace po = boost::program_options;
    po::options_description description;
    description.add(opt_description());
    description.add(Individual::opt_description());
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm["help"].as<bool>()) {
        description.print(std::cout);
        exit(0);
    }
    TOP_DIR = fs::path(vm["top_dir"].as<std::string>());
    CONFIG_STRING = flags_into_string(description, vm);
    prandom().seed(SEED); // TODO: want to read seed?
    if (VERBOSE) {
        std::cout << CONFIG_STRING << std::endl;
    }
    if (vm["test"].as<bool>()) {
        test();
        exit(0);
    }
}

inline std::pair<size_t, size_t> migrate(const size_t row_orig, const size_t col_orig) {
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
    if ((row > 7) | (col > 7)) {row = row_orig; col = col_orig;}
    return {row, col};
}

inline void life_cycle() {
    std::vector<std::vector<Patch> > next_generation(population);
    wtl::Semaphore sem(PPN);
    std::mutex mtx;
    auto patch_task = [&](const size_t row, const size_t col) {
        auto offsprings = population[row][col].mate_and_reproduce();
        for (const auto& child: offsprings) {
                if (prandom().bernoulli(MIGRATION_RATE)) {
                    auto new_coords = migrate(row, col);
                    std::lock_guard<std::mutex> lck(mtx);
                    next_generation[new_coords.first][new_coords.second].append(child);
                } else {
                    std::lock_guard<std::mutex> lck(mtx);
                    next_generation[row][col].append(child);
                }
        }
        sem.unlock();
    };
    std::vector<std::thread> threads;
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            sem.lock();
            threads.emplace_back(patch_task, row, col);
        }
    }
    for (auto& th: threads) {th.join();}
    threads.clear();
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            sem.lock();
            threads.emplace_back([row, col, &sem, &next_generation] {
                next_generation[row][col].viability_selection();
                sem.unlock();
            });
        }
    }
    for (auto& th: threads) {th.join();}
    population.swap(next_generation);
}

template <class Vec2D, class Func> inline
std::string str_population(const Vec2D& m, Func func,
                      const std::string sep_col=" ", const std::string sep_row="\n") {
    std::ostringstream ost;
    for (const auto& row: m) {
        for (const auto& cell: row) {
            ost << func(cell) << sep_col;
        }
        ost << sep_row;
    }
    return ost.str();
}

inline void run() {HERE;
    const std::string now(wtl::strftime("%Y%m%d_%H%M%S"));
    std::ostringstream pid_at_host;
    pid_at_host << ::getpid() << "@" << wtl::gethostname();
    WORK_DIR = TMP_DIR / (now + "_" + LABEL + "_" + pid_at_host.str());
    derr(WORK_DIR << std::endl);
    const fs::path job_dir = TOP_DIR / (LABEL + "_" + now + "_" + pid_at_host.str());
    fs::create_directory(WORK_DIR);
    wtl::cd(WORK_DIR.string());
    wtl::Fout{"program_options.conf"} << CONFIG_STRING;

    population.assign(NUM_ROWS, std::vector<Patch>(NUM_COLS));
    population[0][0] = Patch(INITIAL_PATCH_SIZE);
    for (size_t i=0; i<OBSERVATION_PERIOD; ++i) {
        if (VERBOSE) {
            std::cout << "T = " << i << "\n"
                << str_population(population, [](const Patch& p) {return p.size();})
                << std::endl;
        }
        life_cycle();
    }
    if (VERBOSE) {
        std::cout << "T = " << OBSERVATION_PERIOD << "\n"
            << str_population(population, [](const Patch& p) {return p.size();})
            << std::endl;
    }
    wtl::gzip{wtl::Fout{"population.csv.gz"}}
        << Individual::header()
        << str_population(population, [](const Patch& p) {return p;}, "", "");

    derr(job_dir << std::endl);
    fs::create_directory(TOP_DIR);
    fs::rename(WORK_DIR, job_dir);
    std::cout << wtl::iso8601datetime() << std::endl;
}

#endif /* MAIN_H_ */
