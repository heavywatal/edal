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
#include "cxxwtils/os.hpp"
#include "cxxwtils/multiprocessing.hpp"

#include "individual.h"
#include "patch.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// parameters

size_t num_rows = 8;
size_t num_cols = 8;
size_t initial_patch_size = 20;
size_t observation_period = 1000;
double migration_rate = 1e-1;

unsigned int seed_ = std::random_device{}();
std::string label_;
std::string config_string_;

const char* HOME_ = std::getenv("HOME");
const fs::path HOME_DIR{HOME_};
const fs::path TMP_DIR{HOME_DIR / "tmp"};
fs::path WORK_DIR;
fs::path TOP_DIR;
fs::path LOAD_DIR;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// data members

std::vector<std::vector<Patch> > population;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions


inline void test() {HERE;
    individual_unit_test();
    patch_unit_test();
}

inline void check_flags(int argc, char* argv[]) {HERE;
    namespace po = boost::program_options;
    po::options_description description("Operation");
    description.add_options()
            ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
            ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "verbose output")
            ("test", po::value<bool>()->default_value(false)->implicit_value(true))
            ("label", po::value<std::string>(&label_)->default_value("default"))
            ("top_dir", po::value<std::string>()
                ->default_value((TMP_DIR / wtl::strftime("out%Y%m%d")).string()))
            ("row", po::value<size_t>(&num_rows)->default_value(num_rows))
            ("col", po::value<size_t>(&num_cols)->default_value(num_cols))
            ("time,T", po::value<size_t>(&observation_period)->default_value(observation_period))
            ("migration_rate,m", po::value<double>(&migration_rate)->default_value(migration_rate))
            ("seed", po::value<unsigned int>(&seed_)->default_value(seed_))
            ;
    description.add(Individual::opt_description());
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    if (vm["help"].as<bool>()) {
        description.print(std::cout);
        exit(0);
    }
    TOP_DIR = fs::path(vm["top_dir"].as<std::string>());
    prandom().seed(seed_); // TODO: want to read seed?
    vm.notify();
    config_string_ = flags_into_string(description, vm);
    if (vm["verbose"].as<bool>()) {
        std::cout << description << std::endl;
        std::cout << config_string_ << std::endl;
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
    wtl::Semaphore sem(4);
    std::mutex mtx;
    auto patch_task = [&](const size_t row, const size_t col) {
        auto offsprings = population[row][col].mate_and_reproduce();
        for (const auto& child: offsprings) {
                if (prandom().bernoulli(migration_rate)) {
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
    for (size_t row=0; row<num_rows; ++row) {
        for (size_t col=0; col<num_cols; ++col) {
            sem.lock();
            threads.emplace_back(patch_task, row, col);
        }
    }
    for (auto& th: threads) {th.join();}
    threads.clear();
    for (size_t row=0; row<num_rows; ++row) {
        for (size_t col=0; col<num_cols; ++col) {
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
    WORK_DIR = TMP_DIR / (now + "_" + label_ + "_" + pid_at_host.str());
    derr(WORK_DIR << std::endl);
    const fs::path job_dir = TOP_DIR / (label_ + "_" + now + "_" + pid_at_host.str());
    fs::create_directory(WORK_DIR);
    wtl::cd(WORK_DIR.string());
    wtl::Fout{"program_options.conf"} << config_string_;

    population.assign(num_rows, std::vector<Patch>(num_cols));
    population[0][0] = Patch(initial_patch_size);
    for (size_t i=0; i<observation_period; ++i) {
        std::cout << "T = " << i << std::endl;
        std::cout << str_population(population, [](const Patch& p) {return p.size();});
        std::cout << std::endl;
        life_cycle();
    }
    std::cout << "T = " << observation_period << std::endl;
    std::cout << str_population(population, [](const Patch& p) {return p.size();});
    std::cout << std::endl;
    wtl::Fout{"population.csv"}
        << Individual::header()
        << str_population(population, [](const Patch& p) {return p;}, "", "");

    derr(job_dir << std::endl);
    fs::create_directory(TOP_DIR);
    fs::rename(WORK_DIR, job_dir);
}

#endif /* MAIN_H_ */
