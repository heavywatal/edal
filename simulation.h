// -*- mode: c++; coding: utf-8 -*-
/** @file simulation.h
    @brief Interface of Simulation class
*/
#pragma once
#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <random>

#include "cxxwtils/iostr.hpp"

#include "patch.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

namespace boost {
    namespace program_options {
        class options_description;
    }
}

/** @brief Represents single run
*/
class Simulation {
  public:
    //! Construct with command arguments
    Simulation(int argc, char* argv[]);

    //! Main body of the program
    void run();

  private:
    //! One step forward
    void life_cycle();

    /** @brief Calculate new position
        @param row_orig Original Y position
        @param col_orig Original X position
        @return new position (row, col)
    */
    std::pair<size_t, size_t> migrate(const size_t row_orig, const size_t col_orig);

    //! Program options for this class
    boost::program_options::options_description& opt_description();

    //! apply function for each patch and concatenate output into a string
    template <class Func> inline
    std::string str_population(Func func,
                               const std::string sep_col=" ",
                               const std::string sep_row="\n") {
        std::ostringstream ost;
        for (const auto& row: population) {
            for (const auto& cell: row) {
                ost << func(cell) << sep_col;
            }
            ost << sep_row;
        }
        return ost.str();
    }

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

    //! The number of CPU cores to use
    size_t PPN = 4;
    
    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();
    
    //! Group name of this run such as altered parameter
    std::string LABEL;

    //! Home directory (c_str)
    const char* HOME_CHAR = std::getenv("HOME");

    //! Home directory
    const fs::path HOME_DIR = HOME_CHAR;

    //! Temporary directory
    const fs::path TMP_DIR = HOME_DIR / "tmp";

    //! Working directory to write out temporal results
    fs::path WORK_DIR;

    //! Target directory to which the contents in WORK_DIR are moved
    fs::path OUT_DIR = TMP_DIR / wtl::strftime("out%Y%m%d");

    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data members
    //! Two-dimensional matrix of Patch
    std::vector<std::vector<Patch> > population;

};

#endif /* SIMULATION_H_ */
