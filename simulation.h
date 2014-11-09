// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.h
    @brief Interface of Simulation class
    @defgroup biol_param Biological parameters
    @defgroup biol_proc Biological processes
    @{
        @defgroup habitat_pareference Habitat preference
        @defgroup natural_selection Natural selection
        @defgroup mating Mating
        @defgroup life_cycle Life cycle
    @}
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

/*! @brief Represents single run
*/
class Simulation {
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /*! @addtogroup biol_param
        @{*/

    //! The number of patches on Y axis
    size_t NUM_ROWS = 8;

    //! The number of patches on X axis
    size_t NUM_COLS = 8;

    //! The number of individuals in an original patch
    size_t INITIAL_PATCH_SIZE = 40;

    //! The overall number of generations to observe
    size_t ENTIRE_PERIOD = 1000;

    //! Interval between snapshots
    size_t OBSERVATION_CYCLE = 100;

    /** @} endgroup */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! The number of CPU cores to use
    size_t PPN;

    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Control execution mode
    int MODE = 0;

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

  public:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! Parse command arguments
    Simulation(int argc, char* argv[]);

    //! Top level function that should be called from main()
    void run();

  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /** @addtogroup life_cycle
        @{*/

    //! Call life_cycle() repeatedly
    void evolve();

    /*! @brief One step forward

        > There are two options.
        > (1) viability selection → dispersal → mating and offspring production,
        > and (2) dispersal → viability selection → mating and offspring production.
        > [Q: which one is more appropriate for lizards?
        > Jonathan: I think dispersal before selection is more appropriate,
        > although we actually know little about dispersal;
        > still, the best bet is that it occurs when they are young.]
        > So, we consider the second option only.
    */
    void life_cycle();

    /*! @brief Choose migration desination randomly
        @param row_orig Original Y position
        @param col_orig Original X position
        @return new position (row, col)

        > With probability \f$ m > 0 \f$, each offspring becomes a "migrant."
        > Each migrant goes to one of the 8 neighboring patches.
        > For patches at the boundary,
        > the probability \f$ m \f$ is reduced according to the number of neighbors they have.
    */
    std::pair<size_t, size_t> choose_destination(const size_t row_orig, const size_t col_orig);

    /** @} life_cycle */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    void write_snapshot(const size_t time, std::ostream& ost) const;

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
    // data members
    //! Two-dimensional matrix of Patch
    std::vector<std::vector<Patch> > population;

};

#endif /* SIMULATION_H_ */
