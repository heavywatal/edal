/*! @file simulation.hpp
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
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <iosfwd>
#include <string>
#include <vector>
#include <memory>
#include <random>

namespace boost {namespace program_options {
  class options_description;
  class variables_map;
}}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

class Patch;

/*! @brief Represents single run
*/
class Simulation {
  private:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    /*! @addtogroup biol_param
        @{*/

    //! \f$K_0\f$, the number of individuals in an initial patch
    size_t INITIAL_PATCH_SIZE = 100;

    //! The number of patches on Y axis
    size_t NUM_ROWS = 4;

    //! The number of patches on X axis
    size_t NUM_COLS = 4;

    //! Dimension number of resources/traits
    size_t DIMENSIONS = 2;

    //! The overall number of generations to observe
    size_t ENTIRE_PERIOD = 1000;

    //! Interval between snapshots
    size_t OBSERVATION_CYCLE = 100;

    /** @} endgroup */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! Print extra information
    bool VERBOSE = false;

    //! Seed for random number generator
    unsigned int SEED = std::random_device{}();

    //! Control execution mode
    int MODE = 0;

    //! Set same values of \$h\$ and \$s\$ for two axes
    bool SYMMETRIC = false;

    //! Group name of this run such as altered parameter
    std::string LABEL;

  public:
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! Parse command arguments
    Simulation(int argc, char* argv[]);
    //! destructor in cpp for incomplete type
    ~Simulation();

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

    /** @} life_cycle */
    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

    //! output genotype frequencies for each patch
    void write_snapshot(const size_t time, std::ostream& ost) const;

    boost::program_options::options_description opt_description();

    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data members
    //! optional variables
    std::unique_ptr<boost::program_options::variables_map> vars_;

    //! Two-dimensional matrix of Patch
    std::vector<std::vector<Patch> > population;

};

} // namespace edal

#endif /* SIMULATION_HPP_ */
