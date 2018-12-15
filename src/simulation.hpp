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

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

class Patch;

/*! @brief Represents single run
*/
class Simulation {
  public:
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

    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // data members
    //! Two-dimensional matrix of Patch
    std::vector<std::vector<Patch> > population;

    //! Print verbose output
    bool verbose_ = false;
};

} // namespace edal

#endif /* SIMULATION_HPP_ */
