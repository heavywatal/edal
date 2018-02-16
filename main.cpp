/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/simulation.hpp"

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    edal::Simulation simulation(argc, argv);
    simulation.run();
    return 0;
}
