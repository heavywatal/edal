# Evolutionary Diversification in Anolis Lizards

Simulation program is written in C++.
The results are analyzed in R.

- [Doxygen documentation](https://heavywatal.github.io/edal/)
- [Model description](model.html)
- [Project page on GitHub](https://github.com/heavywatal/edal)

## Requirements

- Unix-like environment (macOS, Linux, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.12.0)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [zlib](https://zlib.net)

## Parameters and functions

- @ref biol_param
- @ref habitat_pareference
- @ref natural_selection
- @ref mating
- @ref life_cycle

## Results

- [Possible Ke values with sojourn time normalization (anolis_v3)](http://meme.biology.tohoku.ac.jp/edal/results/ke_v3.html)
- [Simulation results with sojourn time normalization (anolis_v3)](http://meme.biology.tohoku.ac.jp/edal/results/sim_v3.html)
- [Possible Ke values with habitat preference normalization (anolis_v3a)](http://meme.biology.tohoku.ac.jp/edal/results/ke_v3a.html)
- [Evolutionary branching on simple condition](http://meme.biology.tohoku.ac.jp/edal/results/ad_20141112.html)
- [Evolutionary branching on simple condition (Gaussian form)](http://meme.biology.tohoku.ac.jp/edal/results/ad_20141122.html)
- [Isolation by distance in one-dimensional stepping-stone](http://meme.biology.tohoku.ac.jp/edal/results/ibd_20141203.html)
- [Evolutionary branching on simple condition (non-random mating by Débarre 2012)](http://meme.biology.tohoku.ac.jp/edal/results/ad_20141206.html)
- [Evolutionary branching on simple condition (non-random mating by TPG2013)](http://meme.biology.tohoku.ac.jp/edal/results/ad_20141210.html)
- [Evolutionary branching on simple condition (non-random mating by TPH2011)](http://meme.biology.tohoku.ac.jp/edal/results/ad_20141122-d.html)
- [Mutation-drift balance](http://meme.biology.tohoku.ac.jp/edal/results/mutation_drift.html)
