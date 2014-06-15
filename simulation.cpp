// -*- mode: c++; coding: utf-8 -*-
/*! @file simulation.cpp
    @brief Inplementation of Simulation class
*/
#include "simulation.h"

#include "cxxwtils/iostr.hpp"
#include "cxxwtils/getopt.hpp"
#include "cxxwtils/prandom.hpp"
#include "cxxwtils/os.hpp"
#include "cxxwtils/gz.hpp"
#include "cxxwtils/multiprocessing.hpp"

#include "individual.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
// functions

boost::program_options::options_description& Simulation::opt_description() {HERE;
    namespace po = boost::program_options;
    static po::options_description description("Simulation");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
        ("verbose,v", po::value<bool>(&VERBOSE)
            ->default_value(VERBOSE)->implicit_value(true), "verbose output")
        ("test", po::value<bool>()->default_value(false)->implicit_value(true))
        ("ppn", po::value<size_t>(&PPN)->default_value(PPN))
        ("label", po::value<std::string>(&LABEL)->default_value("default"))
        ("top_dir", po::value<std::string>()->default_value(OUT_DIR.string()))
        ("row", po::value<size_t>(&NUM_ROWS)->default_value(NUM_ROWS))
        ("col", po::value<size_t>(&NUM_COLS)->default_value(NUM_COLS))
        ("time,T", po::value<size_t>(&OBSERVATION_PERIOD)->default_value(OBSERVATION_PERIOD))
        ("migration_rate,m", po::value<double>(&MIGRATION_RATE)->default_value(MIGRATION_RATE))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED))
    ;
    return description;
}

//! Unit test for each class
inline void test() {HERE;
    individual_unit_test();
    patch_unit_test();
}

Simulation::Simulation(int argc, char* argv[]) {HERE;
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
    OUT_DIR = fs::path(vm["top_dir"].as<std::string>());
    const std::string CONFIG_STRING = flags_into_string(description, vm);
    prandom().seed(SEED); // TODO: want to read seed?
    if (VERBOSE) {
        std::cout << CONFIG_STRING << std::endl;
    }
    if (vm["test"].as<bool>()) {
        test();
        exit(0);
    }
    const std::string now(wtl::strftime("%Y%m%d_%H%M%S"));
    std::ostringstream pid_at_host;
    pid_at_host << ::getpid() << "@" << wtl::gethostname();
    WORK_DIR = TMP_DIR / (now + "_" + LABEL + "_" + pid_at_host.str());
    derr("mkdir && cd to " << WORK_DIR << std::endl);
    fs::create_directory(WORK_DIR);
    wtl::cd(WORK_DIR.string());
    fs::create_directory(OUT_DIR);
    OUT_DIR /= (LABEL + "_" + now + "_" + pid_at_host.str());
    wtl::Fout{"program_options.conf"} << CONFIG_STRING;
}

void Simulation::run() {HERE;
    population.assign(NUM_ROWS, std::vector<Patch>(NUM_COLS));
    population[0][0] = Patch(INITIAL_PATCH_SIZE);
    for (size_t i=0; i<OBSERVATION_PERIOD; ++i) {
        if (VERBOSE) {
            std::cout << "T = " << i << "\n"
                << str_population([](const Patch& p) {return p.size();})
                << std::endl;
        }
        life_cycle();
    }
    if (VERBOSE) {
        std::cout << "T = " << OBSERVATION_PERIOD << "\n"
            << str_population([](const Patch& p) {return p.size();})
            << std::endl;
    }
    wtl::gzip{wtl::Fout{"population.csv.gz"}}
        << Individual::header()
        << str_population([](const Patch& p) {return p;}, "", "");
    derr("mv results to " << OUT_DIR << std::endl);
    fs::rename(WORK_DIR, OUT_DIR);
    std::cout << wtl::iso8601datetime() << std::endl;
}

void Simulation::life_cycle() {
    std::vector<std::vector<Patch> > next_generation(population);
    wtl::Semaphore sem(PPN);
    std::mutex mtx;
    auto patch_task = [&](const size_t row, const size_t col) {
        auto offsprings = population[row][col].mate_and_reproduce();
        for (const auto& child: offsprings) {
                if (prandom().bernoulli(MIGRATION_RATE)) {
                    auto new_coords = choose_destination(row, col);
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

std::pair<size_t, size_t> Simulation::choose_destination(const size_t row_orig, const size_t col_orig) {
    size_t row = row_orig;
    size_t col = col_orig;
    switch (prandom().randrange(8)) {
      case 0:        ++col; break;
      case 1: ++row; ++col; break;
      case 2: ++row;        break;
      case 3: ++row; --col; break;
      case 4:        --col; break;
      case 5: --row; --col; break;
      case 6: --row;        break;
      case 7: --row; ++col; break;
    }
    if ((row > 7) | (col > 7)) {  // i.e. 8 or -1
        row = row_orig; col = col_orig;
    }
    return {row, col};
}