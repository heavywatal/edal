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

//! Symbols for the program options can be different from those in equations
/*! @ingroup biol_param
    @return Program options description

    Command line option | Symbol    | Variable
    ------------------- | --------- | ------------------------------------------
    `-k,--patch_size`   | \f$K_0\f$ | Simulation::INITIAL_PATCH_SIZE
    `--row,--col`       | -         | Simulation::NUM_COLS, Simulation::NUM_ROWS
    `-T,--time`         | -         | Simulation::ENTIRE_PERIOD
    `-I,--interval`     | -         | Simulation::OBSERVATION_CYCLE
*/
boost::program_options::options_description& Simulation::opt_description() {HERE;
    namespace po = boost::program_options;
    static po::options_description description("Simulation");
    description.add_options()
        ("help,h", po::value<bool>()->default_value(false)->implicit_value(true), "produce help")
        ("verbose,v", po::value<bool>(&VERBOSE)
            ->default_value(VERBOSE)->implicit_value(true), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1))
        ("mode", po::value<int>(&MODE)->default_value(MODE))
        ("ppn", po::value<size_t>(&PPN)->default_value(wtl::num_threads()))
        ("label", po::value<std::string>(&LABEL)->default_value("default"))
        ("top_dir", po::value<std::string>()->default_value(OUT_DIR.string()))
        ("patch_size,k", po::value<size_t>(&INITIAL_PATCH_SIZE)->default_value(INITIAL_PATCH_SIZE))
        ("row", po::value<size_t>(&NUM_ROWS)->default_value(NUM_ROWS))
        ("col", po::value<size_t>(&NUM_COLS)->default_value(NUM_COLS))
        ("dimensions,D", po::value<size_t>(&DIMENSIONS)->default_value(DIMENSIONS))
        ("time,T", po::value<size_t>(&ENTIRE_PERIOD)->default_value(ENTIRE_PERIOD))
        ("interval,I", po::value<size_t>(&OBSERVATION_CYCLE)->default_value(OBSERVATION_CYCLE))
        ("seed", po::value<unsigned int>(&SEED)->default_value(SEED))
    ;
    return description;
}

//! Unit test for each class
inline void test() {HERE;
    Individual::unit_test();
    Patch::unit_test();
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
    prandom().seed(SEED); // TODO: want to read seed?
    if (DIMENSIONS == 1) {
        std::ostringstream ost;
        ost << "diameter_pref = 1e6\n"
            << "limb_select = 1e6\n"
            << "mutation_mask = 10\n";
        std::istringstream ist(ost.str());
        po::store(po::parse_config_file(ist, description, false), vm);
        vm.notify();
    }
    const std::string CONFIG_STRING = wtl::flags_into_string(description, vm);
    if (VERBOSE) {
        std::cout << CONFIG_STRING << std::endl;
    }
    if (ENTIRE_PERIOD % OBSERVATION_CYCLE > 0) {
        std::cerr << wtl::strprintf(
            "T=%d is not a multiple of I=%d",
            ENTIRE_PERIOD, OBSERVATION_CYCLE) << std::endl;
        exit(1);
    }
    switch (vm["test"].as<int>()) {
      case 0:
        break;
      case 1:
        test();
        exit(0);
      default:
        exit(1);
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
    switch (MODE) {
      case 0:
        evolve();
        break;
      case 1:
        Individual::write_resource_abundance();
        break;
      case 2:
        wtl::gzip{wtl::Fout{"possible_geographic.csv.gz"}} << Individual::possible_geographic();
        wtl::gzip{wtl::Fout{"possible_phenotypes.csv.gz"}} << Individual::possible_phenotypes();
        break;
      default:
        exit(1);
    }
    derr("mv results to " << OUT_DIR << std::endl);
    fs::rename(WORK_DIR, OUT_DIR);
    std::cout << wtl::iso8601datetime() << std::endl;
}

void Simulation::evolve() {HERE;
    assert(ENTIRE_PERIOD % OBSERVATION_CYCLE == 0);
    population.assign(NUM_ROWS, std::vector<Patch>(NUM_COLS));
    if (DIMENSIONS == 1) {
        population[0][0] = Patch(INITIAL_PATCH_SIZE, Individual{{15, 0, 15, 0}});
    } else {
        population[0][0] = Patch(INITIAL_PATCH_SIZE);
    }
    std::ostringstream ost;
    for (size_t t=0; t<=ENTIRE_PERIOD; ++t) {
        if (VERBOSE) {
            std::cout << "\nT = " << t << "\n"
                << str_population([](const Patch& p) {return p.size();})
                << std::flush;
        }
        if (t % OBSERVATION_CYCLE == 0) {
            write_snapshot(t, ost);
        }
        if (t < ENTIRE_PERIOD) {
            life_cycle();
        }
    }
    wtl::gzip{wtl::Fout{"evolution.csv.gz"}} << ost.str();
}

void Simulation::life_cycle() {
    std::vector<std::vector<Patch> > parents(population);
    wtl::Semaphore sem(PPN);
    std::mutex mtx;
    auto patch_task = [&](const size_t row, const size_t col) {
        auto offsprings = parents[row][col].mate_and_reproduce();
        std::lock_guard<std::mutex> lck(mtx);
        for (auto& child: offsprings) {
            size_t dest_row = row;
            size_t dest_col = col;
            choose_patch(&dest_row, &dest_col);
            population[dest_row][dest_col].append(std::move(child));
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
            threads.emplace_back([row, col, &sem, this] {
                population[row][col].viability_selection();
                sem.unlock();
            });
        }
    }
    for (auto& th: threads) {th.join();}
}

void Simulation::choose_patch(size_t* row, size_t* col) const {
    if (!prandom().bernoulli(Individual::MIGRATION_RATE())) {return;}
    size_t r = *row;
    size_t c = *col;
    switch (prandom().randrange(8)) {
      case 0:      ++c; break;
      case 1: ++r; ++c; break;
      case 2: ++r;      break;
      case 3: ++r; --c; break;
      case 4:      --c; break;
      case 5: --r; --c; break;
      case 6: --r;      break;
      case 7: --r; ++c; break;
    }
    // size_t(-1) also makes true
    if ((r >= NUM_ROWS) | (c >= NUM_COLS)) {return;}
    *row = r; *col = c;
}

void Simulation::write_snapshot(const size_t time, std::ostream& ost) const {
    const std::string sep{","};
    if (time == 0) {
        ost << "time" << sep << "row" << sep << "col" << sep
            << "n" << sep << Individual::header() << "\n";
    }
    std::size_t popsize = 0;
    for (size_t row=0; row<population.size(); ++row) {
        for (size_t col=0; col<population[row].size(); ++col) {
            for (const auto& item: population[row][col].summarize()) {
                ost << time << sep << row << sep << col << sep
                    << item.second << sep << item.first << "\n";
                popsize += item.second;
            }
        }
    }
    derr("N = " << popsize << std::endl);
}

