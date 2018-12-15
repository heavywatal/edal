/*! @file simulation.cpp
    @brief Implementation of Simulation class
*/
#include "simulation.hpp"
#include "patch.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/getopt.hpp>
#include <wtl/chrono.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/zlib.hpp>
#include <sfmt.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

namespace fs = wtl::filesystem;
namespace po = boost::program_options;

using URBG = wtl::sfmt19937_64;

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
boost::program_options::options_description Simulation::opt_description() {HERE;
    namespace po = boost::program_options;
    po::options_description description("Simulation");
    description.add_options()
        ("help,h", po::bool_switch(), "produce help")
        ("verbose,v", po::bool_switch(&VERBOSE), "verbose output")
        ("mode", po::value(&MODE)->default_value(MODE))
        ("symmetric", po::bool_switch(&SYMMETRIC))
        ("label", po::value(&LABEL)->default_value("default"))
        ("top_dir", po::value<std::string>()->default_value(wtl::strftime("edal%Y%m%d")))
        ("patch_size,k", po::value(&INITIAL_PATCH_SIZE)->default_value(INITIAL_PATCH_SIZE))
        ("row", po::value(&NUM_ROWS)->default_value(NUM_ROWS))
        ("col", po::value(&NUM_COLS)->default_value(NUM_COLS))
        ("dimensions,D", po::value(&DIMENSIONS)->default_value(DIMENSIONS))
        ("time,T", po::value(&ENTIRE_PERIOD)->default_value(ENTIRE_PERIOD))
        ("interval,I", po::value(&OBSERVATION_CYCLE)->default_value(OBSERVATION_CYCLE))
        ("seed", po::value(&SEED)->default_value(SEED))
    ;
    return description;
}

Simulation::~Simulation() = default;

Simulation::Simulation(int argc, char* argv[])
: vars_(std::make_unique<po::variables_map>()) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    std::vector<std::string> arguments(argv, argv + argc);
    wtl::join(arguments, std::cout, " ") << std::endl;
    std::cout << wtl::iso8601datetime() << std::endl;

    namespace po = boost::program_options;
    po::options_description description;
    description.add(opt_description());
    description.add(Individual::opt_description());
    auto& vm = *vars_;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm["help"].as<bool>()) {
        description.print(std::cout);
        exit(0);
    }
    wtl::sfmt64().seed(SEED);
    if (DIMENSIONS == 1) {
        std::ostringstream ost;
        ost << "diameter_pref = 1e6\n"
            << "limb_select = 1e6\n"
            << "mutation_mask = 10\n";
        std::istringstream ist(ost.str());
        po::store(po::parse_config_file(ist, description, false), vm);
        vm.notify();
    }
    if (SYMMETRIC) {
        std::ostringstream ost;
        ost << "diameter_pref = " << vm["height_pref"].as<double>() << "\n"
            << "limb_select = " << vm["toepad_select"].as<double>() << "\n";
        std::istringstream ist(ost.str());
        po::store(po::parse_config_file(ist, description, false), vm);
        vm.notify();
    }
    const std::string CONFIG_STRING = wtl::flags_into_string(vm);
    if (VERBOSE) {
        std::cout << CONFIG_STRING << std::endl;
    }
    if (ENTIRE_PERIOD % OBSERVATION_CYCLE > 0) {
        std::cerr << "T=" << ENTIRE_PERIOD
                  << " is not a multiple of I="
                  << OBSERVATION_CYCLE << std::endl;
        exit(1);
    }
    fs::path OUT_DIR(vm["top_dir"].as<std::string>());
    fs::create_directory(OUT_DIR);
    const std::string now(wtl::strftime("%Y%m%d_%H%M%S"));
    OUT_DIR /= (LABEL + "_" + now + "_" + std::to_string(::getpid()));
    fs::create_directory(OUT_DIR);
    fs::current_path(OUT_DIR.string());
    wtl::make_ofs("program_options.conf") << CONFIG_STRING;
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
        wtl::zlib::ofstream{"possible_geographic.csv.gz"} << Individual::possible_geographic();
        wtl::zlib::ofstream{"possible_phenotypes.csv.gz"} << Individual::possible_phenotypes();
        break;
      default:
        exit(1);
    }
    std::cout << wtl::iso8601datetime() << std::endl;
}

void Simulation::evolve() {HERE;
    assert(ENTIRE_PERIOD % OBSERVATION_CYCLE == 0);
    population.reserve(NUM_ROWS);
    for (size_t row=0; row<NUM_ROWS; ++row) {
        std::vector<Patch> pop_row;
        pop_row.reserve(NUM_COLS);
        for (size_t col=0; col<NUM_COLS; ++col) {
            pop_row.emplace_back(wtl::sfmt64()());
        }
        population.emplace_back(std::move(pop_row));
    }
    if (DIMENSIONS == 1) {
        population[0][0].assign(INITIAL_PATCH_SIZE, Individual{{15, 0, 15, 0}});
    } else {
        population[0][0].assign(INITIAL_PATCH_SIZE, Individual{});
    }
    std::ostringstream ost;
    for (size_t t=0; t<=ENTIRE_PERIOD; ++t) {
        if (VERBOSE) {
            std::cout << "\nT = " << t << "\n"
                << wtl::str_matrix(population, " ", wtl::make_oss(),
                        [](const Patch& p) {return p.size();})
                << std::flush;
        }
        if (t % OBSERVATION_CYCLE == 0) {
            write_snapshot(t, ost);
        }
        if (t < ENTIRE_PERIOD) {
            life_cycle();
        }
    }
    wtl::zlib::ofstream{"evolution.csv.gz"} << ost.str();
}

void Simulation::life_cycle() {
    static wtl::ThreadPool pool(2u);
    const size_t num_patches = NUM_ROWS * NUM_COLS;
    std::vector<std::future<
      std::pair<std::vector<Individual>, std::vector<std::pair<unsigned int, unsigned int>>>
    >> futures;
    futures.reserve(num_patches);
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            futures.emplace_back(pool.submit([this](size_t r, size_t c){
                auto& patch = population[r][c];
                auto children = patch.mate_and_reproduce();
                return std::pair<std::vector<Individual>, std::vector<std::pair<unsigned int, unsigned int>>>{
                  children, patch.make_destinations(children.size(), r, c, NUM_ROWS, NUM_COLS)
                };
            }, row, col));
        }
    }
    pool.wait();
    for (auto& ftr: futures) {
        auto children_destination = ftr.get();
        auto& children_i = children_destination.first;
        const auto& destinations_i = children_destination.second;
        const size_t num_children = children_i.size();
        for (size_t j=0; j<num_children; ++j) {
            const auto& dst = destinations_i[j];
            population[dst.first][dst.second].emplace_back(std::move(children_i[j]));
        }
    }
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            auto& patch = population[row][col];
            pool.submit([&patch](){patch.viability_selection();});
        }
    }
    pool.wait();
}

void Simulation::write_snapshot(const size_t time, std::ostream& ost) const {
    const char* sep = ",";
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
    DCERR("N = " << popsize << std::endl);
}

} // namespace edal
