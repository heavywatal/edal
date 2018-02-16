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
#include <wtl/zfstream.hpp>
#include <sfmt.hpp>
#include <boost/asio.hpp>
#include <boost/filesystem.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

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
        ("test", po::value<int>()->default_value(0)->implicit_value(1))
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

//! Unit test for each class
inline void test() {HERE;
    Individual::unit_test();
    Patch::unit_test();
}

Simulation::~Simulation() {HERE;}

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
    wtl::sfmt().seed(SEED);
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
    switch (vm["test"].as<int>()) {
      case 0:
        break;
      case 1:
        test();
        exit(0);
      default:
        exit(1);
    }
    fs::path OUT_DIR(vm["top_dir"].as<std::string>());
    fs::create_directory(OUT_DIR);
    const std::string now(wtl::strftime("%Y%m%d_%H%M%S"));
    std::ostringstream pid_at_host;
    pid_at_host << ::getpid() << "@" << boost::asio::ip::host_name();
    OUT_DIR /= (LABEL + "_" + now + "_" + pid_at_host.str());
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
        wtl::ozfstream{"possible_geographic.csv.gz"} << Individual::possible_geographic();
        wtl::ozfstream{"possible_phenotypes.csv.gz"} << Individual::possible_phenotypes();
        break;
      default:
        exit(1);
    }
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
    wtl::ozfstream{"evolution.csv.gz"} << ost.str();
}

/*! @brief Change row/col with probability \f$m\f$ = Individual::MIGRATION_RATE_

    > With probability \f$ m > 0 \f$, each offspring becomes a "migrant."
    > Each migrant goes to one of the 8 neighboring patches.
    > For patches at the boundary,
    > the probability \f$ m \f$ is reduced according to the number of neighbors they have.
*/
inline std::pair<size_t, size_t> choose_patch(size_t row, size_t col) {
    if (!std::bernoulli_distribution(Individual::MIGRATION_RATE())(wtl::sfmt()))
        {return {row, col};}
    switch (std::uniform_int_distribution<size_t>(0, 7)(wtl::sfmt())) {
      case 0:        ++col; break;
      case 1: ++row; ++col; break;
      case 2: ++row;        break;
      case 3: ++row; --col; break;
      case 4:        --col; break;
      case 5: --row; --col; break;
      case 6: --row;        break;
      case 7: --row; ++col; break;
    }
    return {row, col};
}

void Simulation::life_cycle() {
    std::vector<std::vector<Patch> > parents(population);
    auto reproduction = [&](const size_t row, const size_t col) {
        auto offsprings = parents[row][col].mate_and_reproduce(wtl::sfmt());
        std::vector<std::pair<size_t, size_t> > destinations;
        destinations.reserve(offsprings.size());
        for (size_t i=0; i<offsprings.size(); ++i) {
            const auto dst = choose_patch(row, col);
            if ((dst.first >= NUM_ROWS) | (dst.second >= NUM_COLS)) {
                destinations.emplace_back(row, col);
            } else {
                destinations.push_back(std::move(dst));
            }
        }
        for (size_t i=0; i<offsprings.size(); ++i) {
            const auto& dst = destinations[i];
            population[dst.first][dst.second].append(std::move(offsprings[i]));
        }
    };
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            reproduction(row, col);
        }
    }
    for (size_t row=0; row<NUM_ROWS; ++row) {
        for (size_t col=0; col<NUM_COLS; ++col) {
            population[row][col].viability_selection(wtl::sfmt());
        }
    }
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
