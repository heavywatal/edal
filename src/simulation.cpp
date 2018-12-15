/*! @file simulation.cpp
    @brief Implementation of Simulation class
*/
#include "simulation.hpp"
#include "patch.hpp"
#include "individual.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/concurrent.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/zlib.hpp>
#include <clippson/clippson.hpp>
#include <sfmt.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

namespace fs = wtl::filesystem;

//! Global variables mapper of commane-line arguments
nlohmann::json VM;

//! Options description for general purpose
//! Symbols for the program options can be different from those in equations
/*! @ingroup biol_param
    @return Program options description

    Command line option | Symbol    |
    :------------------ | :-------- |
    `-k,--patch_size`   | \f$K_0\f$ |
    `--row,--col`       | -         |
    `-T,--time`         | -         |
    `-I,--interval`     | -         |
*/
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      wtl::option(vm, {"h", "help"}, false, "Print this help"),
      wtl::option(vm, {"version"}, false, "Print version"),
      wtl::option(vm, {"v", "verbose"}, false, "Verbose output")
    ).doc("General:");
}

inline clipp::group simulation_options(nlohmann::json* vm) {
    const std::string TOP_DIR = wtl::strftime("edal%Y%m%d");
    const int seed = static_cast<int>(std::random_device{}()); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"mode"}, 0, "Control execution mode"),
      wtl::option(vm, {"symmetric"}, false, "Set same values of $h$ and $s$ for two axes"),
      wtl::option(vm, {"label"}, "default", "Group name of this run such as altered parameter"),
      wtl::option(vm, {"top_dir"}, TOP_DIR),
      wtl::option(vm, {"k", "patch_size"}, 100u, "$K_0$, the number of individuals in an initial patch"),
      wtl::option(vm, {"nrow"}, 4u, "The number of patches on Y axis"),
      wtl::option(vm, {"ncol"}, 4u, "The number of patches on X axis"),
      wtl::option(vm, {"D", "dimensions"}, 2u, "Dimension number of resources/traits"),
      wtl::option(vm, {"T", "time"}, 1000u, "The overall number of generations to observe"),
      wtl::option(vm, {"I", "interval"}, 100u, "Interval between snapshots"),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Simulation:");
}

//! Symbols for the program options can be different from those in equations
/*! @ingroup biol_param
    @return Program options description

    Command line option      | Symbol         | Variable
    :----------------------- | :------------- | :-------
    `-a,--beta_param`        | \f$\alpha\f$   | IndividualParams::BETA_PARAM
    `-K,--carrying_capacity` | \f$K_{\max}\f$ | IndividualParams::CARRYING_CAPACITY
    `-b,--birth_rate`        | \f$b\f$        | IndividualParams::AVG_NUM_OFFSPINRGS
    `-p,--height_pref`       | \f$h_0\f$      | IndividualParams::HEIGHT_PREFERENCE
    `-P,--diameter_pref`     | \f$h_1\f$      | IndividualParams::DIAMETER_PREFERENCE
    `-s,--toepad_select`     | \f$s_0\f$      | IndividualParams::TOEPAD_SELECTION
    `-S,--limb_select`       | \f$s_1\f$      | IndividualParams::LIMB_SELECTION
    `-c,--pref_compe`        | \f$c_y\f$      | IndividualParams::PREF_COMPETITION
    `-C,--morph_compe`       | \f$c_x\f$      | IndividualParams::MORPH_COMPETITION
    `-f,--mating_sigma`      | \f$\sigma_a\f$ | IndividualParams::MATING_SIGMA
    `-u,--mu_locus`          | -              | IndividualParams::MU_LOCUS
    `-U,--mutation_mask`     | -              | IndividualParams::MUTATION_MASK
    `-m,--migration_rate`    | \f$m\f$        | IndividualParams::MIGRATION_RATE
*/
inline clipp::group individual_options(nlohmann::json* vm, IndividualParams* p) {
    return (
      wtl::option(vm, {"a", "beta_param"}, &p->BETA_PARAM),
      wtl::option(vm, {"K", "carrying_capacity"}, &p->CARRYING_CAPACITY),
      wtl::option(vm, {"b", "birth_rate"}, &p->AVG_NUM_OFFSPINRGS),
      wtl::option(vm, {"p", "height_pref"}, &p->HEIGHT_PREFERENCE),
      wtl::option(vm, {"P", "diameter_pref"}, &p->DIAMETER_PREFERENCE),
      wtl::option(vm, {"s", "toepad_select"}, &p->TOEPAD_SELECTION),
      wtl::option(vm, {"S", "limb_select"}, &p->LIMB_SELECTION),
      wtl::option(vm, {"c", "pref_compe"}, &p->PREF_COMPETITION),
      wtl::option(vm, {"C", "morph_compe"}, &p->MORPH_COMPETITION),
      wtl::option(vm, {"f", "mating_sigma"}, &p->MATING_SIGMA),
      wtl::option(vm, {"u", "mu_locus"}, &p->MU_LOCUS),
      wtl::option(vm, {"U", "mutation_mask"}, &p->MUTATION_MASK),
      wtl::option(vm, {"m", "migration_rate"}, &p->MIGRATION_RATE)
    ).doc("Individual:");
}

Simulation::~Simulation() = default;

Simulation::Simulation(int argc, char* argv[]) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);
    std::vector<std::string> arguments(argv + 1, argv + argc);
    wtl::join(arguments, std::cout, " ") << std::endl;

    nlohmann::json vm_local;
    IndividualParams individual_params;
    auto cli = (
      general_options(&vm_local),
      simulation_options(&VM),
      individual_options(&VM, &individual_params)
    );
    wtl::parse(cli, arguments);
    if (vm_local.at("help")) {
        auto fmt = wtl::doc_format();
        std::cout << "Usage: " << "PROJECT_NAME" << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local.at("version")) {
        std::cout << "PROJECT_VERSION" << "\n";
        throw wtl::ExitSuccess();
    }
    Individual::param(individual_params);
    wtl::sfmt64().seed(VM.at("seed").get<uint32_t>());
    const unsigned dimensions = VM.at("dimensions");
    const bool symmetric = VM.at("symmetric");
    if (dimensions == 1u) {
        // std::ostringstream ost;
        // ost << "diameter_pref = 1e6\n"
        //     << "limb_select = 1e6\n"
        //     << "mutation_mask = 10\n";
        // std::istringstream ist(ost.str());
        // po::store(po::parse_config_file(ist, description, false), vm);
    }
    if (symmetric) {
        // std::ostringstream ost;
        // ost << "diameter_pref = " << vm["height_pref"].as<double>() << "\n"
        //     << "limb_select = " << vm["toepad_select"].as<double>() << "\n";
        // std::istringstream ist(ost.str());
        // po::store(po::parse_config_file(ist, description, false), vm);
    }
    const std::string config = VM.dump(2) + "\n";
    verbose_ = vm_local.at("verbose");
    if (verbose_) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config << std::endl;
    }
    unsigned time = VM.at("time");
    unsigned interval = VM.at("interval");
    if (time % interval > 0u) {
        std::cerr << "T=" << time
                  << " is not a multiple of I="
                  << interval << std::endl;
        exit(1);
    }
    const std::string label = VM.at("label");
    fs::path OUT_DIR(VM.at("top_dir").get<std::string>());
    fs::create_directory(OUT_DIR);
    const std::string now(wtl::strftime("%Y%m%d_%H%M%S"));
    OUT_DIR /= (label + "_" + now + "_" + std::to_string(::getpid()));
    fs::create_directory(OUT_DIR);
    fs::current_path(OUT_DIR.string());
    wtl::make_ofs("config.json") << config;
}

void Simulation::run() {HERE;
    switch (VM.at("mode").get<int>()) {
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
    const unsigned time = VM.at("time");
    const unsigned interval = VM.at("interval");
    const unsigned nrow = VM.at("nrow");
    const unsigned ncol = VM.at("ncol");
    const unsigned dimensions = VM.at("dimensions");
    const unsigned patch_size = VM.at("patch_size");
    std::cout << nrow << "," << ncol << "," << dimensions << "," << patch_size << std::endl;
    assert((time % interval) == 0u);
    population.reserve(nrow);
    for (size_t row=0; row<nrow; ++row) {
        std::vector<Patch> pop_row;
        pop_row.reserve(ncol);
        for (size_t col=0; col<ncol; ++col) {
            pop_row.emplace_back(wtl::sfmt64()());
        }
        population.emplace_back(std::move(pop_row));
    }
    if (dimensions == 1u) {
        population[0][0].assign(patch_size, Individual{{15, 0, 15, 0}});
    } else {
        population[0][0].assign(patch_size, Individual{});
    }
    std::ostringstream ost;
    for (size_t t=0; t<=time; ++t) {
        if (verbose_) {
            std::cout << "\nT = " << t << "\n"
                << wtl::str_matrix(population, " ", wtl::make_oss(),
                        [](const Patch& p) {return p.size();})
                << std::flush;
        }
        if ((t % interval) == 0u) {
            write_snapshot(t, ost);
        }
        if (t < time) {
            life_cycle();
        }
    }
    wtl::zlib::ofstream{"evolution.csv.gz"} << ost.str();
}

void Simulation::life_cycle() {
    static wtl::ThreadPool pool(1u);
    const auto nrow = population.size();
    const auto ncol = population.front().size();
    const size_t num_patches = nrow * ncol;
    auto task = [this, nrow, ncol](size_t r, size_t c){
        auto& patch = this->population[r][c];
        auto children = patch.mate_and_reproduce();
        return std::pair<std::vector<Individual>, std::vector<std::pair<unsigned int, unsigned int>>>{
          children, patch.make_destinations(children.size(), r, c, nrow, ncol)
        };
    };
    std::vector<std::future<
      std::pair<std::vector<Individual>, std::vector<std::pair<unsigned int, unsigned int>>>
    >> futures;
    futures.reserve(num_patches);
    for (size_t row=0; row<nrow; ++row) {
        for (size_t col=0; col<ncol; ++col) {
            futures.emplace_back(pool.submit(task, row, col));
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
    for (size_t row=0; row<nrow; ++row) {
        for (size_t col=0; col<ncol; ++col) {
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
