/*! @file patch.cpp
    @brief Implementation of Patch class
*/
#include "patch.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>
#include <sfmt.hpp>
#include <boost/program_options.hpp>

#include <cmath>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace edal {

Patch::Patch(unsigned int seed)
: engine_(std::make_unique<wtl::sfmt19937_64>(seed)) {}

Patch::Patch(Patch&& other) noexcept
: members_(std::move(other.members_)),
  engine_(std::move(other.engine_)) {}

Patch::~Patch() {}

void Patch::assign(size_t n, const Individual& founder) {
    members_.assign(n, founder);
    change_sex_half(n);
}

void Patch::change_sex_half(size_t n) {
    n /= 2;
    for (auto it=members_.begin()+n; it!=members_.end(); ++it) {
        it->change_sex();
    }
}

std::vector<Individual> Patch::mate_and_reproduce() const {
    std::poisson_distribution<size_t> poisson(Individual::AVG_NUM_OFFSPINRGS());
    std::vector<Individual> offsprings;
    std::vector<size_t> male_indices;
    male_indices.reserve(members_.size());
    for (size_t i=0; i<members_.size(); ++i) {
        if (members_[i].is_male()) {
            male_indices.push_back(i);
        }
    }
    if (male_indices.empty()) {return offsprings;}
    offsprings.reserve(members_.size() * Individual::AVG_NUM_OFFSPINRGS());
    for (const auto& mother: members_) {
        if (mother.is_male()) continue;
        std::vector<double> prefs;
        prefs.reserve(male_indices.size());
        for (const size_t i: male_indices) {
            prefs.push_back(mother.mating_probability_TPG2013(members_[i]));
        }
        std::vector<double> upper_bounds(male_indices.size());
        std::partial_sum(prefs.begin(), prefs.end(), upper_bounds.begin());
        std::uniform_real_distribution<double> uniform(0.0, upper_bounds.back());
        const double dart = uniform(*engine_);
        size_t father_i = 0;
        while (upper_bounds[father_i] < dart) {++father_i;}
        const Individual& father = members_[male_indices[father_i]];
        const size_t num_children = poisson(*engine_);
        for (size_t i=0; i<num_children; ++i) {
            offsprings.emplace_back(
                mother.gametogenesis(*engine_),
                father.gametogenesis(*engine_),
                std::bernoulli_distribution(0.5)(*engine_));
        }
    }
    return offsprings;
}

inline std::pair<unsigned int, unsigned int> choose_patch(size_t row, size_t col, Patch::URBG& engine) {
    thread_local std::bernoulli_distribution bern_mig(Individual::MIGRATION_RATE());
    thread_local std::uniform_int_distribution<unsigned int> unif_int(0, 7);
    if (bern_mig(engine)) {
        switch (unif_int(engine)) {
          case 0:        ++col; break;
          case 1: ++row; ++col; break;
          case 2: ++row;        break;
          case 3: ++row; --col; break;
          case 4:        --col; break;
          case 5: --row; --col; break;
          case 6: --row;        break;
          case 7: --row; ++col; break;
        }
    }
    return {row, col};
}

std::vector<std::pair<unsigned int, unsigned int>>
Patch::make_destinations(size_t n, size_t row, size_t col, size_t num_rows, size_t num_cols) const {
    std::vector<std::pair<unsigned int, unsigned int> > destinations;
    destinations.reserve(n);
    for (size_t i=0; i<n; ++i) {
        const auto dst = choose_patch(row, col, *engine_);
        if ((dst.first >= num_rows) | (dst.second >= num_cols)) {
            destinations.emplace_back(row, col);
        } else {
            destinations.push_back(std::move(dst));
        }
    }
    return destinations;
}

std::vector<double> Patch::effective_num_competitors() const {
    std::vector<double> results(members_.size());
    for (size_t i=0; i<members_.size(); ++i) {
        results[i] += 1.0;
        for (size_t j=0; j<i; ++j) {
            const double value = members_[i].resource_overlap(members_[j]);
            results[i] += value;
            results[j] += value;
        }
    }
    return results;
}

void Patch::viability_selection() {
    std::vector<size_t> indices;
    indices.reserve(members_.size());
    const auto ne = effective_num_competitors();
    for (size_t i=0; i<members_.size(); ++i) {
        const double p = members_[i].survival_probability(ne[i]);
        if (std::bernoulli_distribution(p)(*engine_)) {
            indices.push_back(i);
        }
    }
    std::vector<Individual> tmp;
    tmp.reserve(members_.size());
    for (auto i: indices) {
        tmp.push_back(std::move(members_[i]));
    }
    members_.swap(tmp);
}

std::map<Individual, size_t> Patch::summarize() const {
    std::map<Individual, size_t> genotypes;
    for (const auto& ind: members_) {
        ++genotypes[ind];
    }
    return genotypes;
}

std::ostream& operator<< (std::ostream& ost, const Patch& patch) {
    return ost << patch.summarize();
    // operator<< for std::map is defined in "iostr.hpp"
}

} // namespace edal
