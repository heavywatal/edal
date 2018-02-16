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

void Patch::change_sex_half(size_t n) {
    n /= 2;
    for (auto it=members_.begin()+n; it!=members_.end(); ++it) {
        it->change_sex();
    }
}

std::vector<Individual> Patch::mate_and_reproduce(URNG& rng) const {
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
        std::uniform_real_distribution<> uniform(0.0, upper_bounds.back());
        const double dart = uniform(rng);
        size_t father_i = 0;
        while (upper_bounds[father_i] < dart) {++father_i;}
        const Individual& father = members_[male_indices[father_i]];
        const size_t num_children = poisson(rng);
        for (size_t i=0; i<num_children; ++i) {
            offsprings.emplace_back(
                mother.gametogenesis(rng), father.gametogenesis(rng),
                std::bernoulli_distribution(0.5)(rng));
        }
    }
    return offsprings;
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

void Patch::viability_selection(URNG& rng) {
    std::vector<size_t> indices;
    indices.reserve(members_.size());
    const auto ne = effective_num_competitors();
    for (size_t i=0; i<members_.size(); ++i) {
        const double p = members_[i].survival_probability(ne[i]);
        if (std::bernoulli_distribution(p)(rng)) {
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

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Patch::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    Patch patch(20, Individual({15,0,15,0}));
    std::cerr << patch.size();
    for (size_t i=0; i<10; ++i) {
        for (auto& child: patch.mate_and_reproduce(wtl::sfmt())) {
            patch.append(std::move(child));
        }
        std::cerr << " b " << patch.size();
        patch.viability_selection(wtl::sfmt());
        std::cerr << " d " << patch.size();
    }
    std::cerr << std::endl;
    std::cerr << patch << std::endl;
}

} // namespace edal
