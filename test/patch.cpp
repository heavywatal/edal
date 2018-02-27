#include "patch.hpp"
#include "individual.hpp"

#include <random>
#include <iostream>

int main() {
    edal::Patch patch(std::random_device{}());
    patch.assign(20, edal::Individual({15,0,15,0}));
    std::cout << patch.size();
    for (size_t i=0; i<10; ++i) {
        for (auto& child: patch.mate_and_reproduce()) {
            patch.emplace_back(std::move(child));
        }
        std::cout << " b " << patch.size();
        patch.viability_selection();
        std::cout << " d " << patch.size();
    }
    std::cout << std::endl;
    std::cout << patch << std::endl;
    return 0;
}
