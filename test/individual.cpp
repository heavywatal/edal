#include "individual.hpp"

#include <iostream>

int main() {
    std::cout.precision(15);
    edal::Individual ind;
    std::cout << ind.str_detail();
    return 0;
}
