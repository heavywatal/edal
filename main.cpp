// -*- mode: c++; coding: utf-8 -*-
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream> // cstdio, string

#include "main.h"

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(16);
    std::cerr.precision(6);

//    check_flags(argc, argv);
    run();
    return EXIT_SUCCESS;

    return EXIT_FAILURE;
}

