#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include "define.hpp"

template<typename T>
void DEBUG(T msg) {
    std::cout << "DEBUG> " << msg << std::endl;
}

template<typename T>
STR encap(T input) {
    return "(" + std::to_string(input) + ")";
}

void time_elapsed(SysClock tick,  SysClock tock);

STR  time_stamp(void);

REAL calc_dist(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2);

#endif