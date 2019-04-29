#include <sstream>
#include <ctime>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include "../include/util.hpp"

void time_elapsed(SysClock tick, SysClock tock) {
    auto duration = tock - tick;
    auto seconds =  std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto total_seconds = seconds.count();
    auto n_hours   =  (total_seconds / 3600);
    auto n_minutes =  (total_seconds - n_hours*3600)/60;
    auto n_seconds =  (total_seconds - n_hours*3600 - n_minutes*60);
    std::cout << "EZ-MC> Time elapsed: "
              << encap(n_hours)   << " HR : "
              << encap(n_minutes) << " MIN : "
              << encap(n_seconds) << " SEC" << std::endl;
}

STR time_stamp(void) {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d %X");
    return ss.str();
}

REAL calc_dist(REAL x1, REAL y1, REAL z1, REAL x2, REAL y2, REAL z2) {
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}