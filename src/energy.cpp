#include "../include/define.hpp"
#include "../include/energy.hpp"
#include "../include/util.hpp"

const REAL get_eij(const ResPair& respair) {
    return CONTACT_ENERGY.find(respair)!=CONTACT_ENERGY.end() ? 
           CONTACT_ENERGY.at(respair) : 
           CONTACT_ENERGY.at(ResPair(respair.res2, respair.res1));
}

const REAL calc_energy(const MCsystem& mc_system) {
    REAL result = 0.0;
    for(INT ires=0; ires<mc_system.size()-1; ++ires) {
        const auto typei = mc_system.sequ(ires);
        const auto xi    = mc_system.x(ires);
        const auto yi    = mc_system.y(ires);
        const auto zi    = mc_system.z(ires);
        for(INT jres=ires+1; jres<mc_system.size()-2; ++jres) {
            const auto typej = mc_system.sequ(jres);
            const auto xj    = mc_system.x(jres);
            const auto yj    = mc_system.y(jres);
            const auto zj    = mc_system.z(jres);
            const auto d     = calc_dist(xi, yi, zi, xj, yj, zj);
            result += (d < RCUTOFF) ? get_eij(ResPair(typei, typej)) : 0.0;
        }
    }
    return result;
}