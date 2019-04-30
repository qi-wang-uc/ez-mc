#ifndef MC_HPP
#define MC_HPP

#include <cmath>
#include "config.hpp"
#include "define.hpp"

/* Coefficients for Acceptance Ration Method */
const REAL COEFF_A = 0.672924;
const REAL COEFF_B = 0.064428;

/* Didn't implement too much. 
   Just for cleaner code. */
struct ResVec {
    INT  id;
    REAL x;
    REAL y;
    REAL z;
    ResVec () {}
    ResVec (INT id, REAL x, REAL y, REAL z):id(id), x(x),y(y),z(z) {}
};

class MCsystem {
    protected:
        CORARR X; // Cartesian X
        CORARR Y; // Cartesian Y
        CORARR Z; // Cartesian Z
        STRVEC S; // sequence
        INT    N; // system size
    public:
        MCsystem () {}
        MCsystem (const STR& inp_seq);
        /* getter : for block-wise data transfer */
        const float* x();
        const float* y();
        const float* z();
        float x(const INT& resid) const;
        float y(const INT& resid) const;
        float z(const INT& resid) const;
        /* getter : for regular query */
        const STRVEC& sequ(void) const;
        const STR& sequ(const INT& resid) const;
        const INT&    size(void) const;
        const bool has_overlap(const ResVec& inp_incr) const;
        const bool has_badbond(const ResVec& inp_incr) const;
        /* setter : update coordiante of a certain bead */
        void update(const ResVec& inp_incr);
};

void run_monte_carlo(MCsystem& mc_system, const Config& config);
bool pass_metropolis_crit(const REAL& e_new, const REAL& e_old);
REAL calc_stepsize(const REAL& delta_old, const REAL& p_old, const REAL& p_ideal);

#endif