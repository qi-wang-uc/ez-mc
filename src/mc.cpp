#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include "../include/mc.hpp"
#include "../include/define.hpp"
#include "../include/io.hpp"
#include "../include/util.hpp"
#include "../include/energy.hpp"

MCsystem::MCsystem (const STR& inp_seq) {
    /* Read sequence and resize container */
    const auto size = inp_seq.size();
    this->X.resize(size);
    this->Y.resize(size);
    this->Z.resize(size);
    this->N = size;
    for (const auto& s : inp_seq) {
        this->S.push_back(AMINO_ACID_CODE.at(s));
    }
    std::cout << "MCSYS> Initialized with sequence (" << this->N << "):" << std::endl;
    INT counter = 0;
    for (const auto& s : this->S) {
        std::cout << s << " ";
        if(++counter%10==0) std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    /* Build a linear peptide along X-axis */
    for(counter=0; counter<size; ++counter) {
        this->X[counter] = counter*BONDLEN;
        this->Y[counter] = 0.0;
        this->Z[counter] = 0.0;
    }
}

const STRVEC& MCsystem::sequ(void) const {
    return this->S;
}

const STR& MCsystem::sequ(const INT& resid) const {
    return this->S.at(resid);
}

const INT& MCsystem::size(void) const {
    return this->N;
}

const bool MCsystem::has_overlap(const ResVec& inp_incr) const {
    const REAL xi = this->X[inp_incr.id] + inp_incr.x;
    const REAL yi = this->Y[inp_incr.id] + inp_incr.y;
    const REAL zi = this->Z[inp_incr.id] + inp_incr.z;
    for(INT i=0; i<this->N; ++i) {
        if (inp_incr.id==i) continue;
        if (calc_dist(xi, yi, zi, this->X[i], this->Y[i], this->Z[i])<= 2*RBEAD) {
            return true;
        }
    }
    return false;
}

const bool MCsystem::has_badbond(const ResVec& inp_incr) const {
    bool result;
    const REAL xi = this->X[inp_incr.id] + inp_incr.x;
    const REAL yi = this->Y[inp_incr.id] + inp_incr.y;
    const REAL zi = this->Z[inp_incr.id] + inp_incr.z;
    if(inp_incr.id==0) {
        const REAL d = calc_dist(xi, yi, zi, this->X[1], this->Y[1], this->Z[1]);
        result = (d-BONDMIN)*(d-BONDMAX) < 0 ? false : true;
    } else if (inp_incr.id==(this->N-1)) {
        const REAL d = calc_dist(xi, yi, zi, this->X[N-2], this->Y[N-2], this->Z[N-2]);
        result = (d-BONDMIN)*(d-BONDMAX) < 0 ? false : true;
    } else {
        const REAL d1 = calc_dist(xi, yi, zi, this->X[inp_incr.id-1], this->Y[inp_incr.id-1], this->Z[inp_incr.id-1]);
        const REAL d2 = calc_dist(xi, yi, zi, this->X[inp_incr.id+1], this->Y[inp_incr.id+1], this->Z[inp_incr.id+1]);
        result = ((d1-BONDMIN)*(d1-BONDMAX)<0) && ((d2-BONDMIN)*(d2-BONDMAX)<0) ? false : true;
    }
    return result;
}

void MCsystem::update(const ResVec& inp_incr) {
    this->X[inp_incr.id] += inp_incr.x;
    this->Y[inp_incr.id] += inp_incr.y;
    this->Z[inp_incr.id] += inp_incr.z;
}

const float* MCsystem::x() {
    return this->X.data();
}
const float* MCsystem::y() {
    return this->Y.data();
}
const float* MCsystem::z() {
    return this->Z.data();
}

float MCsystem::x(const INT& resid) const {
    return this->X.at(resid);
}
float MCsystem::y(const INT& resid) const {
    return this->Y.at(resid);
}
float MCsystem::z(const INT& resid) const {
    return this->Z.at(resid);
}

void run_monte_carlo(MCsystem& mc_system, const Config& config) {
    /* Initialize output files with headers */
    std::cout << "RUNMC> Starting Monte Carlo simulation." << std::endl << std::endl;
    std::ofstream log_file(config.log_name);
    log_file << std::setw(12) << "#Step" << " "
             << std::setw(12) << "Trials" << " "
             << std::setw(12) << "Energy" << " "
             << std::setw(12) << "Step size" << " "
             << std::setw(12) << "Accept raio" << " "
             << std::endl;
    std::ofstream dcd_file(config.dcd_name, std::ios::binary);
    write_dcdheader( dcd_file, DCD_Header(
            static_cast<int32_t>(mc_system.size()),
            static_cast<int32_t>(config.nstep/config.fsavc),
            static_cast<int32_t>(0),
            static_cast<int32_t>(config.fsavc),
            static_cast<int32_t>(config.nstep),
            static_cast<int32_t>(0),
            static_cast<float>(1.0),
            static_cast<const char*>(("REMARK CREATED AT " + time_stamp()).c_str()),
            static_cast<const char*>(("REMARK CREATED BY " + STR(getenv("USER"))).c_str()))
    );
    /* Set up random numbers */
    std::mt19937_64 rand_gen(time(0));
    std::uniform_int_distribution<>  ud_int(0, mc_system.size()-1);
    std::uniform_real_distribution<> ud_real(0.0, 1.0);
    /* MC begin */
    REAL e_new, e_old, p_accept, step_size=1.0/MATH_E;
    bool is_adaptive_size = true;
    if (config.nsize > 0.0) {
        is_adaptive_size = false;
        step_size = config.nsize;
    }
    INT  overall_trial_counter = 0;
    for(INT istep=1; istep<=config.nstep; ++istep) {
        std::cout << "\e[A" <<"RUNMC> "
                  << std::fixed << std::setw(6) << std::setprecision(2)
                  << istep*100.0 / config.nstep << "% Done."
                  << std::endl;

        /* Generate a random move for a random bead */
        ResVec incr_vec;
        INT  step_trial_counter = 0;
        bool success_trial = false;
        while(!success_trial) {
            /* Generate a random bead and count its trials */
            INT rand_bead = ud_int(rand_gen); 
            INT bead_trial_counter = 0;
            while(bead_trial_counter++ < MAXTRIAL) {
                auto phi   = 2.0*MATHPI*ud_real(rand_gen);
                auto theta = acos(2.0*ud_real(rand_gen)-1.0);
                incr_vec   = ResVec(rand_bead, step_size*sin(theta)*cos(phi), 
                                               step_size*sin(theta)*sin(phi), 
                                               step_size*cos(theta));
                if (mc_system.has_overlap(incr_vec)) continue;
                if (mc_system.has_badbond(incr_vec)) continue;
                e_new = calc_energy(mc_system);
                if (pass_metropolis_crit(e_new, e_old)) {
                    success_trial = true;
                    break;
                }
            }
            step_trial_counter += bead_trial_counter;
        }
        overall_trial_counter += step_trial_counter;
        e_old = e_new;
        mc_system.update(incr_vec);

        /* Update MC step size */
        if(0==istep%config.fsize && is_adaptive_size) {
            p_accept  = 1.0*config.fsize/overall_trial_counter;
            step_size = calc_stepsize(step_size, p_accept, config.idacc);
            overall_trial_counter = 0;
        }

        /* Write binary coordinates */
        if(0==istep%config.fsavc) {
            write_dcdframe(dcd_file, mc_system.size(),
                           mc_system.x(), mc_system.y(), mc_system.z());
        }

        /* Write log file upon request */
        if(0==istep%config.fsavl) {
            log_file << std::right << std::fixed
                     << std::setw(12) << istep << " "
                     << std::setw(12) << step_trial_counter << " "
                     << std::setprecision(6)
                     << std::setw(12) << e_old << " "
                     << std::setprecision(6)
                     << std::setw(12) << step_size << " "
                     << std::setprecision(6)
                     << std::setw(12) << p_accept  << " "
                     << std::endl;
        }
    }
    log_file.close();
    dcd_file.close();
}

bool pass_metropolis_crit(const REAL& e_new, const REAL& e_old) {
    if (e_new < e_old) return true;
    std::default_random_engine rand_gen;
    std::uniform_real_distribution<> ud_real(0.0, 1.0);
    return std::exp(BETA*(e_new-e_old)) > ud_real(rand_gen);
}

REAL calc_stepsize(const REAL& delta_old, const REAL& p_old, const REAL& p_ideal) {
    return delta_old*(std::log(COEFF_A*p_ideal+COEFF_B)/std::log(COEFF_A*p_old+COEFF_B));
}