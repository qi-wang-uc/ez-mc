#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "define.hpp"

class Config{
public:
    /* input */
    INT  nstep = 0; /* Total MC steps */
    INT  fsize = 0; /* frequency to update MC size */
    INT  fsavc = 0; /* frequency to save coordinate */
    INT  fsavl = 0; /* freqeuncy to save log file */
    STR  fasta; /* input file name for fasta sequence */
    REAL nsize = 0.0; /* Constant step size */
    REAL idacc = 0.5; /* Ideal accept rate for adaptive MC step size */

    /* output */
    STR log_name;
    STR psf_name;
    STR dcd_name;
    
    /* setter */
    Config () {}
    Config (const STR& inp_name);
    void print(void) const;
};

#endif