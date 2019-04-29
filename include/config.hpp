#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "define.hpp"

class Config{
public:
    /* input */
    INT  nstep;
    INT  nsavc;
    INT  nsavl;
    REAL zstep;
    STR  fasta;

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