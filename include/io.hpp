#ifndef IO_HPP
#define IO_HPP
#include <fstream>
#include "define.hpp"

/* FASTA, PSF and DCD files I/O */

struct DCD_Header {
    int32_t natom;
    int32_t nfile;
    int32_t npriv;
    int32_t nsavc;
    int32_t nstep;
    int32_t ifpbc;
    float delta;
    const char* prog_title; 
    const char* user_title; 
    DCD_Header() {}
    DCD_Header(int32_t natom, int32_t nfile, int32_t npriv, 
               int32_t nsavc, int32_t nstep, int32_t ifpbc, 
               float delta, const char* prog_title, const char* user_title ):
               natom(natom),nfile(nfile),npriv(npriv),
               nsavc(nsavc),nstep(nstep),ifpbc(ifpbc),
               delta(delta),prog_title(prog_title), user_title(user_title) {}
};

const struct DCD_Pads {
    const int32_t pad0   = 0;
    const int32_t pad2   = 2;
    const int32_t pad4   = 4;
    const int32_t pad24  = 24;
    const int32_t pad84  = 84;
    const int32_t pad164 = 164;  
} dcd_pads;

STR read_fasta(const STR& inp_name);

void write_psf(const STR& out_name, const STRVEC& inp_seq);

void write_dcdheader(std::ofstream& dcd_file, const DCD_Header& dcd_header);

void write_dcdframe(std::ofstream& dcd_file, const INT& natom, 
                const float* fX, const float* fY, const float* fZ);

#endif