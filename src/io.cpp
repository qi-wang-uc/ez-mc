#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/io.hpp"

STR read_fasta(const STR& inp_name) {
    std::cout << "RDSEQ> Reading sequence from file: " << inp_name << std::endl;
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) std::cerr << "ERROR> Cannot open file." << std::endl;
    STR each_line;
    STR sequence;
    while(getline(inp_file, each_line)) {
        if (each_line.empty()) continue;
        if (each_line.front()=='>') continue;
        sequence += each_line;
    }
    inp_file.close();
    std::cout << std::endl;
    return sequence;
}

void write_psf(const STR& out_name, const STRVEC& inp_seq) {
    std::cout << "PSFIO> Writting dummy PSF file to: " << out_name << std::endl;
    std::ofstream out_file(out_name);
    const auto NATOM = inp_seq.size();
    /* title */
    out_file << "PSF CMAP CHEQ" << std::endl
             << std::endl
             << std::setw(8) << 1 << " !NTITLE" << std::endl
             << "* DUMMY PSF FILE FOR MC TRAJECTORY" << std::endl
             << std::endl
             << std::setw(8) << NATOM << " !NATOM" << std::endl;
    /* natom */
    for(INT i=0; i<NATOM; ++i) {
        out_file << std::right
                 << std::setw(8) << i+1 << " "          // %8d_
                 << std::left
                 << std::setw(4) << "EZMC" << " "       // %-4s_
                 << std::setw(4) << i+1 << " "          // %-4i_
                 << std::setw(4) << inp_seq[i] << " "   // %-4s_                             
                 << std::setw(4) << "CA" << " "         // %-4s_
                 << std::setw(4) << "CA" << " "         // %-4s_
                 << std::right
                 << std::fixed << std::setw(10) 
                 << std::setprecision(6) << 0.0 << " "  // %10.6f_
                 << std::fixed << std::setw(13) 
                 << std::setprecision(4) << 110.0       // %13.6f
                 << "           "                       // 11*' '
                 << std::endl;
    }
    /* nbond */
    out_file << std::endl
             << std::setw(8) << NATOM-1 << " !NBOND" << std::endl;
    for(INT i=1; i<NATOM; ++i) {
        out_file << std::setw(8) << i << std::setw(8) << i+1;
        if(0==i%4) out_file << std::endl;
    }
    out_file  << std::endl;
    std::cout << std::endl;
    out_file.close();

}

void write_dcdheader(std::ofstream& dcd_file, const DCD_Header& dcd_header) {
    /* dcd header part 1. */
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad84), 4);
    dcd_file.write("CORD", 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nfile), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.npriv), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nsavc), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nstep), 4);
    for(int i=0; i<5; ++i) {
        dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad0), 4);
    }
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.delta), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.ifpbc), 4);
    for(int i=0; i<8; ++i) {
        dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad0), 4);
    }
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad24), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad84), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad2), 4);
    /* dcd header part 2. */
    dcd_file.write(dcd_header.prog_title, 80);
    dcd_file.write(dcd_header.user_title, 80);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.natom), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
}

void write_dcdframe(std::ofstream& dcd_file, const INT& natom, 
                const float* fX, const float* fY, const float* fZ) {
    const int32_t pad4N = 4*natom;
    /* write x coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fX), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write y coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fY), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write z coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fZ), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
}