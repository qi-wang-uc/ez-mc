#include <iostream>
#include <fstream>
#include <sstream>
#include "../include/config.hpp"

Config::Config(const STR& inp_name) {
    std::cout << "RDCFG> Reading configure from file: " << inp_name << std::endl;
    std::ifstream inp_file(inp_name);
    STR arg, val;
    STR each_line;
    while(getline(inp_file, each_line)) {
        // Remove leading spaces
        each_line.erase(0, each_line.find_first_not_of(" \t"));
        if('#'==each_line.front()) continue;
        if(each_line.empty()) continue;
        std::stringstream each_stream(each_line);
        each_stream >> arg >> val;
        if("nstep"==arg) this->nstep = std::stoi(val);
        if("nsize"==arg) this->nsize = std::stof(val);
        if("fsavc"==arg) this->fsavc = std::stoi(val);
        if("fsavl"==arg) this->fsavl = std::stoi(val);
        if("fasta"==arg) this->fasta = val;
        if("logname"==arg) this->log_name = val;
        if("psfname"==arg) this->psf_name = val;
        if("dcdname"==arg) this->dcd_name = val;
        if("fsize"==arg) this->fsize = std::stoi(val);
        if("idacc"==arg) this->idacc = std::stof(val);
    }
    inp_file.close();
    this->print();
}

void Config::print(void) const {
    std::cout << "RDCFG> --- Configuration for MC setup --- " << std::endl
              << "RDCFG> Number of MC steps: "   << this->nstep << std::endl
              << "RDCFG> Step size of MC move: " << this->nsize << std::endl
              << "RDCFG> Frequency to save coordinates: " << this->fsavc << std::endl
              << "RDCFG> Frequency to save log: " << this->fsavl << std::endl
              << "RDCFG> Input sequence file name: " << this->fasta << std::endl
              << "RDCFG> Output log file name: " << this->log_name << std::endl
              << "RDCFG> Output psf file name: " << this->psf_name << std::endl
              << "RDCFG> Output dcd file name: " << this->dcd_name << std::endl
              << "RDCFG> ----------------------------------" << std::endl
              << std::endl;
}