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
        if("zstep"==arg) this->zstep = std::stof(val);
        if("nsavc"==arg) this->nsavc = std::stoi(val);
        if("nsavl"==arg) this->nsavl = std::stoi(val);
        if("fasta"==arg) this->fasta = val;
        if("logname"==arg) this->log_name = val;
        if("psfname"==arg) this->psf_name = val;
        if("dcdname"==arg) this->dcd_name = val;
    }
    inp_file.close();
    this->print();
}

void Config::print(void) const {
    std::cout << "RDCFG> --- Configuration for MC setup --- " << std::endl
              << "RDCFG> Number of MC steps: "   << this->nstep << std::endl
              << "RDCFG> Step size of MC move: " << this->zstep << std::endl
              << "RDCFG> Frequency to save coordinates: " << this->nsavc << std::endl
              << "RDCFG> Frequency to save log: " << this->nsavl << std::endl
              << "RDCFG> Input sequence file name: " << this->fasta << std::endl
              << "RDCFG> Output log file name: " << this->log_name << std::endl
              << "RDCFG> Output psf file name: " << this->psf_name << std::endl
              << "RDCFG> Output dcd file name: " << this->dcd_name << std::endl
              << "RDCFG> ----------------------------------" << std::endl
              << std::endl;
}