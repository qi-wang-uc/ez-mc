#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <unordered_map>
#include "define.hpp"
#include "mc.hpp"

struct ResPair{
  STR res1;
  STR res2;
  ResPair(STR res1, STR res2):res1(res1),res2(res2) {}
  bool operator==(const ResPair& rp) const {
    return this->res1==rp.res1 && this->res2==rp.res2;
  }
};

struct HashResPair {
  INT operator()(const ResPair& rp) const {
    std::hash<STR> h1, h2;
    return h1(rp.res1) ^ h2(rp.res2);
  }
};

const std::unordered_map<ResPair, REAL, HashResPair> CONTACT_ENERGY = {
  /* CYS */
  {{"CYS", "CYS"}, -5.44},
  /* MET */
  {{"MET", "CYS"}, -5.05},
  {{"MET", "MET"}, -6.06},
  /* PHE */
  {{"PHE", "CYS"}, -5.63},
  {{"PHE", "MET"}, -6.68},
  {{"PHE", "PHE"}, -6.85},
  /* ILE */
  {{"ILE", "CYS"}, -5.03},
  {{"ILE", "MET"}, -6.33},
  {{"ILE", "PHE"}, -6.39},
  {{"ILE", "ILE"}, -6.22},
  /* LEU */
  {{"LEU", "CYS"}, -5.03},
  {{"LEU", "MET"}, -6.01},
  {{"LEU", "PHE"}, -6.26},
  {{"LEU", "ILE"}, -6.17},
  {{"LEU", "LEU"}, -5.79},
  /* VAL */
  {{"VAL", "CYS"}, -4.46},
  {{"VAL", "MET"}, -5.52},
  {{"VAL", "PHE"}, -5.75},
  {{"VAL", "ILE"}, -5.58},
  {{"VAL", "LEU"}, -5.38},
  {{"VAL", "VAL"}, -4.94},
  /* TRP */
  {{"TRP", "CYS"}, -4.76},
  {{"TRP", "MET"}, -6.37},
  {{"TRP", "PHE"}, -6.02},
  {{"TRP", "ILE"}, -5.64},
  {{"TRP", "LEU"}, -5.50},
  {{"TRP", "VAL"}, -5.05},
  {{"TRP", "TRP"}, -5.42},
  /* TYR */
  {{"TYR", "CYS"}, -3.89},
  {{"TYR", "MET"}, -4.92},
  {{"TYR", "PHE"}, -4.95},
  {{"TYR", "ILE"}, -4.63},
  {{"TYR", "LEU"}, -4.26},
  {{"TYR", "VAL"}, -4.05},
  {{"TYR", "TRP"}, -4.44},
  {{"TYR", "TYR"}, -3.55},
  /* ALA */
  {{"ALA", "CYS"}, -3.38},
  {{"ALA", "MET"}, -3.99},
  {{"ALA", "PHE"}, -4.36},
  {{"ALA", "ILE"}, -4.41},
  {{"ALA", "LEU"}, -3.96},
  {{"ALA", "VAL"}, -3.62},
  {{"ALA", "TRP"}, -3.93},
  {{"ALA", "TYR"}, -2.85},
  {{"ALA", "ALA"}, -2.51},
  /* GLY */
  {{"GLY", "CYS"}, -3.16},
  {{"GLY", "MET"}, -3.75},
  {{"GLY", "PHE"}, -3.72},
  {{"GLY", "ILE"}, -3.65},
  {{"GLY", "LEU"}, -3.43},
  {{"GLY", "VAL"}, -3.06},
  {{"GLY", "TRP"}, -3.37},
  {{"GLY", "TYR"}, -2.50},
  {{"GLY", "ALA"}, -2.15},
  {{"GLY", "GLY"}, -2.17},
  /* THR */
  {{"THR", "CYS"}, -2.88},
  {{"THR", "MET"}, -3.73},
  {{"THR", "PHE"}, -3.76},
  {{"THR", "ILE"}, -3.74},
  {{"THR", "LEU"}, -3.43},
  {{"THR", "VAL"}, -2.95},
  {{"THR", "TRP"}, -3.31},
  {{"THR", "TYR"}, -2.48},
  {{"THR", "ALA"}, -2.15},
  {{"THR", "GLY"}, -2.03},
  {{"THR", "THR"}, -1.72},
  /* SER */
  {{"SER", "CYS"}, -2.86},
  {{"SER", "MET"}, -3.55},
  {{"SER", "PHE"}, -3.56},
  {{"SER", "ILE"}, -3.43},
  {{"SER", "LEU"}, -3.16},
  {{"SER", "VAL"}, -2.79},
  {{"SER", "TRP"}, -2.95},
  {{"SER", "TYR"}, -2.30},
  {{"SER", "ALA"}, -1.89},
  {{"SER", "GLY"}, -1.70},
  {{"SER", "THR"}, -1.59},
  {{"SER", "SER"}, -1.48},
  /* GLN */
  {{"GLN", "CYS"}, -2.73},
  {{"GLN", "MET"}, -3.17},
  {{"GLN", "PHE"}, -3.30},
  {{"GLN", "ILE"}, -3.22},
  {{"GLN", "LEU"}, -3.09},
  {{"GLN", "VAL"}, -2.67},
  {{"GLN", "TRP"}, -3.16},
  {{"GLN", "TYR"}, -2.53},
  {{"GLN", "ALA"}, -1.70},
  {{"GLN", "GLY"}, -1.54},
  {{"GLN", "THR"}, -1.59},
  {{"GLN", "SER"}, -1.37},
  {{"GLN", "GLN"}, -0.89},
  /* ASN */
  {{"ASN", "CYS"}, -2.59},
  {{"ASN", "MET"}, -3.50},
  {{"ASN", "PHE"}, -3.55},
  {{"ASN", "ILE"}, -2.99},
  {{"ASN", "LEU"}, -2.99},
  {{"ASN", "VAL"}, -2.36},
  {{"ASN", "TRP"}, -3.11},
  {{"ASN", "TYR"}, -2.47},
  {{"ASN", "ALA"}, -1.44},
  {{"ASN", "GLY"}, -1.56},
  {{"ASN", "THR"}, -1.51},
  {{"ASN", "SER"}, -1.31},
  {{"ASN", "GLN"}, -1.36},
  {{"ASN", "ASN"}, -1.59},
  /* GLU */
  {{"GLU", "CYS"}, -2.08},
  {{"GLU", "MET"}, -3.19},
  {{"GLU", "PHE"}, -3.51},
  {{"GLU", "ILE"}, -3.23},
  {{"GLU", "LEU"}, -2.91},
  {{"GLU", "VAL"}, -2.56},
  {{"GLU", "TRP"}, -2.94},
  {{"GLU", "TYR"}, -2.42},
  {{"GLU", "ALA"}, -1.51},
  {{"GLU", "GLY"}, -1.22},
  {{"GLU", "THR"}, -1.45},
  {{"GLU", "SER"}, -1.48},
  {{"GLU", "GLN"}, -1.33},
  {{"GLU", "ASN"}, -1.43},
  {{"GLU", "GLU"}, -1.18},
  /* ASP */
  {{"ASP", "CYS"}, -2.66},
  {{"ASP", "MET"}, -2.90},
  {{"ASP", "PHE"}, -3.31},
  {{"ASP", "ILE"}, -2.91},
  {{"ASP", "LEU"}, -2.59},
  {{"ASP", "VAL"}, -2.25},
  {{"ASP", "TRP"}, -2.91},
  {{"ASP", "TYR"}, -2.25},
  {{"ASP", "ALA"}, -1.57},
  {{"ASP", "GLY"}, -1.62},
  {{"ASP", "THR"}, -1.66},
  {{"ASP", "SER"}, -1.46},
  {{"ASP", "GLN"}, -1.26},
  {{"ASP", "ASN"}, -1.33},
  {{"ASP", "GLU"}, -1.23},
  {{"ASP", "ASP"}, -0.96},
  /* HIS */
  {{"HIS", "CYS"}, -3.63},
  {{"HIS", "MET"}, -3.31},
  {{"HIS", "PHE"}, -4.61},
  {{"HIS", "ILE"}, -3.76},
  {{"HIS", "LEU"}, -3.84},
  {{"HIS", "VAL"}, -3.38},
  {{"HIS", "TRP"}, -4.02},
  {{"HIS", "TYR"}, -3.33},
  {{"HIS", "ALA"}, -2.09},
  {{"HIS", "GLY"}, -1.94},
  {{"HIS", "THR"}, -2.31},
  {{"HIS", "SER"}, -1.94},
  {{"HIS", "GLN"}, -1.85},
  {{"HIS", "ASN"}, -2.01},
  {{"HIS", "GLU"}, -2.27},
  {{"HIS", "ASP"}, -2.14},
  {{"HIS", "HIS"}, -2.78},
  /* ARG */
  {{"ARG", "CYS"}, -2.70},
  {{"ARG", "MET"}, -3.49},
  {{"ARG", "PHE"}, -3.54},
  {{"ARG", "ILE"}, -3.33},
  {{"ARG", "LEU"}, -3.15},
  {{"ARG", "VAL"}, -2.78},
  {{"ARG", "TRP"}, -3.56},
  {{"ARG", "TYR"}, -2.75},
  {{"ARG", "ALA"}, -1.50},
  {{"ARG", "GLY"}, -1.68},
  {{"ARG", "THR"}, -1.97},
  {{"ARG", "SER"}, -1.22},
  {{"ARG", "GLN"}, -1.85},
  {{"ARG", "ASN"}, -1.41},
  {{"ARG", "GLU"}, -2.07},
  {{"ARG", "ASP"}, -1.98},
  {{"ARG", "HIS"}, -2.12},
  {{"ARG", "ARG"}, -1.39},
  /* LYS */
  {{"LYS", "CYS"}, -1.54},
  {{"LYS", "MET"}, -3.11},
  {{"LYS", "PHE"}, -2.83},
  {{"LYS", "ILE"}, -2.70},
  {{"LYS", "LEU"}, -2.63},
  {{"LYS", "VAL"}, -1.95},
  {{"LYS", "TRP"}, -2.49},
  {{"LYS", "TYR"}, -2.01},
  {{"LYS", "ALA"}, -1.10},
  {{"LYS", "GLY"}, -0.84},
  {{"LYS", "THR"}, -1.02},
  {{"LYS", "SER"}, -0.83},
  {{"LYS", "GLN"}, -1.02},
  {{"LYS", "ASN"}, -0.91},
  {{"LYS", "GLU"}, -1.60},
  {{"LYS", "ASP"}, -1.32},
  {{"LYS", "HIS"}, -1.09},
  {{"LYS", "ARG"}, -0.06},
  {{"LYS", "LYS"},  0.13},
  /* PRO */
  {{"PRO", "CYS"}, -2.92},
  {{"PRO", "MET"}, -4.11},
  {{"PRO", "PHE"}, -3.73},
  {{"PRO", "ILE"}, -3.47},
  {{"PRO", "LEU"}, -3.06},
  {{"PRO", "VAL"}, -2.96},
  {{"PRO", "TRP"}, -3.66},
  {{"PRO", "TYR"}, -2.80},
  {{"PRO", "ALA"}, -1.81},
  {{"PRO", "GLY"}, -1.72},
  {{"PRO", "THR"}, -1.66},
  {{"PRO", "SER"}, -1.35},
  {{"PRO", "GLN"}, -1.73},
  {{"PRO", "ASN"}, -1.43},
  {{"PRO", "GLU"}, -1.40},
  {{"PRO", "ASP"}, -1.19},
  {{"PRO", "HIS"}, -2.17},
  {{"PRO", "ARG"}, -1.85},
  {{"PRO", "LYS"}, -0.67},
  {{"PRO", "PRO"}, -1.18},
};

const std::unordered_map<char, STR> AMINO_ACID_CODE = {
  {'A', "ALA"},
  {'C', "CYS"},
  {'D', "ASP"},
  {'E', "GLU"},
  {'F', "PHE"},
  {'G', "GLY"},
  {'H', "HIS"},
  {'I', "ILE"},
  {'K', "LYS"},
  {'L', "LEU"},
  {'M', "MET"},
  {'N', "ASN"},
  {'P', "PRO"},
  {'Q', "GLN"},
  {'R', "ARG"},
  {'S', "SER"},
  {'T', "THR"},
  {'V', "VAL"},
  {'W', "TRP"},
  {'Y', "TYR"}
};

const REAL get_eij(const ResPair& respair);

const REAL calc_energy(const MCsystem& mc_system);

#endif