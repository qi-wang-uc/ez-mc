#ifndef DEFINE_HPP
#define DEFINE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

using INT  = size_t;

using REAL = float;

using STR  = std::string;

using STRVEC = std::vector<STR>;

using CORARR = std::vector<REAL>;

using TimeStamp = std::chrono::system_clock;

using SysClock  = std::chrono::time_point<std::chrono::system_clock>;

/* Average, Minimum and Maximum  values of bond length 
   obtained from 1 ns equilibration of example protein
   1PGB.pdb at atomistic resolution. */
const REAL BONDMIN  = 3.6;
const REAL BONDLEN  = 3.8;
const REAL BONDMAX  = 4.0;

/* [?] Bead radius, needs to be verified. */
const REAL RBEAD    = 1.8;

/* Cut-off radius to calculate energy. Same value 
   as reference paper for residue-wise contact energy. */
const REAL RCUTOFF = 6.5;

/* Maximum iteration for a move. Going beyond this number 
   will result in a new choice of bead */
const INT  MAXTRIAL = 1000000;

/* Other math and physics constants */
const REAL MATHPI  = 3.14159;
const REAL MATH_E  = 2.71828;
const REAL BETA    = 1.0;

#endif