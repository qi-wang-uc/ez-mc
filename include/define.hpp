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

/* Bond length */
const REAL BONDLEN  = 3.8;

/* [?] Minimum allowed bond length */
const REAL BONDMIN  = 2.8;

/* [?] Maximum allowed bond length */
const REAL BONDMAX  = 4.8;

/* [?] Bead radius */
const REAL RBEAD    = 1.8;

/* Cut-off radius to calculate energy */
const REAL RCUTOFF = 6.5;

/* Math Pi */
const REAL MATHPI  = 3.14159;

/* Maximum iteration for a move. Going beyond this number 
   will result in a new choice of bead */
const INT  MAXTRIAL = 1000000;

#endif