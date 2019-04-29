#include <iostream>
#include "../include/define.hpp"
#include "../include/io.hpp"
#include "../include/mc.hpp"
#include "../include/util.hpp"
#include "../include/config.hpp"

int main(int argc, char* argv[]) {
    
    std::ios::sync_with_stdio(false);

    std::cout << "===== EZ-MC version 1.0 ====="
              << std::endl << std::endl;

    if (1==argc) {
        std::cerr << "ERROR> Missing input file." << std::endl;
        return -1;
    }

    auto tick = TimeStamp::now();

    const auto config = Config(argv[1]);

    MCsystem mc_system(read_fasta(config.fasta));

    write_psf(config.psf_name, mc_system.sequ());

    run_monte_carlo(mc_system, config);

    auto tock = TimeStamp::now();

	time_elapsed(tick, tock);
    
    return 0;
}