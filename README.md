# ez-mc
Quick and easy Monte Carlo simulations

Example input file (`config.cfg`):
```
# MC steps (int)
nstep       100000

# Step size of a MC move, need to tweak (float)
zstep       0.1

# freqeuncy to save coordinates (int)
nsavc       1000

# freqeuncy to write log file (int)
nsavl       1000

# Input fasta sequence file. (str)
fasta       1pgb.fasta

# Output file names (str)
logname     test.log
psfname     test.psf
dcdname     test.dcd
```

Once build the program with `make` command in the root folder, run with `./ez-mc config.cfg`
