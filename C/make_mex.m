% This will compile all mex files
mex mex/FastHadamard.c fht.c
mex mex/SparseFHT_mex.c sfht.c fht.c gf2.c
mex CFLAGS="\$CFLAGS -O2 -Wall" mex/HadamardBenchmark.c sfht.c fht.c gf2.c
