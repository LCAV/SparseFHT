#!/bin/bash

# This will run all the scripts to reproduce the figures of the paper
# A fast Hadamard transform for signals sparse in the tranform domain
# by Robin Scheibler, Saeid Haghighatshoar, and Martin Vetterli
# 
# This bash file was written by Robin Scheibler, October 2016
# License: MIT

# Config
########

# Enable test mode. This will run a single loop of each
# script. This can be used to test if the configuration
# is correct
ENABLE_TEST=0

# Enable serial mode. This deactivate the use of parallel
# workers. The code runs in a straight loop.
ENABLE_SERIAL=0

# This sets the number of nodes in the ipyparallel cluster
# If no cluster is used, this can be set to zero to run
# in serial mode (super slow though...)
# This number can be set to the number of available threads of
# your CPU minus 1. Usually, the number of threads is twice
# the number of cores.
PARALLEL_WORKERS=4

# Show help function
####################
function show_help {
  echo "$1 [OPTS]"
  echo "Options:"
  echo "  -t    Runs a single loop only for test purpose"
  echo "  -s    Runs all the code in a simple for loop. No parallelism"
  echo "  -n x  Runs the loops in parallel using x workers. This option is ignored if -s is used"
}

# Process arguments
###################

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?tsn:" opt; do
  case "$opt" in
    h|\?)
      show_help $0
      exit 0
      ;;
    t)  ENABLE_TEST=1
      ;;
    s)  ENABLE_SERIAL=1
      ;;  
    n)  PARALLEL_WORKERS=$OPTARG
      ;;
  esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# Process SERIAL flag
if [ $ENABLE_SERIAL -eq 1 ]; then
  PARALLEL_WORKERS=0
fi

# Run all the scripts
#####################

# Prepare parallel processing
if [ $PARALLEL_WORKERS -gt 0 ]; then
  echo "Starting ${PARALLEL_WORKERS} ipyparallel workers."
  ipcluster start -n ${PARALLEL_WORKERS} --daemonize
  echo "Wait for 30 seconds to give time to engines to start..."
  sleep 30
  SERIAL_FLAG=
else
  SERIAL_FLAG=-s
fi

# Process test flag
if [ $ENABLE_TEST -eq 1 ]; then
  TEST_FLAG=-t
else
  TEST_FLAG=
fi

# Make some folders
mkdir -p figures
mkdir -p data

FLAGS="${SERIAL_FLAG} ${TEST_FLAG}"
echo "Running with flags ${FLAGS}"

# Run all the scripts and get the output data file name
echo 'Running Monte-Carlo simulation for the probability of error...'
FILE1=`python error_sim.py ${FLAGS} | grep 'Saved data to file:' | awk '{ print $5 }'`
echo 'Running Monte-Carlo for the less sparse regime...'
FILE2=`python less_sparse_sim.py ${FLAGS} | grep 'Saved data to file:' | awk '{ print $5 }'`
echo 'Running the runtime benchmark...'
FILE3=`python timing_sim.py ${FLAGS} | grep 'Saved data to file:' | awk '{ print $5 }'`

echo "All processing done! The data was saved in files:"
echo "  ${FILE1}"
echo "  ${FILE2}"
echo "  ${FILE3}"

# Now produce all the plots
echo 'Creating all figures...'
python error_plot.py -f $FILE1
python less_sparse_plot.py -f $FILE2
python timing_plot.py -f $FILE3
python density_evolution.py

if [ $PARALLEL_WORKERS -gt 0 ]; then
  echo 'Stopping parallel processing now.'
  ipcluster stop
fi

echo 'All done. See you soon...'

