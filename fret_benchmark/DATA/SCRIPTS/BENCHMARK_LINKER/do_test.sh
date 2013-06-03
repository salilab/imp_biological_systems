## EDIT HERE
# IMP location 
export IMP=/netapp/sali/mbonomi/build/IMP-git/imp-fast/setup_environment.sh
# Root directory
export ROOT=/netapp/sali/mbonomi/working/FRET_BENCHMARK/

## DO NOT CHANGE BELOW
# Directory where tests are executed and results are stored
export RESULTDIR=${ROOT}/RESULTS_LINKER
# Data directory (with useful stuff)
export DATADIR=${ROOT}/DATA
# Gaussian Mixture Models for model of linker flexibility
export GMMDIR=${DATADIR}/GMM
# PDB location
export PDBDIR=${DATADIR}/PDBs
# SCRIPTs needed by the benchmark
export SCRIPTDIR=${DATADIR}/SCRIPTS/BENCHMARK_LINKER

# JOB ID is given in input
ID=$1

# create working dir
mkdir ${RESULTDIR}/testdir${ID}

# create command lines for test #ID from global list
sed -n ${ID}p GENERATE_COMMAND_LIST > ${RESULTDIR}/testdir${ID}/go.sh
sed -n ${ID}p SAMPLE_COMMAND_LIST  >> ${RESULTDIR}/testdir${ID}/go.sh                                   

# go into directory
cd ${RESULTDIR}/testdir${ID}

# run 
IMMPY=$IMP SCRIPTDIR=$SCRIPTDIR PDBDIR=$PDBDIR GMMDIR=$GMMDIR bash go.sh

# go back
cd $SCRIPTDIR
