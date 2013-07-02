#!/bin/sh

# Root directory
export ROOT=`cd ../../../ && pwd`

# Directory where tests are executed and results are stored
export RESULTDIR=${ROOT}/RESULTS_NOLINKER
# Data directory (with useful stuff)
export DATADIR=${ROOT}/DATA
# PDB location
export PDBDIR=${DATADIR}/PDBs
# SCRIPTs needed by the benchmark
export SCRIPTDIR=${DATADIR}/SCRIPTS/BENCHMARK_NOLINKER

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
SCRIPTDIR=$SCRIPTDIR PDBDIR=$PDBDIR bash go.sh

# go back
cd $SCRIPTDIR
