#!/bin/sh
#$ -S /bin/sh
#$ -N FRET-BENCHMARK-GMM 
#$ -o /netapp/sali/mbonomi/working/FRET_BENCHMARK/OUT_ERR 
#$ -e /netapp/sali/mbonomi/working/FRET_BENCHMARK/OUT_ERR 
#$ -r n
#$ -j n
#$ -l netapp=1G
#$ -l h_rt=24:00:00
#$ -l mem_free=2G
#$ -q lab.q
#$ -l xe5430=true
#$ -l arch=linux-x64 
#$ -R yes
#$ -V
#$ -cwd
#$ -t 1-53 

# IMP location
export IMP=/netapp/sali/mbonomi/build/IMP-git/imp-fast/setup_environment.sh
# Root directory for the benchmark
export ROOT=/netapp/sali/mbonomi/working/FRET_BENCHMARK/
# Data directory
export DATADIR=${ROOT}/DATA
# Gaussian Mixture models for linker flexibility model
export GMMDIR=${DATADIR}/GMM
# PDB Location of the PDBs
export PDBDIR=${DATADIR}/PDBs
# SCRIPTs 
export SCRIPTDIR=${DATADIR}/SCRIPTS/BENCHMARK_LINKER/FIT_GMM

# JOB ID
ID=$SGE_TASK_ID

# parameters
nprobes=5000                                                    
rprobe=10.0
ncomp=10

# create command lines
PREFIX=`sed -n ${ID}p ${PDBDIR}/pdblist | cut -d. -f1` 

# get N terminus resid and chain                                                 
N=`grep  ATOM ${PDBDIR}/${PREFIX}.pdb | head -n 1  | awk '{print $6}'`
CN=`grep ATOM ${PDBDIR}/${PREFIX}.pdb | head -n 1  | awk '{print $5}'`
# get C terminus resid and chain                                                 
C=`grep  ATOM ${PDBDIR}/${PREFIX}.pdb | tail -n 1  | awk '{print $6}'`
CC=`grep ATOM ${PDBDIR}/${PREFIX}.pdb | tail -n 1  | awk '{print $5}'`

# go to GMMDIR
cd $GMMDIR

#run
$IMP python ${SCRIPTDIR}/fit_GMM.py --pdb ${PDBDIR}/${PREFIX}.pdb  --nprobes $nprobes --rprobe $rprobe --ncomp $ncomp --histoN ${SCRIPTDIR}/GFP_to_N.histo --histoC ${SCRIPTDIR}/GFP_to_C.histo --out ${PREFIX}.out  --residN $N --residC $C --chainN $CN --chainC $CC
