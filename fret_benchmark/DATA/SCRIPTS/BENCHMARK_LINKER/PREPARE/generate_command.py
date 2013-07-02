#!/bin/env python

import sys

sigma0=[0.001,0.01]
ndata=[0.5,1.0]
ntrials=20

# read pdbs from file
complexes3=[]
for lines in open(sys.argv[1],"r").readlines():
    riga=lines.strip().split()
    complexes3.append((riga[0],riga[1],riga[2]))

# complexes with 3 SUBUNITS
for c in complexes3:
    for s in sigma0:
        for n in ndata:
            for i in range(0,ntrials):
                print "${SCRIPTDIR}/generate_DATA.py --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --name %s --name %s --name %s --ntrials 50 --R0 49.0 --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --Ida 6.0 --Kda 16.0 --Sigma0 %lf --sensitivity 0.0001" % (c[0],c[1],c[2],c[0],c[1],c[2],c[0],c[1],c[2],s)
                print "${SCRIPTDIR}/sample.py --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --name %s --name %s --name %s --refpdb ${PDBDIR}/%s.pdb --refpdb ${PDBDIR}/%s.pdb --refpdb ${PDBDIR}/%s.pdb --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --datanumber %lf --indatafile data.out --outdatafile selected-data.out --fixonerigidbody --Box 100" % (c[0],c[1],c[2],c[0],c[1],c[2],c[0],c[1],c[2],c[0],c[1],c[2],n)

# read pdbs from file
complexes4=[]
for lines in open(sys.argv[2],"r").readlines():
    riga=lines.strip().split()
    complexes4.append((riga[0],riga[1],riga[2],riga[3]))

# complexes with 4 SUBUNITS
for c in complexes4:
    for s in sigma0:
        for n in ndata:
            for i in range(0,ntrials):
                print "${SCRIPTDIR}/generate_DATA.py --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --name %s --name %s --name %s --name %s --ntrials 50 --R0 49.0 --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --Ida 6.0 --Kda 16.0 --Sigma0 %lf --sensitivity 0.0001" % (c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],s)
                print "${SCRIPTDIR}/sample.py --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --pdb ${PDBDIR}/%s.pdb --name %s --name %s --name %s %s --refpdb ${PDBDIR}/%s.pdb --refpdb ${PDBDIR}/%s.pdb --refpdb ${PDBDIR}/%s.pdb --refpdb ${PDBDIR}/%s.pdb --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --GMM ${GMMDIR}/%s.out --datanumber %lf --indatafile data.out --outdatafile selected-data.out --fixonerigidbody --Box 100" % (c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],n)
