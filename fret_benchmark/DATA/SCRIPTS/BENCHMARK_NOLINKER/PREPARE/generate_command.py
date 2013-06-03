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
                print "$IMMPY python ${SCRIPTDIR}/generate_DATA.py --pdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --names %s %s %s --ntrials 50 --R0 49.0 --Ida 6.0 --Kda 7.5 --Sigma0 %lf --sensitivity 0.0001" % (c[0],c[1],c[2],c[0],c[1],c[2],s)
                print "$IMMPY python ${SCRIPTDIR}/sample.py --pdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --names %s %s %s --refpdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --R0 49.0 --datanumber %lf --indatafile data.out --outdatafile selected-data.out --fixonerigidbody --Box 100" % (c[0],c[1],c[2],c[0],c[1],c[2],c[0],c[1],c[2],n)

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
                print "$IMMPY python ${SCRIPTDIR}/generate_DATA.py --pdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --names %s %s %s %s --ntrials 50 --R0 49.0 --Ida 6.0 --Kda 7.5 --Sigma0 %lf --sensitivity 0.0001" % (c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],s)
                print "$IMMPY python ${SCRIPTDIR}/sample.py --pdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --names %s %s %s %s --refpdbs ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb ${PDBDIR}/%s.pdb --R0 49.0 --datanumber %lf --indatafile data.out --outdatafile selected-data.out --fixonerigidbody --Box 100" % (c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],c[0],c[1],c[2],c[3],n)
