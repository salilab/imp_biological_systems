#! /usr/bin/env python                                                           
import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf                                                                   
import RMF
import random
import argparse                                                                  
import numpy as np
import sklearn
import sklearn.mixture
from sklearn.mixture import GMM                                                  
                                                                                 
parser = argparse.ArgumentParser(description='Gaussian Mixture Model of a GFP fluorophore density around a protein')
parser.add_argument('--pdb',       dest="pdbfile",     help="Pdb file of the protein")
parser.add_argument('--nprobes',   dest="nprobes",     help="Number of probe points")
parser.add_argument('--rprobe',    dest="rprobe",      help="GFP probe radius")
parser.add_argument('--histoN',    dest="histoNfile",  help="Radial distribution of a GFP fluorophore attached to the N terminus")
parser.add_argument('--histoC',    dest="histoCfile",  help="Radial distribution of a GFP fluorophore attached to the C terminus")
parser.add_argument('--ncomp',     dest="ncomp",       help="Maximum number of components for GMM")
parser.add_argument('--covtype',   dest="covtype",     default="spherical", help="Covariance type: spherical, diagonal, full")
parser.add_argument('--out',       dest="outfile",     default="data.out", help="Output data file")
parser.add_argument('--writermf',  dest="writermf",    action="store_true", default=False, help="Write rmf with probes. Default False.")
parser.add_argument('--rmf',       dest="rmffile",     default="probes.rmf", help="Output RMF")
parser.add_argument('--residN',    dest="residN",      help="Residue number for N terminus")
parser.add_argument('--chainN',    dest="chainN",      default="", help="Chain id for N terminus")
parser.add_argument('--residC',    dest="residC",      help="Residue number for C terminus")
parser.add_argument('--chainC',    dest="chainC",      default="", help="Chain id for C terminus")

result=parser.parse_args()

NPROBES_=int(result.nprobes)
RPROBE_=float(result.rprobe)
PDBFILE_=result.pdbfile
HISTONFILE_=result.histoNfile
HISTOCFILE_=result.histoCfile                                                     
NCOMP_=int(result.ncomp)
COVTYPE_=result.covtype
WRITERMF_=result.writermf
RMFFILE_=result.rmffile
OUTFILE_=result.outfile

# create model
m=IMP.Model()                                                                    

def read_histo_file(filename):
    x=[]; y=[]
    for line in open(filename).readlines():                                       
        riga=(line.strip()).split()                                                    
        x.append(float(riga[0]))                                                    
        y.append(float(riga[1]))                                                    
    return x,y
                                                                                 
def load_pdb(name,filename):                                                     
    prot=IMP.atom.read_pdb(filename,m,IMP.atom.CAlphaPDBSelector())              
    prot.set_name(name)                                                          
    return prot                                                                  

def get_atoms(prot):
    atoms=[]                                                              
    for atom in IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE):                  
        residue=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent())               
        restype=residue.get_residue_type()                                       
        vol=IMP.atom.get_volume_from_residue_type(restype)                       
        radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)                   
        IMP.core.XYZR(atom).set_radius(radius)                                   
        atoms.append(atom)                                                
    return atoms

def get_excluded_volume(atoms,kappa):
    lsc=IMP.container.ListSingletonContainer(atoms)
    evr=IMP.core.ExcludedVolumeRestraint(lsc,kappa)
    return evr 

def create_rigid_body(atoms):
    prb=IMP.Particle(m)                                                          
    rb=IMP.core.RigidBody.setup_particle(prb,atoms)                              
    return rb                                                                    

def get_cumulative(x,y):
    tot=0.
    cumul=[]
    cumul.append(0.)
    for i in range(1,len(x)):
        dx=x[i]-x[i-1]
        tot+=(y[i]+y[i-1])/2.0*dx
        cumul.append(tot) 
    # normalize
    for i in range(0,len(cumul)):
        cumul[i]/=tot
    return cumul

def get_random_radius(x,y_cumul):
    a=random.random() 
    closest=min([(abs(a-z), z) for z in y_cumul])[1]                                       
    i=y_cumul.index(closest)                                                 
    return x[i]

def get_probes(prot_atoms,pterminus,x,y_cumul,name):
    # create hierarchy
    probes=IMP.atom.Hierarchy(IMP.Particle(m))                                      
    probes.set_name(name)
    # get position of the terminus
    center=IMP.core.XYZ(pterminus).get_coordinates() 
    counter=0
    while (counter < NPROBES_):
       # new particle
       p=IMP.Particle(m)
       # decorate as atom 
       atom=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)
       # get radius from distribution
       radius=get_random_radius(x,y_cumul)
       position=IMP.algebra.get_random_vector_on(IMP.algebra.Sphere3D(center,radius))
       # and decorate as XYZR
       IMP.core.XYZR.setup_particle(p,IMP.algebra.Sphere3D(position,RPROBE_)) 
       # check for clashes
       evr=get_excluded_volume(prot_atoms+[atom],1.0)
       if(evr.evaluate(False)<=0.):
         probes.add_child(atom)
         counter+=1
    return probes


# read pdb file, coarse grained with CA
prot=load_pdb("protein",PDBFILE_)

# get list of atoms                                                              
prot_atoms=get_atoms(prot)

# make it a rigid body
prot_rb=create_rigid_body(prot_atoms)

# read histogram of fluorophore-Nterminus distance
(xN,yN)=read_histo_file(HISTONFILE_)
# and get cumulative distribution
yN_cumul=get_cumulative(xN,yN)

# read histogram of fluorophore-Cterminus distance                               
(xC,yC)=read_histo_file(HISTOCFILE_)
# and get cumulative distribution
yC_cumul=get_cumulative(xC,yC)

p_ter={}
# get N terminus of the protein                                                  
s=IMP.atom.Selection(prot, chain=result.chainN, residue_index=int(result.residN))
p_ter["N"]=s.get_selected_particles()[0]
# get C terminus of the protein
s=IMP.atom.Selection(prot, chain=result.chainC, residue_index=int(result.residC))     
p_ter["C"]=s.get_selected_particles()[0]

# hierarchy for probes at N terminus
probesN=get_probes(prot_atoms,p_ter["N"],xN,yN_cumul,"ProbesN")

# hierarchy for probes at C terminus                                             
probesC=get_probes(prot_atoms,p_ter["C"],xC,yC_cumul,"ProbesC")

# print particles to rmf
if(WRITERMF_):
  # open rmf for writing                                                           
  rh = RMF.create_rmf_file(RMFFILE_)                             
  # add coordinates                                                                
  IMP.rmf.add_hierarchy(rh, prot)                   
  # add Nprobe particles                                                              
  IMP.rmf.add_hierarchy(rh,probesN)
  # add Cprobe particles                                                              
  IMP.rmf.add_hierarchy(rh,probesC)

output={}
# add general data
output["PDB"]=PDBFILE_
output["Nprobes"]=NPROBES_
output["Rprobe"]=RPROBE_
output["Covtype"]=COVTYPE_

# GMM fit
for h in [(probesN,"N"), (probesC,"C")]:
    # create numpy arrays
    atoms=[]                                                                        
    for a in IMP.atom.get_leaves(h[0]):                                           
        atoms.append(IMP.core.XYZR(a).get_coordinates())                            
    npa=np.array(atoms)

    bic=[]; classifiers=[]
    for i in range(NCOMP_,NCOMP_+1): 
        classifier=GMM(n_components=i, covariance_type=COVTYPE_, n_iter=1000)
        classifier.fit(npa)
        classifiers.append(classifier)
        bic.append(classifier.bic(npa))
    # find model that minimize BIC
    ncomp=bic.index(min(bic))

    output[h[1]+"-Bic"]=bic[ncomp]
    output[h[1]+"-Ncomp"]=ncomp+1
    # write list of sigmas
    sigmas=[]
    for c in classifiers[ncomp].covars_: 
        sigmas.append(c[0])       
    output[h[1]+"-Sigmas"]=sigmas
    # write list of weights
    weights=[] 
    for w in classifiers[ncomp].weights_: 
        weights.append(w)
    output[h[1]+"-Weights"]=weights
    # write centers in local coordinates
    rf=prot_rb.get_reference_frame()
    centers=[]
    for c in classifiers[ncomp].means_:
        c_loc=rf.get_local_coordinates(IMP.algebra.Vector3D(c[0],c[1],c[2])) 
        centers.append(c_loc)
    output[h[1]+"-Centers"]=centers
    # write coordinates of terminus in local coordinates
    ter=IMP.core.XYZ(p_ter[h[1]]).get_coordinates()
    output[h[1]+"-Terminus"]=rf.get_local_coordinates(ter)

    if(WRITERMF_):
      # write GMM center to rmf
      probes_GMM=IMP.atom.Hierarchy(IMP.Particle(m))
      probes_GMM.set_name(h[1]+"-GMM")
      for i in range(0,len(centers)):
          p=IMP.Particle(m)                                                      
          # decorate as atom                                                        
          atom=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)                    
          position=rf.get_global_coordinates(centers[i])
          r=sigmas[i]
          # and decorate as XYZR                                                    
          IMP.core.XYZR.setup_particle(p,IMP.algebra.Sphere3D(position,r)) 
          probes_GMM.add_child(atom)
      # add hierarchy to rmf                                                           
      IMP.rmf.add_hierarchy(rh,probes_GMM) 
      # sample from the GMM
      sample_np=classifiers[ncomp].sample(NPROBES_)
      # create hierarchy                                                           
      probes_sampled=IMP.atom.Hierarchy(IMP.Particle(m))                                   
      probes_sampled.set_name(h[1]+"-sampled")                                                        
      for s in sample_np:
          p=IMP.Particle(m)                                                         
          # decorate as atom                                                        
          atom=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)                       
          position=IMP.algebra.Vector3D(s[0],s[1],s[2])
          # and decorate as XYZR                                                    
          IMP.core.XYZR.setup_particle(p,IMP.algebra.Sphere3D(position,2.0))
          probes_sampled.add_child(atom)
      # add hierarchy to rmf
      IMP.rmf.add_hierarchy(rh,probes_sampled)

# save frame                                                                     
if(WRITERMF_): IMP.rmf.save_frame(rh,0)

# write dictionary to file
log=open(OUTFILE_,'w')                                                              
log.write("%s \n" % output)
