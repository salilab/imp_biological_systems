#! /usr/bin/env python                                                           
import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.isd
import math
import random
import argparse                                                                  
                                                                                 
parser = argparse.ArgumentParser(description='Gaussian Mixture Model of a GFP fluorophore density around a protein')
parser.add_argument('--pdbs',       dest="pdbs",         nargs="+", help="Pdb files")
parser.add_argument('--refpdbs',    dest="refpdbs",      nargs="+", help="Reference Pdb files")
parser.add_argument('--names',      dest="names",        nargs="+", help="Name of the proteins")
parser.add_argument('--pbl',        dest="pbl",          help="Photobleaching probability", default="1.0")
parser.add_argument('--R0',         dest="R0",           help="R0 value", default="49.0")
parser.add_argument('--Ida',        dest="Ida",          help="Ida value",    default="6.0")
parser.add_argument('--errIda',     dest="errIda",       help="Error on Ida", default="2.0")       
parser.add_argument('--norandom',   dest="norandom",action="store_true", default=False, help="No random initial conformation. Default False.")
parser.add_argument('--datanumber', dest="datanumber",   help="Fraction of total data.")
parser.add_argument('--indatafile', dest="indatafile",   help="Input data filename.")
parser.add_argument('--outdatafile', dest="outdatafile", help="Output data filename after random data selection.")
parser.add_argument('--fixonerigidbody', dest="fixonerigidbody",action="store_true", default=False, help="Fix one of the rigid bodies. Default False.")
parser.add_argument('--Box',        dest="Box",          help="Box size",    default="100.0")

result=parser.parse_args()

# check stuff
if(len(result.pdbs)!=len(result.names)):
  print "one name per pdb file is needed"
  exit()

# parameters
PBL_=float(result.pbl)
R0_=float(result.R0)
IDA_=float(result.Ida)
errIDA_=float(result.errIda)
NORANDOM_=result.norandom
DATAFILE_=result.indatafile
DATANUMBER_=float(result.datanumber)
BOX_=float(result.Box)
FIXONERIGIDBODY_=result.fixonerigidbody
# Sigma is sampled
SIGMA0_=0.002
SIGMA0MIN_=0.001
SIGMA0MAX_=0.01
# Mc parameters
MAXTRANS_=1.0
MAXROT_=0.3
TEMPMIN_=1.0
TEMPMAX_=5.0
NITER_=300000
NOPT_=100                                                                        

# create model
m=IMP.Model()                                                                    

# prepare restraint set
rst_dict={}

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
        xyzr=IMP.core.XYZR(atom)
        xyzr.set_radius(radius)                                   
        xyzr.set_coordinates_are_optimized(True)
        atoms.append(xyzr)                                                
    return atoms

def get_excluded_volume(atoms,kappa):
    lsc=IMP.container.ListSingletonContainer(atoms)
    evr=IMP.core.ExcludedVolumeRestraint(lsc,kappa)
    return evr 

def create_rigid_body(atoms):
    prb=IMP.Particle(m)                                                          
    rb=IMP.core.RigidBody.setup_particle(prb,atoms)                              
    rb.set_coordinates_are_optimized(True)
    return rb

def get_grid(gmin,gmax,ngrid,boundaries):
    grid=[]
    dx = ( gmax - gmin ) / float(ngrid)
    for i in range(0,ngrid+1):
        if(not boundaries and i==0): continue
        if(not boundaries and i==ngrid): continue
        grid.append( gmin + float(i) * dx )
    return grid

def get_log_grid(gmin,gmax,ngrid):
   grid=[]
   for i in range(0,ngrid+1):
       grid.append( gmin*math.exp(float(i)/ngrid*math.log(gmax/gmin)) )
   return grid

def get_drmsd_pairs(prot):
    pairs=[]
    for i in range(0,len(prot)-1):
        for j in range(i+1,len(prot)):
            pairs.append((prot[i],prot[j]))
    return pairs

def get_drmsd_reference(pairs):
    ref_dist=[]
    for p in pairs:
        dist=IMP.core.get_distance(IMP.core.XYZ(p[0]),IMP.core.XYZ(p[1]))
        ref_dist.append(dist)
    return ref_dist 

def get_drmsd(pairs, ref_dist):
    drmsd=0.
    for i,p in enumerate(pairs): 
        dist=IMP.core.get_distance(IMP.core.XYZ(p[0]),IMP.core.XYZ(p[1]))
        drmsd+=(dist-ref_dist[i])**2
    return math.sqrt(drmsd/float(len(pairs)))

def get_rb_movers(rblist,tmax,rmax,id_fix):
    mvs=[]
    for i,rb in enumerate(rblist):
        if (FIXONERIGIDBODY_ and i==id_fix): continue
        mv= IMP.core.RigidBodyMover(rb, tmax, rmax)
        mvs.append(mv)
    return mvs

def shuffle_configuration(rb,c,bbl):
    hbbl = bbl / 2.0
    ub = IMP.algebra.Vector3D(c[0]-hbbl, c[1]-hbbl, c[2]-hbbl)
    lb = IMP.algebra.Vector3D(c[0]+hbbl, c[1]+hbbl, c[2]+hbbl)
    bb = IMP.algebra.BoundingBox3D(ub, lb)
    translation    = IMP.algebra.get_random_vector_in(bb)
    rotation       = IMP.algebra.get_random_rotation_3d()
    transformation = IMP.algebra.Transformation3D(rotation, translation)
    rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))

def recenter_rb(rb):
    rotation = IMP.algebra.get_random_rotation_3d()
    translation = IMP.algebra.Vector3D(0.,0.,0.)
    transformation = IMP.algebra.Transformation3D(rotation, translation) 
    rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))
 
def get_external_barrier(rblist,bbl,c,kappa):
    hbbl = bbl / 2.0
    ub= IMP.core.HarmonicUpperBound(hbbl, kappa)
    ss= IMP.core.DistanceToSingletonScore(ub, c)
    lsc= IMP.container.ListSingletonContainer(rblist)
    rst= IMP.container.SingletonsRestraint(ss, lsc)
    return rst
 
def temp_simulated_annealing(istep,ncold,nhot):
     if istep%(ncold+nhot)< ncold:
        value=0.0
     else:
        value=1.0
     temp=TEMPMIN_+(TEMPMAX_-TEMPMIN_)*value
     return temp

def get_coordinates(ps):
    coords=[]
    for p in ps:
        coords.append(IMP.core.XYZ(p).get_coordinates())
    return coords

def get_placement_distance_angle(placement_atoms,placement_ref_atoms):

    # dictionaries
    coords_dict={}; ref_coords_dict={};
    # and lists
    coords_list=[]; ref_coords_list=[]

    # get coordinates 
    for key in placement_atoms:
        # get particles
        at     = placement_atoms[key]
        ref_at = placement_ref_atoms[key]
        # fill dictionaries
        coords_dict[key]     = get_coordinates(at)
        ref_coords_dict[key] = get_coordinates(ref_at)
        # extend lists
        coords_list     += coords_dict[key]
        ref_coords_list += ref_coords_dict[key]

    # find optimal transformation
    optimal = IMP.algebra.get_transformation_aligning_first_to_second(coords_list, ref_coords_list) 

    # initialize counters
    distance = 0.
    angle = 0.
    counter = 0.

    for key in placement_atoms:
        # increase counter
        counter += 1.
        # get aligned coordinates
        coords=[];
        for c in coords_dict[key]:
            coords.append(optimal.get_transformed(c))
        # get centroids
        model_centroid  = IMP.algebra.get_centroid(coords)
        native_centroid = IMP.algebra.get_centroid(ref_coords_dict[key])
        # and distance between them
        translation_vector = native_centroid - model_centroid
        distance += translation_vector.get_magnitude()
        # get angle
        TT = IMP.algebra.get_transformation_aligning_first_to_second(coords,ref_coords_dict[key])
        P  = IMP.algebra.get_axis_and_angle( TT.get_rotation() )
        angle += P.second

    # calculate average
    distance /= counter
    angle /= counter

    return distance,angle

# read pdb files, coarse grained with CA
molecules={}
for i in range(0,len(result.pdbs)):
    prot=load_pdb(result.names[i],result.pdbs[i])
    molecules[result.names[i]]=prot

# read reference positions
ref_molecules={}
for i in range(0,len(result.refpdbs)):                                              
    prot=load_pdb(result.names[i],result.refpdbs[i])                                
    ref_molecules[result.names[i]]=prot

# make rigid bodies
rbs={}; all_atoms=[]; rb_list=[]
# global DRMS lists
drms_atoms=[]; ref_drms_atoms=[];
# positional DRMS lists
drms_ter_atoms=[]; ref_drms_ter_atoms=[];
n_atoms=[]
# dictionary of termini particles for FRET
ps={}
# dictionay for placement scores
placement_atoms={}; placement_ref_atoms={}

for key in molecules:
    # get list of atoms                                                              
    prot_atoms=get_atoms(molecules[key])
    ref_prot_atoms=get_atoms(ref_molecules[key])
    # collect list of all atoms for evr
    all_atoms+=prot_atoms
    # initialize lists for placement scores
    placement_atoms[key]=[]
    placement_ref_atoms[key]=[]
    # collect 10 atoms for global DRMS
    for i in range(0,len(prot_atoms),len(prot_atoms)/10): 
        drms_atoms.append(prot_atoms[i])
        ref_drms_atoms.append(ref_prot_atoms[i])
        # and also for placement scores
        placement_atoms[key].append(prot_atoms[i])
        placement_ref_atoms[key].append(ref_prot_atoms[i])
    # collect termini for positional DRMS
    drms_ter_atoms.append(prot_atoms[0])
    drms_ter_atoms.append(prot_atoms[-1])
    ref_drms_ter_atoms.append(ref_prot_atoms[0])
    ref_drms_ter_atoms.append(ref_prot_atoms[-1])
    # make it a rigid body
    rb=create_rigid_body(prot_atoms)
    rbs[key+"-N"]=rb; rbs[key+"-C"]=rb
    rb_list.append(rb)
    # count number of atoms
    n_atoms.append(len(prot_atoms))
    # store termini
    ps[key+"-N"]=[prot_atoms[0].get_particle()] 
    ps[key+"-C"]=[prot_atoms[-1].get_particle()] 

# get the pairs of atoms for global and terminal drms calculation
drms_pairs=    get_drmsd_pairs(drms_atoms)
drms_ter_pairs=get_drmsd_pairs(drms_ter_atoms)
# get refence pairs and distances for drms
ref_drms_pairs=    get_drmsd_pairs(ref_drms_atoms)
ref_drms_ter_pairs=get_drmsd_pairs(ref_drms_ter_atoms)
ref_dist=    get_drmsd_reference(ref_drms_pairs)
ref_ter_dist=get_drmsd_reference(ref_drms_ter_pairs)

# excluded volume
rst_dict["Excluded_Volume"]=get_excluded_volume(all_atoms,1.0)

# read FRET data
for line in open(DATAFILE_, "r").readlines():
    data_dict=eval(line)
# list of all data
fret_list=[]
for key in data_dict:
    if "Data_Point_Average" in key:
       fret_list.append(data_dict[key])
# pick random data from fret_list
ndata=int((DATANUMBER_)*len(fret_list))
selected_fret_list=random.sample(fret_list,ndata)
# save all this info to file
output={}
output["NoRandom"]=NORANDOM_
output["FixOneRigidBody"]=FIXONERIGIDBODY_
output["DataFile"]=DATAFILE_
output["DataNumber"]=DATANUMBER_
output["Pbl"]=PBL_
output["R0"]=R0_
output["Ida"]=IDA_
output["errIDA"]=errIDA_
for f in selected_fret_list:
    output["Data_Point_Average|"+f[0][0]+"_"+f[0][1]]=f
log_output=open(result.outdatafile,"w")
log_output.write(str(output))

# create particles
# kda particle                                                                   
kda=IMP.isd.Scale.setup_particle(IMP.Particle(m),15.0)                           
kda.set_lower(1.0)                                                          
kda.set_upper(30.0)                                                          
kda.set_is_optimized(kda.get_nuisance_key(),True)                                
# maximum step per mover                                                         
MAX_KDA_MOVE_=0.3
# Ida particle                                                                   
Ida=IMP.isd.Scale.setup_particle(IMP.Particle(m),IDA_)                          
Ida.set_lower(1.0)
Ida.set_upper(10.0)
Ida.set_is_optimized(Ida.get_nuisance_key(),True)                                
# maximum step per mover                                                         
MAX_IDA_MOVE_=0.3                                                               
# add Gaussian restraint on Ida                                                  
rst_dict["Ida_Score"] = IMP.isd.GaussianRestraint(Ida,IDA_,errIDA_)
# Sigma0 particle                                                                   
Sigma0=IMP.isd.Scale.setup_particle(IMP.Particle(m),SIGMA0_)                          
Sigma0.set_lower(SIGMA0MIN_)                                                               
Sigma0.set_upper(SIGMA0MAX_)                                                  
Sigma0.set_is_optimized(Sigma0.get_nuisance_key(),True)                                
# maximum step per mover                                                         
MAX_SIGMA0_MOVE_=0.001
# Pbleaching particle
Pbl=IMP.isd.Scale.setup_particle(IMP.Particle(m),PBL_)
Pbl.set_lower(PBL_)                                                               
Pbl.set_upper(PBL_)                                                             
Pbl.set_is_optimized(Ida.get_nuisance_key(),False)
# R0 particle
R0=IMP.isd.Scale.setup_particle(IMP.Particle(m),R0_)                          
R0.set_lower(R0_)                                                               
R0.set_upper(R0_)                                                             
R0.set_is_optimized(Ida.get_nuisance_key(),False)

# create restraints
fret_set=IMP.RestraintSet("Fret")
fret_dict={}
for f in selected_fret_list:
    d=f[0][0]; a=f[0][1]; fexp=f[1]
    rst=IMP.isd.FretRestraint(ps[d],ps[a],kda,Ida,R0,Sigma0,Pbl,fexp)
    fret_dict[d+"_"+a]=rst
    fret_set.add_restraint(rst) 

rst_dict["Fret_Score"]=fret_set

# id of the biggest rigid body
id_fix=n_atoms.index(max(n_atoms))

if not NORANDOM_:
# Randomize initial configuration inside a bounding box
# centered on the biggest rigid body, whose
# position is not randomized
   for i,rb in enumerate(rb_list):
       if(i==id_fix): continue
       shuffle_configuration(rb,rb_list[id_fix].get_coordinates(),BOX_)

# set barrier centered on the position of the biggest rigid body
rst_dict["External_barrier"]=get_external_barrier(rb_list,BOX_,rb_list[id_fix].get_coordinates(),10.)

# add restraints
for rst in rst_dict:
    m.add_restraint(rst_dict[rst])

# Movers
mvs=[]                                                                           
mvs.append(IMP.core.NormalMover([kda],IMP.FloatKeys([IMP.FloatKey("nuisance")]),MAX_KDA_MOVE_))
mvs.append(IMP.core.NormalMover([Ida],IMP.FloatKeys([IMP.FloatKey("nuisance")]),MAX_IDA_MOVE_))
mvs.append(IMP.core.NormalMover([Sigma0],IMP.FloatKeys([IMP.FloatKey("nuisance")]),MAX_SIGMA0_MOVE_))
# rigid body movers
mvs+=get_rb_movers(rb_list,MAXTRANS_,MAXROT_,id_fix)
# SerialMover                                                                    
smv=IMP.core.SerialMover(mvs)

# setting up MonteCarlo                                                          
mc=IMP.core.MonteCarlo(m)                                                        
mc.set_return_best(False)                                                        
mc.set_kt(1.0)                                                                   
mc.add_mover(smv)

# reset lowest score
score_min = 100000000.0

# Hierarchy for pdb
h=IMP.atom.Hierarchy(IMP.Particle(m))
# add chain from molecules (in the order provided in input)
for i in range(0,len(result.names)):
    h.add_child(molecules[result.names[i]].get_children()[0])

# Sampling
for istep in range(0,NITER_):                                                    

    # add annealing schedule
    temp=temp_simulated_annealing(istep,1000,300)
    mc.set_kt(temp)

    # Monte Carlo
    mc.optimize(NOPT_)

    # get score 
    score = m.evaluate(False)

    # printout stats only for best posterior model (save space)
    # this can be easily changed
    if(score<score_min):
        # new minimum
        score_min = score
        # prepare printout
        output={}
        output["Total_Score"]=score
        output["Step_Number"]=istep                                              
        output["Temperature"]=mc.get_kt()                                        
        output["Acceptance"]=mc.get_number_of_forward_steps()/float(NOPT_)       
        for key in rst_dict: 
           output[key]=rst_dict[key].evaluate(False)                        
        output["Ida"]=Ida.get_scale()
        output["kda"]=kda.get_scale()
        output["Sigma0"]=Sigma0.get_scale()
        for key in fret_dict:
            output["Experimental-"+key]=fret_dict[key].get_experimental_value()
            output["Model-"+key]=fret_dict[key].get_model_fretr()
            output["StdError-"+key]=fret_dict[key].get_standard_error()

        output["DRMS"]=get_drmsd(drms_pairs,ref_dist)
        output["DRMS-ter"]=get_drmsd(drms_ter_pairs,ref_ter_dist)

        pl_dist,pl_angle = get_placement_distance_angle(placement_atoms,placement_ref_atoms)

        output["Placement_distance"]=pl_dist
        output["Placement_angle"]=pl_angle

        # open log file
        log=open("best_posterior.dat","w")
        # write dictionary
        log.write("%s \n" % output)
        # close log file
        log.close()

        # print all information to pdb 
        IMP.atom.write_pdb(h,"best_posterior.pdb") 
