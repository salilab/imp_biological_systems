#!/usr/bin/env python
import IMP
import IMP.core
import IMP.em
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import sys


               


class run_class():


        #####################################################

    def create_rigid_body(self,rigid_body_list,prot):
        rigidbody_list=[]
        for element in rigid_body_list:
            atoms=[]
            for interval in element:
            #rinterval upper bound is incremented by one because
            #range() cuts the upper edge
                rinterval=(interval[0],interval[1]+1)
                if (interval[0]==-1 or interval[1]==-1):
                    s=IMP.atom.Selection(prot,chains=interval[2])
                else:
                    s=IMP.atom.Selection(prot,chains=interval[2],
                                         residue_indexes=range(rinterval))
                for p in s.get_selected_particles():
                    atoms.append(IMP.core.XYZR(p))
            prb=IMP.Particle(self.m)
            rb=IMP.core.RigidBody.setup_particle(prb,atoms)
            rigidbody_list.append(rb)
        return rigidbody_list

        #####################################################

    def link_domains(self,prot, resrangelist, kappa):
        rs = IMP.RestraintSet('linker')
        for pair in resrangelist:
            try:
                s0=IMP.atom.Selection(prot, chains=pair[2], residue_index=pair[0])
                p0=s0.get_selected_particles()[0]
            except:
                "error"
                continue
            try:
                s1=IMP.atom.Selection(prot, chains=pair[2], residue_index=pair[1])
                p1=s1.get_selected_particles()[0]
            except:
                "error"
                continue
            dist0=float(pair[1]-pair[0])*4.0

            h=IMP.core.HarmonicUpperBound(dist0, kappa)
            dps=IMP.core.DistancePairScore(h)
            pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p0,p1))
            rs.add_restraint(pr)

        self.m.add_restraint(rs)
        self.rest_set[('link',prot)]=rs

        #####################################################

    def add_excluded_volume(self,prot,kappa):
        rs = IMP.RestraintSet('excluded_volume')
        atoms=IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE)
        for atom in atoms:
            restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
            vol=IMP.atom.get_volume_from_residue_type(restype)
            radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
            IMP.core.XYZR(atom).set_radius(radius)
        lsa=IMP.container.ListSingletonContainer(self.m)
        lsa.add_particles(atoms)
        evr=IMP.core.ExcludedVolumeRestraint(lsa,kappa)
        rs.add_restraint(evr)
        self.m.add_restraint(rs)
        self.rest_set[('exvo',prot)]=rs

        #####################################################

    def add_symmetry_excluded_volume(self,prot_ref,prot_symm_list,kappa):
        rs = IMP.RestraintSet('symm_excluded_volume')
        
        atoms_ref=IMP.atom.get_by_type(prot_ref, IMP.atom.ATOM_TYPE)
        ls_ref=IMP.container.ListSingletonContainer(self.m)
        ls_ref.add_particles(atoms_ref)


        ls_symm=IMP.container.ListSingletonContainer(self.m)        
        for prot_symm in prot_symm_list:
           
           atoms_symm=IMP.atom.get_by_type(prot_symm, IMP.atom.ATOM_TYPE)
           ls_symm.add_particles(atoms_symm)
           for atom in atoms_symm:
             restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
             vol=IMP.atom.get_volume_from_residue_type(restype)
             radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
             IMP.core.XYZR(atom).set_radius(radius)
        
           cbpc=IMP.container.CloseBipartitePairContainer(ls_ref,ls_symm,1.0,10.0)
           ssps=IMP.core.SoftSpherePairScore(kappa)
           evr3=IMP.container.PairsRestraint(ssps,cbpc)
           rs.add_restraint(evr3)
        
        self.m.add_restraint(rs)
        self.rest_set[('exvo_symm',prot_ref)]=rs

        #####################################################

    def add_external_barrier(self,rad,prot):
        rs = IMP.RestraintSet('barrier')
        c3= IMP.algebra.Vector3D(0,0,0)
        ub3= IMP.core.HarmonicUpperBound(rad, 10.0)
        ss3= IMP.core.DistanceToSingletonScore(ub3, c3)
        lsc= IMP.container.ListSingletonContainer(self.m)
        IMP.atom.get_by_type

        lsc.add_particles(IMP.atom.get_leaves(prot))
        r3= IMP.container.SingletonsRestraint(ss3, lsc)
        rs.add_restraint(r3)
        self.m.add_restraint(rs)
        self.rest_set[('barr',prot)]=rs

        #####################################################

    def cross_link_ms_simple_symmetrized_create(self,prot1,prot2,restraints_file):
     """read crosslink restraints between two residue 
     of different chains from an external text file
     sintax: part_name_1 part_name_2 distance error
     example:     0 1 1.0 0.1"""
     rs=IMP.RestraintSet('xlms')
     self.pairs=[]
     
     crosslinker="BS3" 
    
     #this is an harmonic potential with mean 12 and sigma=5, therefore
     #k=1/sigma^2
     hf=IMP.core.TruncatedHarmonicBound(12.0,1.0/25.0,15.0,5)
     dps=IMP.core.DistancePairScore(hf)
  
     index=0

     addedd_pairs_list=[]     
     for line in open(restraints_file):

        #force_restraint=True makes all intra rigid body restraint to be accepted
        force_restraint=True
        tokens=line.split()

        #skip character
        if (tokens[0]=="#"): continue
        r1=int(tokens[0])
        c1=tokens[1]
        r2=int(tokens[2])
        c2=tokens[3]
        crosslinker=tokens[4]

        #two restraints with the same index will be ambiguous
        index=int(tokens[5])

        #force restraint even if it belong to the same rigid body, use it for ambiguous restraints
        if (tokens[len(tokens)-1]=="F"): force_restraint=True

        print "attempting to add restraint between residue %d of chain %s and residue %d of chain %s" % (r1,c1,r2,c2)

        p1s=[]
        p2s=[]
        
        
        if (c1==c2):
          
          #apply the cross-link to the main copy
          try:
            s1=IMP.atom.Selection(prot1, residue_index=r1, atom_type=IMP.atom.AT_CA)
            p1=(s1.get_selected_particles()[0])
          except:
            print "WARNING> residue %d of chain %s is not there" % (r1,c1)
            continue
          try:
            s2=IMP.atom.Selection(prot1, residue_index=r2, atom_type=IMP.atom.AT_CA)
            p2=(s2.get_selected_particles()[0])
          except:
            print "WARNING> residue %d of chain %s is not there" % (r2,c2)
            continue

        if (c1!=c2):
          
          #apply the cross-link to the main copy and the symmetric copy
          try:
            s1=IMP.atom.Selection(prot1, residue_index=r1, atom_type=IMP.atom.AT_CA)
            p1=(s1.get_selected_particles()[0])
          except:
            print "WARNING> residue %d of chain %s is not there" % (r1,c1)
            continue
          try:
            s2=IMP.atom.Selection(prot2, residue_index=r2, atom_type=IMP.atom.AT_CA)
            p2=(s2.get_selected_particles()[0])
          except:
            print "WARNING> residue %d of chain %s is not there" % (r2,c2)
            continue


        print "attempting to add restraint between residue %d of chain %s and residue %d of chain %s" % (r1,c1,r2,c2)
                              
        #check whether the atom pair belongs to the same rigid body          
        if(IMP.core.RigidMember.particle_is_instance(p1) and IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember(p1).get_rigid_body() == IMP.core.RigidMember(p2).get_rigid_body() and not force_restraint): 
               print "WARNING> residue %d of chain %s and residue %d of chain %s belong to the same rigid body" % (r1,c1,r2,c2)                
               continue

            #this list contain the list of simmetric pairs to avoid duplications
        if (p1,p2,crosslinker) in addedd_pairs_list: 
                print "WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)                
                continue
        if (p2,p1,crosslinker) in addedd_pairs_list: 
                print "WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)                
                continue
            
        print "added pair %d %s %d %s" % (r1,c1,r2,c2)   
        index+=1
        addedd_pairs_list.append((p1,p2,crosslinker))

        rs_name='restraint_'+str(index)

        ln=IMP.core.PairRestraint(dps,IMP.ParticlePair(p1,p2))       
        ln.set_weight(0.5)                    
        rs.add_restraint(ln)        
        
       
        self.pairs.append((p1,  p2,  crosslinker,  rs_name,  100,  100,  (r1,c1),  (r2,c2), crosslinker, ln))
    

     self.m.add_restraint(rs) 
     self.rest_set["xlms"]=rs    
     
    ##############################################

    def add_em_map(self,filename,domains,resname):
        for icopy in range(0,self.ncopies):
             
             particles=[]
             for interval in domains:
               #rinterval upper bound is incremented by one because
               #range() cuts the upper edge
                  rinterval=(interval[0],interval[1]+1)
                  if (interval[0]==-1 or interval[1]==-1):
                      s=IMP.atom.Selection(interval[3],chains=interval[2])
                  else:
                      s=IMP.atom.Selection(interval[3],chains=interval[2],
                                           residue_indexes=range(rinterval))
                  particles+=s.get_selected_particles()
             #read the map
             map = IMP.em.read_map(filename,IMP.em.MRCReaderWriter())                  
             

             #determine the density threshold  
             resnum=len(particles)
             mass = IMP.atom.get_mass_from_number_of_residues(resnum)
             density_threshold = IMP.em.get_threshold_for_approximate_mass(map, 1.25*mass)
             efr = IMP.em.EnvelopeFitRestraint(particles, map, 11.0, 5.0)
             efr.set_weight(6.)    
             self.m.add_restraint(efr)
             self.rest_set[(resname,self.prot[icopy])]=efr

        #####################################################
             
    def add_symmetry_restraint(self,(rigid_bodies_ref,rigid_bodies_copy),transformation):
         
        sm=IMP.core.TransformationSymmetry(transformation)
        lc=IMP.container.ListSingletonContainer(self.m)        
        for i in range(len(rigid_bodies_ref)):
           IMP.core.Reference.setup_particle(rigid_bodies_copy[i],rigid_bodies_ref[i])
           lc.add_particle(rigid_bodies_copy[i])
        c=IMP.container.SingletonsConstraint(sm,None,lc)
        self.m.add_score_state(c)

        #####################################################

    def add_template_restraint(self,ps1,ps2,label):
       rset=IMP.RestraintSet('template_restraint')   
       for p1 in  ps1:
           for p2 in ps2:
               #check that the two particles are not in the same rigid body
               if(IMP.core.RigidMember.particle_is_instance(p1) and IMP.core.RigidMember.particle_is_instance(p2) and
               IMP.core.RigidMember(p1).get_rigid_body() == IMP.core.RigidMember(p2).get_rigid_body()): continue
               d0=IMP.core.XYZ(p1)
               d1=IMP.core.XYZ(p2)
               dist=IMP.core.get_distance(d0,d1)             
               if dist <= 6.5:
                  hf=IMP.core.Harmonic(dist,1.0)
                  dps=IMP.core.DistancePairScore(hf)                  
                  pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p1,p2))
                  rset.add_restraint(pr)

       self.m.add_restraint(rset)
       self.rest_set[('Template_Restraint_'+label)]=rset

    
    #####################################################

    def setup_rigid_body_MonteCarlo(self,mc_kt,mc_dx=0.3,mc_dang=0.1,bm_dr=0.5,set_return_best=False):
        mc=IMP.core.MonteCarlo(self.m)
        mc.set_return_best(set_return_best)
        mc.set_kt(mc_kt)
        #create a list of rigid body movers for a serial mover

        for icopy in range(0,self.ncopies):
            mvs=[]
            for rb in self.rigidbody_list[icopy]:
                mvs.append(IMP.core.RigidBodyMover(rb, mc_dx, mc_dang))

            for ps in IMP.atom.get_leaves(self.prot[icopy]):

                if (not IMP.core.RigidMember.particle_is_instance(ps)):
                    mvs.append(IMP.core.BallMover([ps],bm_dr))

            mc.add_mover(IMP.core.SerialMover(mvs))

        return mc

    ##############################################
    
    def temp_simulated_annealing(self,istep,ncold,nhot):
        if istep%(ncold+nhot)< ncold:
            value=0.0
        else:
            value=1.0

        temp=self.tempmin+(self.tempmax-self.tempmin)*value
        return temp
 


        #####################################################

    def init_stats(self):
        statfile = 'output/models_stats.out'
        flstat=open(statfile,'w')
        flstat.close()

        #####################################################

    def write_stats(self):
        statfile='output/models_stats.out'
        flstat=open(statfile,'a')

        output={}
        output["Nframe"]=self.nframe
        output["Total_Score"]=self.m.evaluate(False)

        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            ln=self.pairs[i][9]
            resid1=self.pairs[i][6][0]
            chain1=self.pairs[i][6][1]
            resid2=self.pairs[i][7][0]
            chain2=self.pairs[i][7][1]
            
            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["xlms_"+label]=ln.evaluate(False)


            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["Distance_"+label]=IMP.core.get_distance(d0,d1)

            
        for i  in range(len(self.prot)):

            output["Excluded_Volume_"+str(i)] = self.rest_set[('exvo',self.prot[i])].evaluate(False)
            output["Excluded_Volume_Symm_"+str(i)] = self.rest_set[('exvo_symm',self.prot[i])].evaluate(False)
            output["Domain_Linking_"+str(i)] = self.rest_set[('link',self.prot[i])].evaluate(False)
            output["em_all_"+str(i)] = self.rest_set[('em_all',self.prot[i])].evaluate(False)

        output["xlms_Score"] = self.rest_set['xlms'].evaluate(False)
        output["Template_Restraint_Inter"] = self.rest_set['Template_Restraint_Inter'].evaluate(False)         
        output["Template_Restraint_Intra"] = self.rest_set['Template_Restraint_Intra'].evaluate(False)
        output["Temperature"]= self.temp

        flstat.write("%s \n" % output)
        flstat.close()


    #####################################################

    def test_scores(self):
        """This function tests whether there the scores terms
        are similar to the ones calculated using the original IMP version"""
        
        from numpy.testing import assert_approx_equal as aae
        
        #open the file containing the target initial score values
        inputfile=open("data/initial_stats.out","r")
        
        for l in inputfile:
            initialscores=eval(l)       
        
        output={}
        output["Total_Score"]=self.m.evaluate(False)

        aae(float(initialscores["Total_Score"]),output["Total_Score"],7,"Total_Score: test failed")

        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            ln=self.pairs[i][9]
            resid1=self.pairs[i][6][0]
            chain1=self.pairs[i][6][1]
            resid2=self.pairs[i][7][0]
            chain2=self.pairs[i][7][1]
            
            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["xlms_"+label]=ln.evaluate(False)
            
            aae(float(initialscores["xlms_"+label]),
                                                output["xlms_"+label],7,"xlms1: test failed")

            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["Distance_"+label]=IMP.core.get_distance(d0,d1)
            
            aae(float(initialscores["Distance_"+label]),
                                        output["Distance_"+label],7,"Distance: test failed")
            
        for i  in range(len(self.prot)):

            output["Excluded_Volume_"+str(i)] = self.rest_set[('exvo',self.prot[i])].evaluate(False)
            output["Excluded_Volume_Symm_"+str(i)] = self.rest_set[('exvo_symm',self.prot[i])].evaluate(False)
            output["Domain_Linking_"+str(i)] = self.rest_set[('link',self.prot[i])].evaluate(False)
            output["em_all_"+str(i)] = self.rest_set[('em_all',self.prot[i])].evaluate(False)

            aae(float(initialscores["Excluded_Volume_"+str(i)]),
                               output["Excluded_Volume_"+str(i)],7,"Excluded_Volume: test failed")
            aae(float(initialscores["Excluded_Volume_Symm_"+str(i)]),
                          output["Excluded_Volume_Symm_"+str(i)],7,"Excluded_Volume_Symm: test failed")
            aae(float(initialscores["Domain_Linking_"+str(i)]),
                                      output["Domain_Linking_"+str(i)],7,"Domain_Linking: test failed")
            aae(float(initialscores["em_all_"+str(i)]),
                                                      output["em_all_"+str(i)],7,"em_all: test failed")                                    

        output["xlms_Score"] = self.rest_set['xlms'].evaluate(False)
        output["Template_Restraint_Inter"] = self.rest_set['Template_Restraint_Inter'].evaluate(False)         
        output["Template_Restraint_Intra"] = self.rest_set['Template_Restraint_Intra'].evaluate(False)
        output["Temperature"]= self.temp

        aae(float(initialscores["xlms_Score"]),
                                                             output["xlms_Score"],7,"xlms2: test failed") 
        aae(float(initialscores["Template_Restraint_Inter"]),
                           output["Template_Restraint_Inter"],7,"Template_Restraint_Inter: test failed") 
        aae(float(initialscores["Template_Restraint_Intra"]),
                           output["Template_Restraint_Intra"],7,"Template_Restraint_Intra: test failed") 
        aae(float(initialscores["Temperature"]),
                                                     output["Temperature"],7,"Temperature: test failed")


    #####################################################

    def read_pdbs(self,list_pdb_file):
        """read pdbs from an external list file
        create a simplified representation"""

        chains=[] #list of chains
        rigidbodies=[] # list of rigid bodies: more that one chain can belong to different rigid bodies
                       # each line in list_pdb_file will become a rigid body
        chain_id={} #transform chains to id
        id_chain={} #transform id to chains

        hier=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m)) #create an empty hierarchy
        for pdb in list_pdb_file:

            p=IMP.atom.read_pdb(pdb, self.m,
                         IMP.atom.CAlphaPDBSelector())
            ch=IMP.atom.get_by_type(p, IMP.atom.CHAIN_TYPE)

            hier.add_child(p) #add read chains into hierarchy
            chains+=ch
            #rigidbodies.append(ch)

        for c in chains:
            chain_id[c]=IMP.atom.Chain(c).get_id()
            id_chain[IMP.atom.Chain(c).get_id()]=c

        return hier


        #####################################################

    def init_pdb(self):
        self.pdbs=[]
        self.pdbs_symm=[]
        for i in range(self.ncopies):
            pdbfile = "output/trajectory."+str(i)+".pdb"
            self.pdbs.append(pdbfile)
            flpdb=open(pdbfile,'w')
            flpdb.close()
            
            empty_list=[]
            for s in range(self.symmetric_copies):
                pdbfile = "output/trajectory."+str(i)+".symmetry."+str(s)+".pdb"
                empty_list.append(pdbfile)
                flpdb=open(pdbfile,'w')
                flpdb.close()
            self.pdbs_symm.append(empty_list)
                     

        #####################################################

    def write_pdb(self):
        "append models to a pdb file"
        for i in range(self.ncopies):
            flpdb=open(self.pdbs[i],'a')
            IMP.atom.write_pdb(self.prot[i],flpdb)
            
            flpdb.close()            
            for s in range(self.symmetric_copies):
                flpdb=open(self.pdbs_symm[s][i],'a')
                IMP.atom.write_pdb(self.prot_scopies[s][i],flpdb)            
            
                flpdb.close()

        #####################################################
        #####################################################
        #####################################################

    def setup_and_run(self,test=False):

        #####################################################
        #arrays
        #####################################################
        self.rigidbody_list=[]           #list of rigid bodies
        self.prot=[]                     #hierarchies list
        self.rest_set={}                 #restraint set to contain different score terms
        self.ncopies=1                   #number of copies

        
        #####################################################
        #the constants
        #####################################################

        self.delta_r=0.5                 #monte carlo rigid body parameters
        self.delta_a=0.1
        self.delta_dr=0.1
        self.tempmin=0.5
        self.tempmax=2.0
        self.symmetric_copies=1          #number of symmetric copies


        self.m = IMP.Model()

        #####################################################
        #read pdb
        #####################################################


        #for copy in range(0,self.ncopies):
        prot=self.read_pdbs(["data/homology_model_A.pdb"])
        self.prot.append(prot)
        
        #setup the symmetric copies for each state copy
        self.prot_scopies=[]
        for scopy in range(self.symmetric_copies):
          empty_list=[]
          for copy in range(0,self.ncopies):
            prot=self.read_pdbs(["data/homology_model_A.pdb"])
            empty_list.append(prot)
          self.prot_scopies.append(empty_list)


        #####################################################
        #setup representation and sampling
        #####################################################

        #Example two chains and two rigid bodies
        #rigid_body_definition=[ [(-1,-1,"A")],[(-1,-1,"B")] ]

        rigid_body_definition_ref=[  [(1,71,"A")],[(72,203,"A")],
                [(204,251,"A")],[(252,413,"A")],[(414,455,"A")],[(463,815,"A")]]   
        rigid_body_definition_symm=[ [(1,71,"A")],[(72,203,"A")],
                [(204,251,"A")],[(252,413,"A")],[(414,455,"A")],[(463,815,"A")]]  

        self.rigidbody_list.append(self.create_rigid_body(rigid_body_definition_ref,self.prot[0]))
        
        #setup montecarlo movers
        
        nrbs=len(self.rigidbody_list[0])*self.ncopies      #total number of rigid bodies
        mc=self.setup_rigid_body_MonteCarlo(0.1,self.delta_r,self.delta_a,self.delta_dr,False)
        
        #setup the symmetric rigid bodies for each state copy

        self.symmetric_rigid_body_list=[]
        for scopy in range(self.symmetric_copies):
          empty_list=[]
          for copy in range(0,self.ncopies):
                  l=self.create_rigid_body(rigid_body_definition_symm,self.prot_scopies[scopy][copy])
                  empty_list.append(l)      
          self.symmetric_rigid_body_list.append(empty_list)

        #apply symmetries

        dimer=self.read_pdbs(["data/homology_model_A_A.pdb"])  
        
        vsfirst= [IMP.core.XYZ(d).get_coordinates() for d in 
                   IMP.atom.Selection(dimer,chains="A").get_selected_particles()]
        vssecond=[IMP.core.XYZ(d).get_coordinates() for d in 
                   IMP.atom.Selection(dimer,chains="B").get_selected_particles()]        
        
        transformation=[IMP.algebra.get_transformation_aligning_first_to_second(vsfirst,vssecond)]
        
        for icopy in range(0,self.ncopies):
           for scopy in range(0,self.symmetric_copies):
               self.add_symmetry_restraint((self.rigidbody_list[icopy],
                       self.symmetric_rigid_body_list[scopy][icopy]),transformation[scopy])
        
        self.init_stats()
        self.init_pdb()
        
        self.m.update()


        #####################################################
        #adding energy terms
        #####################################################

        #linking the domains

        for icopy in range(0,self.ncopies):
            self.link_domains(self.prot[icopy],[(71,72,"A"),(203,204,"A"),
                           (251,252,"A"),(413,414,"A"),(455,463,"A")],1.0)



        #setting up the cross-link restraints
        
        self.cross_link_ms_simple_symmetrized_create(self.prot[0],
                    self.prot_scopies[0][0],"data/restraints.symm.txt")     
        
        #setting up the em-restraint         
                 
        res=self.add_em_map('data/map.mrc',[(-1,-1,'A',self.prot[0]),
                  (-1,-1,'A',self.prot_scopies[0][0])],'em_all')                 

        #setting up the excluded volume and external barrier

        for icopy in range(0,self.ncopies):
            self.add_excluded_volume(self.prot[icopy],0.1)
            
            slist=[]
            for scopy in range(0,self.symmetric_copies):
                slist.append(self.prot_scopies[scopy][icopy])
            
            self.add_symmetry_excluded_volume(self.prot[icopy],slist,0.1)
        
        #setting up the template restraint
        
        self.m.update()
        s1=IMP.atom.Selection(self.prot[0],chains="A",
                              residue_indexes=range(71,453))
        ps1=s1.get_selected_particles()
        s2=IMP.atom.Selection(self.prot_scopies[0][0],chains="A",
                              residue_indexes=range(71,453))
        ps2=s2.get_selected_particles()              
        allps=ps1+ps2
        self.add_template_restraint(ps1,ps2,"Inter")
        self.add_template_restraint(ps1,ps1,"Intra")


        #####################################################
        #running sampling
        #####################################################
        
        self.nframe=0
        self.temp=0
        #self.write_stats()
        
        if test:
           nstepmax=3
           nsteplowtemp=1
           nstephightemp=1 
           self.test_scores()          
        else:
           nstepmax=10001
           nsteplowtemp=150
           nstephightemp=50
                  
        
        for i in range(1,nstepmax):
            self.nframe=i

            self.temp=self.temp_simulated_annealing(i,nsteplowtemp,nstephightemp)
            mc.set_kt(self.temp)
            mc.optimize(100*nrbs)
            
            #update the orientation to the em map
            #before saving coordinates
            
            
            self.rest_set[('em_all',self.prot[0])].apply_transformation()
            
            
            dict_trans={}
            for r1 in self.rigidbody_list:
               	for r in r1:
               	    # get the initial reference frame
                    rf0 = r.get_reference_frame()
                    dict_trans[r]=rf0
                    r.set_reference_frame_from_members(r.get_member_particle_indexes())
                    
            for r1 in self.symmetric_rigid_body_list:
            
              for s in r1:
               	for r in s:
               	    # get the initial reference frame
                    rf0 = r.get_reference_frame()
                    dict_trans[r]=rf0                    
                    r.set_reference_frame_from_members(r.get_member_particle_indexes())

            
            self.write_pdb()
            
            #reapply the transformations
            for r1 in self.rigidbody_list:
               	for r in r1:                              	
                    r.set_reference_frame(dict_trans[r])
            
            for r1 in self.symmetric_rigid_body_list:
            
              for s in r1:
               	for r in s:            	
                    r.set_reference_frame(dict_trans[r])
            
            #write output file
            self.write_stats()        

if __name__ == "__main__":

    try:
      if sys.argv[1]=="test":
        test=True 
    except:
      test=False
    
    IMP.set_log_level(IMP.TERSE)
    rc=run_class()
    
    rc.setup_and_run(test)
    
