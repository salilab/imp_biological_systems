#!/bin/env python

import math
import sys
import os


# ROOT benchmark directory
dir=sys.argv[1]

# list of sub-directories with test
test_dirs = [x for x in os.listdir(dir) if x.startswith('testdir')]

# initialize dictionaries for stats
ave_drms={};     ave2_drms={}
ave_drms_ter={}; ave2_drms_ter={};
ave_pa={};       ave2_pa={};
ave_pd={};       ave2_pd={};
ave_kda={};      ave2_kda={};
ave_Ida={};      ave2_Ida={};
ave_sigma0={};   ave2_sigma0={};
ntests={}

# log file
log=open(sys.argv[2],"w")

# cycle on tests
for test_dir in test_dirs:
 for line in open(dir+"/"+test_dir+"/best_posterior.dat"):
     d=eval(line)
     drms=d["DRMS"]
     drms_ter=d["DRMS-ter"]
     sigma0=d["Sigma0"]
     pa=d["Placement_angle"]
     pd=d["Placement_distance"]
     kda=d["kda"]
     Ida=d["Ida"]
 for line in open(dir+"/"+test_dir+"/data.out"):
     d=eval(line)
     nunits=len(d["PDBs"])
     pdb=d["PDBs"][0].split("/")[-1].split("_")[0]
     sigma0_ref=d["Sigma0"]
     kda_ref=d["Kda"]
     Ida_ref=d["Ida"]
 for line in open(dir+"/"+test_dir+"/selected-data.out"):
     d=eval(line)
     ndata=d["DataNumber"]

# label for the test
 id=(pdb,nunits,sigma0_ref,ndata)
 if id in ave_drms:
    ave_drms[id]+=drms;          ave2_drms[id]+=drms**2
    ave_drms_ter[id]+=drms_ter;  ave2_drms_ter[id]+=drms_ter**2
    ave_pa[id]+=pa;              ave2_pa[id]+=pa**2
    ave_pd[id]+=pd;              ave2_pd[id]+=pd**2
    ave_kda[id]+=(kda-kda_ref)**2
    ave_Ida[id]+=(Ida-Ida_ref)**2
    ave_sigma0[id]+=(sigma0-sigma0_ref)**2
    ntests[id]+=1.0
 else:
    ave_drms[id]=drms;           ave2_drms[id]=drms**2
    ave_drms_ter[id]=drms_ter;   ave2_drms_ter[id]=drms_ter**2
    ave_pa[id]=pa;               ave2_pa[id]=pa**2
    ave_pd[id]=pd;               ave2_pd[id]=pd**2
    ave_kda[id]=(kda-kda_ref)**2
    ave_Ida[id]=(Ida-Ida_ref)**2
    ave_sigma0[id]=(sigma0-sigma0_ref)**2
    ntests[id]=1.0

# print result in a nice order
for nunits in [3,4]:
 # get uniq list of pdb with units in alphabetic order
 pdbs=[]
 for key in ave_drms:
  if(key[1]==nunits): pdbs.append(key[0])
 for pdb in sorted(set(pdbs)):
  for sigma in [0.001,0.01]:
   for ndata in [0.5,1.0]:
    id=(pdb,nunits,sigma,ndata)
    if id not in ave_drms:
        continue
    log.write("%6s %4d %6.4lf %2.1f   " % (id[0],id[1],id[2],id[3]))
    log.write("dRMS %6.3lf "        % (ave_drms[id]/ntests[id]))
    log.write("(%6.3lf )  "         % (math.sqrt(ave2_drms[id]/ntests[id]-(ave_drms[id]/ntests[id])**2)/math.sqrt(ntests[id])))
    log.write("dRMD_ter %6.3lf "    % (ave_drms_ter[id]/ntests[id]))
    log.write("(%6.3lf )  "         % (math.sqrt(ave2_drms_ter[id]/ntests[id]-(ave_drms_ter[id]/ntests[id])**2)/math.sqrt(ntests[id]))) 
    log.write("pdist  %6.3lf "      % (ave_pd[id]/ntests[id]))
    log.write("(%6.3lf )  "         % (math.sqrt(ave2_pd[id]/ntests[id]-(ave_pd[id]/ntests[id])**2)/math.sqrt(ntests[id])))
    log.write("pangle %8.3lf "      % (ave_pa[id]/ntests[id]))
    log.write("(%6.3lf )  "         % (math.sqrt(ave2_pa[id]/ntests[id]-(ave_pa[id]/ntests[id])**2)/math.sqrt(ntests[id])))
    log.write("kda %6.3lf  "        % (math.sqrt(ave_kda[id]/ntests[id])))
    log.write("Ida %6.3lf  "        % (math.sqrt(ave_Ida[id]/ntests[id])))
    log.write("sigma0 %12.6lf  "    % (math.sqrt(ave_sigma0[id]/ntests[id])))
    log.write("ntests %4d\n"        % (ntests[id]))
