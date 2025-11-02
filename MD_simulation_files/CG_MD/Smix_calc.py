import numpy as np
import MDAnalysis as mda
import MDAnalysis.lib.distances 
import MDAnalysis.lib.NeighborSearch as NS
from scipy.spatial import Voronoi
import sys, os
from p_tqdm import p_map
import pickle
import argparse

#########defining the arguements needed to run the code and assigning them to variables as input by the user#########
parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='top',help='topology file',required=True)
parser.add_argument('-f',dest='xtc',help='trajectory file',required=True)
parser.add_argument('-leaflet',dest='leaf',help='leaflet:(upper/lower)',required=True)
parser.add_argument('-skip',dest='skip',help='frames to skip',required=True,type=int)
parser.add_argument('-nt',dest='nprocs',help='number of processor to use',default=8,type=int)
parser.add_argument('-protein',dest="prot",help="is protein(APP C99) present in simulation:YES/NO",default="NO")
args = parser.parse_args()
top=args.top  #"../solute.pdb"
traj=args.xtc  #"../prod.xtc"
side=args.leaf #"upper"
# top="../mem.pdb"
# traj="../prod_nopbc.xtc"
# side="upper"
skip=args.skip 
nprocs=args.nprocs
protein=args.prot
prot_martini_resname=["TRP","TYR","PHE","HIS","HIH","GLN","ASN","SER","THR","ARG","LYS","ASP","GLU","CYS","ILE","LEU","MET","PRO","HYP","VAL","ALA","GLY"]
lipid1 ='DPPC'
lipid2 ='DIPC'
lipid3 ='CHOL'

#########Trajectory Loading and selecting the system for analysis #############################
u = mda.Universe(top,traj)
start_frame=0
end_frame=u.trajectory.n_frames
frames=np.arange(start_frame,end_frame)[::skip]

########Selecting the upper leaflet and lower leaflet residues#################################
lipidheadgroups=u.select_atoms("resname DPPC DIPC and name PO4")
cholheadgroups=u.select_atoms("resname CHOL and name ROH")
midz=np.mean(lipidheadgroups.positions[:,2])

if side=="upper":
    DPPChead=u.select_atoms(f"resname DPPC and name PO4 and prop z > {midz}")
    DIPChead=u.select_atoms(f"resname DIPC and name PO4 and prop z > {midz}")
elif side=="lower":
    DPPChead=u.select_atoms(f"resname DPPC and name PO4 and prop z < {midz}")
    DIPChead=u.select_atoms(f"resname DIPC and name PO4 and prop z < {midz}")
DPPC_resid_str=" ".join(DPPChead.residues.resids.astype('str'))
DIPC_resid_str=" ".join(DIPChead.residues.resids.astype('str'))

def voronoi_tessel(frame):
    u = mda.Universe(top,traj)
    u.trajectory[int(frame)]
    if protein=="YES":
        systemsel='(name BB and (resid 16-19 or resid 32-35 or resid 56-59 or resid 72-75)) or (resname DPPC DIPC and name PO4 C2A C2B D2A D2B) or (resname CHOL and name ROH R2)'
        system=u.select_atoms(systemsel)
        if side=="upper":
            proai_sel = 'resid 16-19 and name BB'
            probi_sel = 'resid 56-59 and name BB'
        elif side=="lower":
            proai_sel = 'resid 32-35 and name BB'
            probi_sel = 'resid 72-75 and name BB'
        proai=system.select_atoms(proai_sel)
        probi=system.select_atoms(probi_sel)

        #print(proai_sel)
    elif protein=="NO":
        systemsel='(resname DPPC DIPC and name PO4 C2A C2B D2A D2B) or (resname CHOL and name ROH R2)'
    system=u.select_atoms(systemsel)
    cholheadgroups=system.select_atoms("resname CHOL and name ROH")
    box=u.atoms.dimensions
    DPPChead=system.select_atoms(f"resname DPPC and name PO4 and resid {DPPC_resid_str}")
    DIPChead=system.select_atoms(f"resname DIPC and name PO4 and resid {DIPC_resid_str}")
    DPPC_DIPC_head=DPPChead+DIPChead 
    DPPC_tail=system.select_atoms(f"resname DPPC and name C2A C2B and resid {DPPC_resid_str}")
    DIPC_tail=system.select_atoms(f"resname DIPC and name D2A D2B and resid {DIPC_resid_str}")
    ns_lipids=NS.AtomNeighborSearch(cholheadgroups, box=box)
    CHOLsel = ns_lipids.search(DPPC_DIPC_head,12.0)
    CHOL_resid_str=" ".join(CHOLsel.residues.resids.astype('str'))
    CHOL_tail=system.select_atoms(f"resname CHOL and name R2 and resid {CHOL_resid_str}")
    #proai=system.select_atoms(proai_sel)
    #probi=system.select_atoms(probi_sel)
    if protein=="YES":
        lipids=DPPC_tail+DIPC_tail+CHOL_tail+proai+probi
    elif protein=="NO":
        lipids=DPPC_tail+DIPC_tail+CHOL_tail
    

    ############Extracting coordinates#############
    coord=lipids.positions
    ###########Extracting xy coordinates and residue names##########
    xy_resname=[[coord[i][0],coord[i][1],lipids[i].resname] for i in range(len(coord))]
   
    ############Creating multiple PBC box around the system#######
    x_box=box[0]
    y_box=box[1]
    xplus=[[xy_resname[i][0]+x_box,xy_resname[i][1],xy_resname[i][2]] for i in range(len(xy_resname))]
    xminus=[[xy_resname[i][0]-x_box,xy_resname[i][1],xy_resname[i][2]] for i in range(len(xy_resname))]
    xy_resname_x=xy_resname+xplus+xminus
    xyplus=[[xy_resname_x[i][0],xy_resname_x[i][1]+y_box,xy_resname_x[i][2]] for i in range(len(xy_resname_x))]
    xyminus=[[xy_resname_x[i][0],xy_resname_x[i][1]-y_box,xy_resname_x[i][2]] for i in range(len(xy_resname_x))]
    xy_resname_pbc=xy_resname_x+xyplus+xyminus
    ############create a list of only x and y coordinate from the replicated pbc list########
    xy_pbc=[[xy_resname_pbc[i][0],xy_resname_pbc[i][1]] for i in range(0,len(xy_resname_pbc))]
    #######Constructing voronoi tesselation##############
    vor=Voronoi(xy_pbc)
    vertices=vor.vertices
    ridge_points=vor.ridge_points
    r_length=len(ridge_points)
    ######lipid-lipid contact count############
    Lpd1_Lpd1_I = 0
    Lpd2_Lpd2_I = 0
    Lpd3_Lpd3_I = 0
    Lpd1_Lpd2_I = 0
    Lpd1_Lpd3_I = 0
    Lpd2_Lpd3_I = 0

    Lpd1_Lpd1_E = 0
    Lpd2_Lpd2_E = 0
    Lpd3_Lpd3_E = 0
    Lpd1_Lpd2_E = 0
    Lpd1_Lpd3_E = 0
    Lpd2_Lpd3_E = 0

    for k in range (0,r_length) :
        ridge_k = ridge_points[k]
        Li = xy_resname_pbc[int(ridge_k[0])]
        Lj = xy_resname_pbc[int(ridge_k[1])]

        if 0 < Li[0] < x_box and 0 < Li[1] < y_box and 0 < Lj[0] < x_box and 0 < Lj[1] < y_box :
        
            if Li[2] == lipid1 and Lj[2] == lipid1:
                Lpd1_Lpd1_I = Lpd1_Lpd1_I + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid2:
                Lpd2_Lpd2_I = Lpd2_Lpd2_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid3:
                Lpd3_Lpd3_I = Lpd3_Lpd3_I + 1
                
            if Li[2] == lipid1 and Lj[2] == lipid2:
                Lpd1_Lpd2_I  = Lpd1_Lpd2_I + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid1:
                Lpd1_Lpd2_I = Lpd1_Lpd2_I + 1

            if Li[2] == lipid1 and Lj[2] == lipid3:
                Lpd1_Lpd3_I  = Lpd1_Lpd3_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid1:
                Lpd1_Lpd3_I = Lpd1_Lpd3_I + 1

            if Li[2] == lipid2 and Lj[2] == lipid3:
                Lpd2_Lpd3_I  = Lpd2_Lpd3_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid2:
                Lpd2_Lpd3_I = Lpd2_Lpd3_I + 1                
                
#Lipids at the EDGE of the box                
                
        if 0 <= Li[0] < x_box and 0 <= Li[1] < y_box or 0 <= Lj[0] < x_box and 0 <= Lj[1] < y_box :

            if Li[2] == lipid1 and Lj[2] == lipid1:
                Lpd1_Lpd1_E = Lpd1_Lpd1_E + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid2:
                Lpd2_Lpd2_E = Lpd2_Lpd2_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid3:
                Lpd3_Lpd3_E = Lpd3_Lpd3_E + 1
                
            if Li[2] == lipid1 and Lj[2] == lipid2:
                Lpd1_Lpd2_E  = Lpd1_Lpd2_E + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid1:
                Lpd1_Lpd2_E = Lpd1_Lpd2_E + 1

            if Li[2] == lipid1 and Lj[2] == lipid3:
                Lpd1_Lpd3_E  = Lpd1_Lpd3_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid1:
                Lpd1_Lpd3_E = Lpd1_Lpd3_E + 1

            if Li[2] == lipid2 and Lj[2] == lipid3:
                Lpd2_Lpd3_E  = Lpd2_Lpd3_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid2:
                Lpd2_Lpd3_E = Lpd2_Lpd3_E + 1    
        #Total = LipidsInside + (Lipids including EDGES - Lipids Inside)/2 -----> Correction for over counting the lipids in periodic images
    Lpd1_Lpd1 = Lpd1_Lpd1_I + (Lpd1_Lpd1_E - Lpd1_Lpd1_I)/2
    Lpd2_Lpd2 = Lpd2_Lpd2_I + (Lpd2_Lpd2_E - Lpd2_Lpd2_I)/2
    Lpd3_Lpd3 = Lpd3_Lpd3_I + (Lpd3_Lpd3_E - Lpd3_Lpd3_I)/2
    Lpd1_Lpd2 = Lpd1_Lpd2_I + (Lpd1_Lpd2_E - Lpd1_Lpd2_I)/2
    Lpd1_Lpd3 = Lpd1_Lpd3_I + (Lpd1_Lpd3_E - Lpd1_Lpd3_I)/2
    Lpd2_Lpd3 = Lpd2_Lpd3_I + (Lpd2_Lpd3_E - Lpd2_Lpd3_I)/2

    sum_bonds = Lpd1_Lpd1 + Lpd2_Lpd2 + Lpd3_Lpd3 + Lpd1_Lpd2 + Lpd1_Lpd3 + Lpd2_Lpd3

#Considering only Similar Lipid (SL) and Dissimilar Lipid (DL) Bonds
    # SL = Lpd1_Lpd1 + Lpd2_Lpd2 + Lpd3_Lpd3
    # DL = Lpd1_Lpd2 + Lpd1_Lpd3 + Lpd2_Lpd3

    # #Calculating Fractions    
    # X_SL = float(SL)/float(sum_bonds)  #Similar Lipid
    # X_DL  = float(DL)/float(sum_bonds) #Dissimilar Lipid
    sum_bonds=Lpd1_Lpd1+Lpd2_Lpd2+Lpd1_Lpd2
    SL = Lpd1_Lpd1 + Lpd2_Lpd2 
    DL = Lpd1_Lpd2

    #Calculating Fractions    
    X_SL = float(SL)/float(sum_bonds)  #Similar Lipid
    X_DL  = float(DL)/float(sum_bonds) #Dissimilar Lipid

#Mixing Entropy
    mix_entropy = -((X_SL * np.log2(X_SL)) +( X_DL * np.log2(X_DL)))

#Calculating Averages
    sum_bonds = Lpd1_Lpd1 + Lpd1_Lpd2 + Lpd2_Lpd2
    avg_bonds = float(sum_bonds)/float(len(xy_resname))

    return Lpd1_Lpd1, Lpd2_Lpd2, Lpd3_Lpd3, Lpd1_Lpd2, Lpd1_Lpd3, Lpd2_Lpd3, sum_bonds, avg_bonds, mix_entropy


results = p_map(voronoi_tessel, frames, num_cpus=nprocs)
pickle.dump( results, open(f'lipid_lipid_contact_Smix_{side}_skip{skip}.p', "wb") )


