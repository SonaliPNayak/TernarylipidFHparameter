import MDAnalysis as mda
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from  MDAnalysis.lib import distances
import sys
import os
import multiprocessing as mp
import numpy as np
import argparse
# %pylab inline

# INPUTS
print("Loading inputs...")

# # Topology and trajectory files
# top = '../../mem.gro'         # Path to the topology file
# traj = '../../prod_mem.xtc'   # Path to the trajectory file
# side = 'up'             # Side of the membrane ('up' or 'down')
# bead_type_ref = 'P04'
# bead_type_other = 'P04'
# ref_mol = 'DPPC'
# other_mol = 'DIPC'
# example python rdf-xy.py -s ../../mem.gro -f ../../prod_mem.xtc -leaflet upper -skip 1 -ref_mol DPPC -other_mol DIPC -bead_type_ref PO4
#-bead_type_other PO4 -binwidth

#########defining the arguements needed to run the code and assigning them to variables as input by the user#########
parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='top',help='topology file',required=True)
parser.add_argument('-f',dest='xtc',help='trajectory file',required=True)
parser.add_argument('-leaflet',dest='leaf',help='leaflet:(upper/lower)',required=True)
parser.add_argument('-skip',dest='skip',help='frames to skip',required=True,type=int)
parser.add_argument('-nt',dest='nprocs',help='number of processor to use',default=8,type=int)

parser.add_argument('-ref_mol',dest='sel_ref_mol',help='ref_mol:(DPPC/DIPC/CHOL/GM1)',required=True)
parser.add_argument('-other_mol',dest='sel_other_mol',help='ref_mol:(DPPC/DIPC/CHOL/GM1)',required=True)

parser.add_argument('-bead_type_ref',dest='sel_atom_ref',help='bead_type_ref:(NC3/PO4/GL1/GL2/C1A/C2A/C3A/C4A/D2A(DIPC))',required=True)
parser.add_argument('-bead_type_other',dest='sel_atom_other',help='bead_type_ref:(NC3/PO4/GL1/GL2/C1A/C2A/D2A(DIPC)/C3A/C4A)',required=True)
parser.add_argument('-binwidth',dest='bin_size',help='dr of g(r)',default=0.1,type=float)
args = parser.parse_args()
top=args.top  #"../solute.pdb"
traj=args.xtc  #"../prod.xtc"
side=args.leaf #"upper"
skip=args.skip 

ref_mol=args.sel_ref_mol
other_mol=args.sel_other_mol

bead_type_ref=args.sel_atom_ref
bead_type_other=args.sel_atom_other

test_bin_size=args.bin_size
nprocs = args.nprocs

print("parser done")


# ---- Bead Mapping ---- #
bead_mapping = {
    "DPPC": {
    "NC3": ["N", "C13", "H13A", "H13B", "H13C", "C14", "H14A", "H14B", "H14C", "C15", "H15A", "H15B", "H15C"],
    "PO4": ["C12", "H12A", "H12B", "C11", "H11A", "H11B", "P", "O11", "O12", "O13", "O14"],
    "GL1": ["C1", "HA", "HB", "C2", "HS"],
    "GL2": ["C3", "HX", "HY"],
    "C1A": ["O21", "C21", "C22", "H2R", "H2S", "C23", "H3R", "H3S", "C24", "H4R", "H4S", "C25", "H5R", "H5S"],
    "C2A": ["C26", "H6R", "H6S", "C27", "H7R", "H7S", "C28", "H8R", "H8S", "C29", "H9R", "H9S"],
    "C3A": ["C210", "H10R", "H10S", "C211", "H11R", "H11S", "C212", "H12R", "H12S", "C213", "H13R", "H13S"],
    "C4A": ["C214", "H14R", "H14S", "C215", "H15R", "H15S", "C216", "H16R", "H16S", "H16T"],
    "C1B": ["O31", "C31", "C33", "H3X", "H3Y", "C34", "H4X", "H4Y", "C35", "H5X", "H5Y"],
    "C2B": ["C36", "H6X", "H6Y", "C37", "H7X", "H7Y", "C38", "H8X", "H8Y", "C39", "H9X", "H9Y"],
    "C3B": ["C310", "H10X", "H10Y", "C311", "H11X", "H11Y", "C312", "H12X", "H12Y", "C313", "H13X", "H13Y"],
    "C4B": ["C314", "H14X", "H14Y", "C315", "H15X", "H15Y", "C316", "H16X", "H16Y", "H16Z"]
},
    "DUPC": {
    "NC3": ["N", "C13", "H13A", "H13B", "H13C", "C14", "H14A", "H14B", "H14C", "C15", "H15A", "H15B", "H15C"],
    "PO4": ["C12", "H12A", "H12B", "C11", "H11A", "H11B", "P", "O11", "O12", "O13", "O14"],
    "GL1": ["C1", "HA", "HB", "C2", "HS"],
    "GL2": ["C3", "HX", "HY"],
    "C1A": ["O21", "C21", "C22", "H2R", "H2S", "C23", "H3R", "H3S", "C24", "H4R", "H4S", "C25", "H5R", "H5S"],
    "C2A": ["C26", "H6R", "H6S", "C27", "H7R", "H7S", "C28", "H8R", "H8S", "C29", "H9R"],
    "C3A": ["C210", "H10R", "C211", "H11R", "H11S", "C212", "H12R", "C213", "H13R"],
    "C4A": ["C214", "H14R", "H14S", "C215", "H15R", "H15S", "C216", "H16R", "H16S", "C217", "H17R", "H17S", "C218", "H18R", "H18S", "H18T"],
    "C1B": ["O31", "C31", "C33", "H3X", "H3Y", "C34", "H4X", "H4Y", "C35", "H5X", "H5Y"],
    "C2B": ["C36", "H6X", "H6Y", "C37", "H7X", "H7Y", "C38", "H8X", "H8Y", "C39", "H9X"],
    "C3B": ["C310", "H10X", "C311", "H11X", "H11Y", "C312", "H12X", "C313", "H13X"],
    "C4B": ["C314", "H14X", "H14Y", "C315", "H15X", "H15Y", "C316", "H16X", "H16Y", "H16Z"]
}
}

def get_atom_names(resname, bead):
    return bead_mapping.get(resname, {}).get(bead, [bead])


# Load the MDAnalysis Universe
u = mda.Universe(top, traj)

# Frames to be calculated
# skip = 1
start_f = 0
end_f = u.trajectory.n_frames
#end_f = u.trajectory.n_frames
print(end_f)
frames = np.arange(start_f, end_f, skip)

# Atom selection: Selecting the atoms of interest
#systemsel = f'(resname {ref_mol} and name NC3 PO4 GL1 GL2 C1A C2A C3A C4A D2A D2B D3A D3B) or (resname {other_mol} and name NC3 PO4 GL1 GL2 C1A PO4 D2A D2B D2B D3A D3B C4A)'

# Dynamically generate mapped atom names for selection

mapped_atoms_ref = get_atom_names(ref_mol, bead_type_ref)
mapped_atoms_other = get_atom_names(other_mol, bead_type_other)

systemsel = (
    f'(resname {ref_mol} and name {mapped_atoms_ref} P) or (resname {other_mol} and name {mapped_atoms_other} P)'
)


system = u.select_atoms(systemsel)

# Selecting lipid head groups
lipidheadsgroup = system.select_atoms(f'resname DPPC DUPC and name P')


# Calculating the midpoint (mid-z) of the bilayer
midz = np.mean(lipidheadsgroup.positions[:, 2])  # Average z-coordinate of lipid head groups

# Upper and lower leaflet lipid selection


upperlipid1heads = system.select_atoms(f'resname {ref_mol} and name P and prop z > {midz}')
lowerlipid1heads = system.select_atoms(f'resname {ref_mol} and name P and prop z < {midz}')
upperlipid2heads = system.select_atoms(f'resname {other_mol} and name P and prop z > {midz}')
lowerlipid2heads = system.select_atoms(f'resname {other_mol} and name P and prop z < {midz}')


# Residue lists for each leaflet
upperlipid1residlist = ' '.join(map(str, upperlipid1heads.residues.resids))
lowerlipid1residlist = ' '.join(map(str, lowerlipid1heads.residues.resids))
upperlipid2residlist = ' '.join(map(str, upperlipid2heads.residues.resids))
lowerlipid2residlist = ' '.join(map(str, lowerlipid2heads.residues.resids))


#rdf related inputs
#print(len(frames))
# box dimentions
box_dim = system.dimensions
test_box_Lx = box_dim[0]
test_box_Ly = box_dim[1]

#predefine bin edges

#test_bin_size = 0.01  # Binwidth for RDF 0.1-0.5 AA 
#max_dist = 0.5 * np.sqrt(box_dim[1]**2) #max_range_for rdf -- should approx half of the y axis box length i.e test_box_Ly
max_dist = 0.5 * np.sqrt(box_dim[0]**2 + box_dim[1]**2)
bins = np.arange(0, max_dist + test_bin_size, test_bin_size) # predefined constand bining -for consitent output
#max_dist = 0.5 * np.sqrt(box_dim[1]**2) #max_range_for rdf


def compute_distances(DPPC, DUPC):
    # Extract 2D positions
    g1_coord = DPPC.positions[:, 0:2]  # 2D position coord as input for 2d rdf ( coordinates of ref particle - g1)
    g2_coord = DUPC.positions[:, 0:2]  # 2D position for 2d rdf ( coord of other particle -g2 )

    # Add a third column with z = 0 to create 3D positions
    g1_with_z = np.hstack((g1_coord, np.zeros((g1_coord.shape[0], 1))))  # 3D position with z = 0
    g2_with_z = np.hstack((g2_coord, np.zeros((g2_coord.shape[0], 1))))  # 3D position with z = 0

    # Calculate capped distances
    pairs, dist = distances.capped_distance(
        g1_with_z,  # 3D position as per module requirement
        g2_with_z,  # 3D position
        max_dist,
        box=box_dim 
    )

    return pairs, dist

def histogram_distances(distance_list, bins):
    # this is the list of bins in which to calculate
    #bins = np.arange(0, max_dist+bin_size, bin_size) # bins predefined
    hist, bin_edges = np.histogram( distance_list, bins=bins )
    return hist, bin_edges


def get_gofr(hist,bin_edges,rho,N_l1):
    
    bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.0
    dr = bin_edges[1]-bin_edges[0]
    num_ref_particle = N_l1 # must be reference particle -- criteria for proper averaging
    denominator = 2.*np.pi*bin_centers*dr*rho*(num_ref_particle) # rho must be calculated for particle whose rdf we are calcuating around the ref particle 
                                          #example- if we want rdf of DIPC around DPPC then we will calculate the rho of DIPC only 
                                          #-- because its distribution is what we are seeking to normalize
    gofr = (hist/denominator)
    
    return gofr, bin_centers


def RDF2d(frame):
    """Calculate 2D RDF for a given frame."""
    u = mda.Universe(top, traj)
    u.trajectory[frame]
    box = u.atoms.dimensions
    system = u.select_atoms(systemsel)
    #print(f'resname {ref_mol} and name {bead_type_ref} and resid {upperlipid1residlist}')
    # Select lipid headgroups based on side
    if side == 'up':
        mapped_l1_atoms = get_atom_names(ref_mol, bead_type_ref)
        mapped_l2_atoms = get_atom_names(other_mol, bead_type_other)
        l1 = system.select_atoms(f'resname {ref_mol} and name {mapped_l1_atoms} and resid {upperlipid1residlist}')
        l2 = system.select_atoms(f'resname {other_mol} and name {mapped_l2_atoms} and resid {upperlipid2residlist}')
    else:
        mapped_l1_atoms = get_atom_names(ref_mol, bead_type_ref)
        mapped_l2_atoms = get_atom_names(other_mol, bead_type_other)
        l1 = system.select_atoms(f'resname {ref_mol} and name {mapped_l1_atoms} and resid {lowerlipid1residlist}')
        l2 = system.select_atoms(f'resname {other_mol} and name {mapped_l2_atoms} and resid {lowerlipid2residlist}')

    #print(l1) # ref mol
    l2 # other mol

    N_l1 = len(l1.indices)  # number of ref beads
    #print(N_l1)
    N_l2 = len(l2.indices) # number of other beads
    N = N_l1 + N_l2 
    rho = N_l2 / (test_box_Lx * test_box_Ly)  # avg density of particle whose rdf we are calcuating around the ref particle
    #print(rho)
    #print(N_DPPC)
    pairs, dist = compute_distances(l1, l2)
    dist_hist, bin_edges = histogram_distances(dist, bins)
    gofr, bin_centers = get_gofr(dist_hist, bin_edges, rho, N_l1)
    
    return gofr, bin_centers

if __name__ == '__main__':
    # Parameters
    nprocs = nprocs  # Number of processors
    gofr_list = []  # List to store RDFs
    frames = np.arange(start_f, end_f, skip)  # Frame range

    print(f'Initiating multiprocessing with {nprocs} processors...')
    
    # Create a pool of workers
    with mp.Pool(processes=nprocs) as pool:
        results = pool.map(RDF2d, frames)

    # Unpack results
    gofr, bin_centers = zip(*results)

    # Calculate average and standard deviation
    g_r_mean = np.mean(gofr, axis=0)
    g_r_std = np.std(gofr, axis=0)
    radii = np.mean(bin_centers, axis=0)  # Bin centers from the first result

    # Save results to a file
    data = np.array([radii, g_r_mean, g_r_std]).T
    np.savetxt(f'g_r_{bead_type_ref}_{bead_type_other}_{ref_mol}_{other_mol}_{test_bin_size}.dat', data, header='Radius RDF_mean RDF_std')
    
    print(f"RDF computation completed. Results saved to 'g_r_{bead_type_ref}_{bead_type_other}_{ref_mol}_{other_mol}_{test_bin_size}.dat'")
