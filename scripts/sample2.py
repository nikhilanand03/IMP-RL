from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
#import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
#import IMP.saxs
import os
import sys
import warnings 
warnings.filterwarnings("ignore")




def recenter(mol):
    "recenter the system using the coordinate of a chain specified by its chain_id"
    m=mol.get_model()
    sel = IMP.atom.Selection(mol,resolution=IMP.atom.ALL_RESOLUTIONS)
    ps = sel.get_selected_particles()
    rb=IMP.atom.create_rigid_body(ps)
    rbcoord=rb.get_coordinates()
    rot=IMP.algebra.get_identity_rotation_3d()
    tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
    trans=tmptrans.get_inverse()
    IMP.core.transform(rb,trans)
    IMP.core.RigidBody.teardown_particle(rb)
    m.remove_particle(rb.get_particle_index())
    # change this back!Now rigid member coords need to be optimized as no longer in rigid body
    for p in ps:
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)
    return(mol)

def get_radius_center_residue(mol):
    sel = IMP.atom.Selection(mol,resolution=1)
    ps = sel.get_selected_particles()
    max_dist=0.0
    min_dist=100000.0
    for p in ps:
        crd= IMP.core.XYZ(p).get_coordinates()
        dist=IMP.algebra.get_distance(crd,IMP.algebra.Vector3D(0.,0.,0.))
        if dist<min_dist:
            min_dist=dist
            center_res=p
        if dist>max_dist:
            max_dist=dist
    return(max_dist,center_res)

def save_rmf(root_hier, i):
    IMP.atom.show_with_representations(root_hier)
    # output to RMF
    fname = '/home/varun/PrISM_3/two_protein_complex/modeling/test_' + str(i) + '.rmf'
    rh = RMF.create_rmf_file(fname)
    IMP.rmf.add_hierarchy(rh, root_hier)
    IMP.rmf.save_frame(rh)


runType = sys.argv[1] # Specify test or prod
name = sys.argv[2]  # Name of complex
runID = sys.argv[3]   # Specify the run number
run_output_dir = '/home/varun/PrISM_3/two_protein_complex/modeling/' + name + '/run_' + str(runID)
rex_max_temp = 2.5


if runType == "test":
    num_frames = 1000
elif runType == "prod":
    num_frames = 8000

# Identify data files
xl_data = '/home/varun/PrISM_3/two_protein_complex/inputs/' + name + '.xlinks'
import pandas as pd


# Restraint weights
xl_weight = 10

# Topology File
topology_file = "/home/varun/PrISM_3/two_protein_complex/inputs/" + name + "_topology.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)
print(bs.execute_macro())
# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans= 2.0,
                                  max_rb_rot= 0.1)
                                  
molecules = t.get_components()


save_rmf(root_hier, 0)

receptor = root_hier.get_children()[0].get_children()[0]
ligand = root_hier.get_children()[0].get_children()[1]
mols = [receptor,ligand]

'''
receptor=recenter(receptor) #HACK
ligand=recenter(ligand)
'''
(receptor_radius,receptor_center_residue)=get_radius_center_residue(receptor)
(ligand_radius,ligand_center_residue)=get_radius_center_residue(ligand)





prot_name = 'A'
fixed_particles=[]
fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot_name).get_selected_particles()

fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,IMP.core.BallMover,
                                          IMP.pmi.TransformMover])



#####################################################
##################### RESTRAINTS ####################
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are connected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of all molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=[root_hier])
evr.set_weight(0.03) 
evr.add_to_model()
output_objects.append(evr)

print("Excluded volume restraint applied")





# -------------------------
# %%%%% CROSSLINKING RESTRAINT
#
# Restrains two particles via a distance restraint based on
# an observed crosslink.
#
# First, create the crosslinking database from the input file
# The "standard keys" correspond to a crosslink csv file of the form:
#
# Protein1,Residue1,Protein2,Residue2
# A,18,G,24
# A,18,G,146
# A,50,G,146
# A,50,G,171
# A,50,G,189
#
# This restraint allows for ambiguity in the crosslinked residues,
# a confidence metric for each crosslink and multiple states.
# See the PMI documentation or the MMB book chapter for a
# full discussion of implementing crosslinking restraints.

# This first step is used to translate the crosslinking data file.
# The KeywordsConverter maps a column label from the xl data file
# to the value that PMI understands.
# Here, we just use the standard keys.
# One can define custom keywords using the syntax below.
# For example if the Protein1 column header is "prot_1"
# xldbkc["Protein1"]="prot_1"

# The CrossLinkDataBase translates and stores the crosslink information
# from the file "xl_data" using the KeywordsConverter.


xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                              converter=xldbkc)
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb, # The crosslink database.
                length=8,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.05,          # This adds a linear term to the scoring function
                label="adh")       # Scaling factor for the restraint score.



output_objects.append(xlr)


# add distance to point restraint from receptor to ligand, to keep ligand floating away

ligand_center_residue_index=IMP.atom.Residue(ligand_center_residue).get_index()
receptor_center_residue_index=IMP.atom.Residue(receptor_center_residue).get_index()
dsr= IMP.pmi.restraints.basic.DistanceRestraint(root_hier=root_hier,tuple_selection1=(ligand_center_residue_index,ligand_center_residue_index,'B',0),tuple_selection2=(receptor_center_residue_index,receptor_center_residue_index,'A',0),distancemin=5,distancemax=15)
dsr.add_to_model()
output_objects.append(dsr)


print("Cross-linking restraint applied")


#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))
# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(root_hier,max_translation=50)
                                    # hierarchies_included_in_collision=fixed_set1_core)
                                    # hierarchies_included_in_collision=fixed_set2)
# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(100)


# Now, add all of the other restraints to the scoring function to start sampling
evr.add_to_model()
xlr.add_to_model()



#print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        crosslink_restraints=[xlr],
        # This allows viewing the crosslinks in Chimera. Also, there is not inter-protein ADH crosslink available. Hence it is not mentioned in this list
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory=run_output_dir,      # The output directory for this sampling run.
        output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)            # Total number of frames to run / write to the RMF file.
        #test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)



# Ok, now we finally do the sampling!
rex.execute_macro()

# Outputs are then analyzed in a separate analysis script.
