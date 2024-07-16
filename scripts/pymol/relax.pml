#import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)

# Define variables
native_structure = "native.pdb"
relaxcst_structure = "relax_cst.pdb"
relax_structure = "relax_nocst.pdb"
fragment_chain_id = "B"
interface_contact_distance = 8.0

# Viz params
cmd.set('ray_shadows',0)
cmd.set('cartoon_flat_sheets',0)
cmd.set_color('nat_green',[107/255,203/255,119/255])
cmd.set_color('colab_blue',[77/255,150/255,255/255])

# Load structures
print(native_structure)
print(relaxcst_structure)
print(relax_structure)

cmd.load(native_structure)
native_name = native_structure[:-4]
cmd.load(relaxcst_structure)
relaxcst_name = relaxcst_structure[:-4]
cmd.load(relax_structure)
relax_name = relax_structure[:-4]

# Make selections
cmd.select("native_fragment",f"{native_name} and chain {fragment_chain_id}")
cmd.select("relaxcst_fragment",f"{relaxcst_name} and chain {fragment_chain_id}")
cmd.select("relax_fragment",f"{relax_name} and chain {fragment_chain_id}")

# Set colors and view CA spheres
cmd.color("white","*")
cmd.color("nat_green","native_fragment")
cmd.color("purple","relaxcst_fragment")
cmd.color("colab_blue","relax_fragment")

# seems to help for some reason...
refresh
objs=cmd.get_object_list()
print(objs)

# Assign to grid slots
cmd.set('grid_mode',1)
cmd.set('grid_slot', 1, native_name)
cmd.set('grid_slot', 2, relaxcst_name)
cmd.set('grid_slot', 3, relax_name)

cmd.center("native_fragment or relaxcst_fragment or relax_fragment")

deselect