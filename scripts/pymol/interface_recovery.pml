import sys
sys.path.append("/Users/sebastianswanson/Keating/utilities/repos/inhibitory_fragments_structure_prediction/inhib_frag_pred/pymol_scripts")
from helpers import string2sel,drawContact

# Define variables
native_structure = NATIVE_STRUCTURE
colabfold_structure = COLABFOLD_STRUCTURE
fragment_chain_id = CHAIN_ID
native_contacts_string = NATIVE_CONTACTS
colabfold_contacts_string = COLABFOLD_CONTACTS
interface_contact_distance = 8.0

# Viz params
cmd.set('ray_shadows',0)
cmd.set('cartoon_flat_sheets',0)
cmd.set('grid_mode',1)
cmd.set_color('nat_green',[107/255,203/255,119/255])
cmd.set_color('colab_blue',[77/255,150/255,255/255])

nat_dash_length = 1.5
nat_dash_gap = 0.4
nat_dash_radius = 0.2
nat_dash_color = 'yellow'

nonnat_dash_length = 1.5
nonnat_dash_gap = 0.2
nonnat_dash_radius = 0.05
nonnat_dash_color = 'orange'

# Load structures
print(native_structure)
print(colabfold_structure)
print(native_contacts_string)
print(colabfold_contacts_string)

cmd.load(native_structure)
native_name = native_structure[:-4]
cmd.load(colabfold_structure)
colabfold_name = colabfold_structure[:-4]

# Select the interface region
print(f"{native_name} and chain B")
print(f"{colabfold_name} and chain B")
cmd.select("native_fragment",f"{native_name} and chain {fragment_chain_id}")
cmd.select("native_protein",f"{native_name} and not chain {fragment_chain_id}")
cmd.select("colabfold_fragment",f"{colabfold_name} and chain {fragment_chain_id}")
cmd.select("colabfold_protein",f"{colabfold_name} and not chain {fragment_chain_id}")

# Set colors and view CA spheres
cmd.color("white","native_protein")
cmd.color("white","colabfold_protein")
cmd.color("nat_green","native_fragment")
cmd.color("colab_blue","colabfold_fragment")

# seems to help for some reason...
refresh
objs=cmd.get_object_list()
print(objs)

cmd.select("native_fragment_interface",f"br. native_fragment nto. {interface_contact_distance} of native_protein")
cmd.select("native_protein_interface",f"br. native_protein nto. {interface_contact_distance} of native_fragment")
cmd.select("colabfold_fragment_interface",f"br. colabfold_fragment nto. {interface_contact_distance} of colabfold_protein")
cmd.select("colabfold_protein_interface",f"br. colabfold_protein nto. {interface_contact_distance} of colabfold_fragment")

# cmd.show("spheres","native_fragment_interface and n. CA")
# cmd.show("spheres","native_protein_interface and n. CA")
# cmd.show("spheres","colabfold_fragment_interface and n. CA")
# cmd.show("spheres","colabfold_protein_interface and n. CA")
# set sphere_scale, 0.3

# Now define the contacts and draw lines
native_contacts_set = set(native_contacts_string.split(','))
colabfold_contacts_set = set(colabfold_contacts_string.split(','))

print(native_contacts_set)

python
native_contacts_list = []
for contact in native_contacts_set:
    name = drawContact(contact,native_name)
    cmd.hide('labels', name)
    cmd.set('dash_radius', nat_dash_radius, name)
    cmd.set('dash_length', nat_dash_length, name)
    cmd.set('grid_slot', 1, name)
print(native_contacts_list)
python end

print(colabfold_contacts_set)

python
native_contacts_list = []
for contact in colabfold_contacts_set:
    name = drawContact(contact,colabfold_name)
    cmd.hide('labels', name)
    if contact in native_contacts_set:
        cmd.set('dash_radius', nat_dash_radius, name)
        cmd.set('dash_gap', nat_dash_gap, name)
        cmd.set('dash_length', nat_dash_length, name)
        cmd.set('dash_color', nat_dash_color, name)
        cmd.set('grid_slot', 2, name)
    else:
        cmd.set('dash_radius', nonnat_dash_radius, name)
        cmd.set('dash_gap', nonnat_dash_gap, name)
        cmd.set('dash_length', nonnat_dash_length, name)
        cmd.set('dash_color', nonnat_dash_color, name)
        cmd.set('grid_slot', 2, name)
print(native_contacts_list)
python end

cmd.center("native_fragment_interface or colabfold_fragment_interface or native_protein_interface or colabfold_protein_interface")

deselect