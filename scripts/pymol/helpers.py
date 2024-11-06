from pymol import cmd

def string2sel(string):
    split = string.split('_')
    assert len(split) == 2
    return f"chain {split[0]} and resid {int(split[1])}"

def drawContact(string,obj):
    split = string.split('-')
    assert len(split) == 2
    sel1,sel2 = map(string2sel,split)
    name = f"contact_{obj}_{string}"
    print("drawing: ",name)
    print(sel1)
    print(sel2)
    cmd.distance(name=name,selection1=f"{obj} and {sel1} and name CA",selection2=f"{obj} and {sel2} and name CA")
    return name