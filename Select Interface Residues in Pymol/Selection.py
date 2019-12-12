
def selection(structure_name):
    from pymol import cmd
    import string
    cmd.load("enter desired PDB file directory" % structure_name, 'protein')
    
    counter = 0
    chain_list_A=[]
    chain_list_B = []
    
    
    A_file = open('chain residues directory' % structure_name)
    filecontents_A = A_file.read().split()

    for val in filecontents_A:
            chain_list_A.append(val)
    
    A_file.close() 
    
         
    array = []
    for val in chain_list_A:
            array.append(int(val))
    for i in array:
        counter = counter + 1
        cmd.select("residuedata", "resi %s%s" % (i,' and chain A'))
        #cmd.show('sphere', "residuedata")
        cmd.color ('red', "residuedata")
    print(counter)
    array = []
    
    
