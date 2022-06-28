def read_headered_data(file_in : str):
    with open(file_in) as f:
        lines=f.readlines()
    labels=lines[0].split()
    mydict={}
    for i, val in enumerate(labels):
        mydict[val] = []
    for i, val in enumerate(lines[1:-1]):
        for j, val2 in enumerate(labels):
            mydict[val2].append(val.split()[j])
    return mydict
