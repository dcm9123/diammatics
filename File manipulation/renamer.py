import os

path = "/bulk/sycuro_bulk/daniel/diabetes/UC_UT_collaboration/16S_data/16S_whole_from_four_primers/fixed_assembly_benchling/"
new_file=""
for file in os.listdir(path):
    if file.endswith(".csv"):
        name = file.split(".")[0]
        if len(name)==5:
            name = name.split("seq")[1]
            name = "seq0"+name
            new_file = name+".csv"
        elif len(name)==4:
            name = name.split("seq")[1]
            name = "seq00"+name
            new_file = name+".csv"
        else:
            continue
        os.rename(os.path.join(path,file),os.path.join(path,new_file))
    else:
        continue
