import pandas as pd
import glob
import sys
import os

def get_data( filename ):
    line = open(filename).readline().replace(", ","").strip().split()
    if len(line) == 4:
        pdb  = os.path.split(line[2])[-1].split(".")[0].split("_")
        return [
            os.path.split(filename)[-1].replace(".pdb",""),
            os.path.split(os.path.split(filename)[0])[-1],
            line[1], pdb[0], pdb[1], line[-1].replace(",","-")
        ]
    else:
        return ["error","e",0.0,"eeee","e","none"]

if __name__ == '__main__':
    workfolder = sys.argv[1]
    files = glob.glob(os.path.join(workfolder,"structures","*","*"))
    sys.stdout.write("listed {} files\n".format(len(files)))
    data = {"str":[], "cluster":[], "rmsd":[], "pdb":[], "chain":[],"range":[]}
    for f in files:
        sys.stdout.write(f+"\n")
        d = get_data(f)
        if d[1] == "e":
            sys.stdout.write("error\n")
        sys.stdout.flush()
        data["str"].append(d[0])
        data["cluster"].append(d[1])
        data["rmsd"].append(d[2])
        data["pdb"].append(d[3])
        data["chain"].append(d[4])
        data["range"].append(d[5])
    df = pd.DataFrame(data).sort_values(["rmsd"])
    df.to_csv(os.path.join(workfolder, "master_search.csv"), index=False )
