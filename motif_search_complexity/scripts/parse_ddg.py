import pandas as pd
import glob
import sys
import os

def get_data( filename ):
    try:
        with open(filename) as fd:
            header = 0
            for line in fd:
                if line.startswith("SEQUENCE"):continue
                if line.startswith("SCORE: total_score"):
                    l = line.split();
                    for i, ii in enumerate(l):
                        if ii == "bb_clash":
                            header = i
                            break
                    continue
                if header == 0: header = 2
                l = line.strip().split();
                return [
                    group,
                    float(l[header]),
                    l[-1].split("_")[0]
                ]
    except KeyboardInterrupt:
        sys.exit()
    except:
        return ["e",100000,"eeee"]

if __name__ == '__main__':
    workfolder = sys.argv[1]
    files = glob.glob(os.path.join(workfolder,"ddgout","*","*"))
    sys.stdout.write("listed {} files\n".format(len(files)))
    data = {}
    for f in files:
        d = get_data(f)
        data.setdefault("cluster", []).append(d[0])
        data.setdefault("ddg", []).append(d[1])
        data.setdefault("str", []).append(d[2])
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(workfolder,"ddg_match.csv"), index=False)
