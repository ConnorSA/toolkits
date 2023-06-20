import re
import sys
# file="/home/physics/phspnt/project_work/local_castep/multiphase_ti/tf_benchmarks/hcp_tf_CASTEP_k0p010/tf_castep_k0p010.castep"

file=sys.argv[1]
with open(file) as f:
    all_lines = f.readlines()

start_lines=[]
stop_lines=[]
for i,l in enumerate(all_lines):
    if re.search('CCC   AA    SSS  TTTTT  EEEEE  PPPP', l):
        start_lines.append(i-2)
    if re.search('Overall parallel efficiency rating', l):
        stop_lines.append(i+4)
# split files
for i, (start, stop) in enumerate(zip(start_lines,stop_lines)):
    lines=all_lines[start:stop]
    with open(f"castep_split_{str(i).zfill(3)}.castep", "w", encoding = 'utf-8') as f:
        f.write("".join(lines))


