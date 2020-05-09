import numpy as np
import pdb
from tqdm import tqdm

#path to xhmm .xcnv file
xncv_path = "/home/furkan/deepXCNV/XHMM/DATA_chaisson_na19240.xcnv"
output_path = "/home/furkan/deepXCNV/XHMM/data_liftover_chaisson_na19240/"

#read all data line-by-line
with open(xncv_path) as f:
    content = f.readlines()

#SAMPLE  CNV     INTERVAL        KB      CHR     MID_BP  TARGETS NUM_TARG        Q_EXACT Q_SOME  Q_NON_DIPLOID   Q_START Q_STOP  MEAN_RD MEAN_ORIG_RD
print(content[1].split("\t"))

#traverse .xcnv file line-by-line to parse cnv data
rows_for_sample = []
for i in tqdm(range(1,len(content))):
    splitted_line = content[i].split("\t")
    sample_name = splitted_line[0]
    detected_cnv = splitted_line[1]
    interval_info = splitted_line[2]
    chr = interval_info.split(":")[0]
    start_index = interval_info.split(":")[1].split("-")[0]
    end_index = interval_info.split(":")[1].split("-")[1]
    line_to_write = chr + "," + start_index + "," + end_index + "," +detected_cnv

    if i == len(content)-1:
        next_name = "end_of_file"
    else:
        next_name = content[i+1].split("\t")[0]

    if next_name == sample_name:
        rows_for_sample.append(line_to_write)
    else:
        with open(output_path+sample_name+'.xhmm_data.txt', 'w') as f:
            for line in rows_for_sample:
                f.write("%s\n" % line)
        rows_for_sample = []










