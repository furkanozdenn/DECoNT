import numpy as np
from os import listdir
import pdb
import pandas as pd
from tqdm import tqdm
'''
xhmm data type is of:
chr2,88861845,89177166,DUP

cnvnator data type is of:
chrY	10409901	CNVnator_dup_8625	C	<DUP>	.	PASS	END=10434500;SVTYPE=DUP;SVLEN=24600;IMPRECISE;natorRD=3.28605;natorP1=0;natorP2=2.00376e+09;natorP3=0;natorP4=2.06321e+09;natorQ0=1
'''





#function to resolve breakpoints chromosome-wise
def match_vectors(chrn_xhmm, chrn_cnvnator):
    mapping_list = []

    for index_xhmm in range(len(chrn_xhmm)):
        for index_cnvnator in range(len(chrn_cnvnator)):
            #continue if start base of xhmm entry is greater than end base of cnvnator entry
            #print(int(chrn_xhmm[index_xhmm].split(",")[1]), int(chrn_cnvnator[index_cnvnator].split("\t")[7].split(";")[0].split("=")[1]))
            if(int(chrn_xhmm[index_xhmm].split(",")[1]) >= int(chrn_cnvnator[index_cnvnator].split("\t")[7].split(";")[0].split("=")[1])):
                continue
            #add to list if regions are overlapping

            if(int(chrn_xhmm[index_xhmm].split(",")[2]) >= int(chrn_cnvnator[index_cnvnator].split("\t")[1])):
                mapping_list.append((index_xhmm,index_cnvnator))
            
    current_xhmm_index = 0
    current_cnvnator_index = 0

    if not mapping_list:
        return chrn_xhmm

    last_xhmm_index = mapping_list[0][0]
    temp_cnv_list = []

    for i in range(0, len(mapping_list)):

        current_xhmm_index = mapping_list[i][0]
        current_cnvnator_index = mapping_list[i][1]

        if(current_xhmm_index == last_xhmm_index):
            temp_cnv_list.append(chrn_cnvnator[current_cnvnator_index].split("\t")[4])
            
        
        else:
            if not temp_cnv_list:
                temp_cnv_list.append("NAN")
            resulting_vote = max(set(temp_cnv_list), key=temp_cnv_list.count)
            temp_cnv_list = []
            temp_cnv_list.append(chrn_cnvnator[current_cnvnator_index].split("\t")[4])
            chrn_xhmm[last_xhmm_index] = chrn_xhmm[last_xhmm_index].rstrip() + "," + resulting_vote + "\n"
            last_xhmm_index = current_xhmm_index

            if i == len(mapping_list) - 1:
                if not temp_cnv_list:
                    temp_cnv_list.append("NAN")
                resulting_vote = max(set(temp_cnv_list), key=temp_cnv_list.count)
                chrn_xhmm[last_xhmm_index] = chrn_xhmm[last_xhmm_index].rstrip() + "," + resulting_vote + "\n"

            continue

        if i == len(mapping_list) - 1:
            if not temp_cnv_list:
                temp_cnv_list.append("NAN")
            resulting_vote = max(set(temp_cnv_list), key=temp_cnv_list.count)
            chrn_xhmm[last_xhmm_index] = chrn_xhmm[last_xhmm_index].rstrip() + "," + resulting_vote + "\n"

        last_xhmm_index = current_xhmm_index
    
    if sample_name=="HG00309":
        #pdb.set_trace()
        pass
    return chrn_xhmm
        


chr_wise_data_path = "/home/furkan/deepXCNV/XHMM/chr_wise_data_liftover_xhmm_hiseq4000_na12878/"
path_to_cnvnator = "/home/furkan/deepXCNV/CNVNATOR_WGS_CALLS/CNVNATOR_WGS_DATA/"
path_to_xhmm = "/home/furkan/deepXCNV/XHMM/data_liftover_hiseq4000_na12878/"


files_in_xhmm = listdir(path_to_xhmm)
files_in_cnvnator = listdir(path_to_cnvnator)



for i in tqdm(range(len(files_in_xhmm))):
    sample_name = files_in_xhmm[i].split(".")[0]
    print(sample_name, i)

    #read xhmm data line-by-line
    with open(path_to_xhmm+files_in_xhmm[i]) as f:
        xhmm_data = f.readlines()

    #read cnvnator data line-by-line
    with open(path_to_cnvnator+sample_name+".cnvnator.illumina_high_coverage.20190825.sites.vcf") as f:
        cnvnator_data = f.readlines()

    chrwise_xhmm_data = [[] for i in range(24)]
    chrwise_cnvnator_data = [[] for i in range(24)]

    for i in range(24):
        if i == 22:
            chrname= "chrX"
            chrwise_cnvnator_data[i] = [x for x in cnvnator_data if x.split("\t")[0]==chrname]
            chrwise_xhmm_data[i] = [x for x in xhmm_data if x.split(",")[0]==chrname]

        elif i == 23:
            chrname= "chrY"
            chrwise_cnvnator_data[i] = [x for x in cnvnator_data if x.split("\t")[0]==chrname]
            chrwise_xhmm_data[i] = [x for x in xhmm_data if x.split(",")[0]==chrname]

        else:
            chrname= "chr"+str(i+1)
            chrwise_cnvnator_data[i] = [x for x in cnvnator_data if x.split("\t")[0]==chrname]
            chrwise_xhmm_data[i] = [x for x in xhmm_data if x.split(",")[0]==chrname]
    #print(chrwise_xhmm_data[4])
    #print(match_vectors(chrwise_xhmm_data[4], chrwise_cnvnator_data[4]))


    for i in range(24):
        if i == 22:
            chrname= "chrX"
            labeled = match_vectors(chrwise_xhmm_data[i], chrwise_cnvnator_data[i])
            with open(chr_wise_data_path + chrname + "_" + sample_name + ".txt", "a") as myfile:
                for row in labeled:
                    myfile.write(row)

        elif i == 23:
            chrname= "chrY"
            labeled = match_vectors(chrwise_xhmm_data[i], chrwise_cnvnator_data[i])
            with open(chr_wise_data_path + chrname + "_" + sample_name + ".txt", "a") as myfile:
                for row in labeled:
                    myfile.write(row)

        else:
            chrname= "chr"+str(i+1)
            labeled = match_vectors(chrwise_xhmm_data[i], chrwise_cnvnator_data[i])
            with open(chr_wise_data_path + chrname + "_" + sample_name + ".txt", "a") as myfile:
                for row in labeled:
                    myfile.write(row)  
#print(chrwise_xhmm_data[4])



