# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules
import csv
import math
import numpy as np
import json

# Load the configuration file
with open("test_input/3_spreading_band_detection.json") as json_conf : 
    conf = json.load(json_conf)

db = conf["db"]
db_r = conf["db_r"]
#value determining the range of host GC to consider spreading
spreading_GC_value = conf["spreading_GC_value"]
spreading_GC_dev = conf["spreading_GC_dev"]
GC_bin = conf["GC_bin"]

#open the db and create a dictionary with the following structure
#COG = {hit:mechanism}
db_dict = {}
items, total_COGs = [], []
with open(db, 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        if row["ARO-COG ID"]:
            aro = row["ARO-COG ID"].split(";")
            aro.sort()
            aro = ";".join(aro)
            if aro not in db_dict.keys():
                db_dict[aro] = {}
            db_dict[aro][row["Hit_accession"]] = row["CARD_Best_Hit_Resistance_Mechanism"]
            items.append(row["CARD_Best_Hit_Resistance_Mechanism"])
            items = list(set(items))
            total_COGs.append(aro)
            total_COGs = list(set(total_COGs))

#open the db_r and create a dictionary with the
#following structure "Hit_COG_accessions:[ARG_GC, Host_GC]"
plsdb_hits = {}
with open(db_r, 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        if row:
            if row["ARO-COG ID"] != "":
                family_name = row["ARO-COG ID"].split(";")
                family_name.sort()
                family_name_def = ";".join(family_name)
                if family_name_def not in plsdb_hits.keys():
                    plsdb_hits[family_name_def] = []
                if row["ARG_GC"] != " " and row["Host_GC"] != " " and row["ARG_GC"] != "" and row["Host_GC"] != "":
                    plsdb_hits[family_name_def].append([row["ARG_GC"], row["Host_GC"], row["Hit_accession"]])

#sort ARG_GC values in the desired range of ARG_GC% and save them as lists in a new 
#dictionary
#create a new dictionary with the total number of values, the mean and SD of 
#number values per ranges of ARG_GC%
plsdb_hits_2 = {}
plsdb_hits_2_mean = {}
ranges = list(np.arange(0, 100+GC_bin, GC_bin))
for key in plsdb_hits.keys():
    #sort plsdb_hits[key] considering the first value of the sublists:
    #[ARG_GC, Host_GC] 
    plsdb_hits[key] = sorted(plsdb_hits[key], key = lambda x: x[0])
    for i in range(1,len(ranges)):
        subset = []
        for gc_value in plsdb_hits[key]:
            gc_value = [float(gc_value[0]),float(gc_value[1]), str(gc_value[2])]
            if gc_value[0] > ranges[i-1] and gc_value[0] <= ranges[i]:
                subset.append(gc_value)
        if key not in plsdb_hits_2.keys():
            plsdb_hits_2[key], plsdb_hits_2_mean[key] = [], []
        if len(subset) >= 1:
            plsdb_hits_2[key].append(subset)
            plsdb_hits_2_mean[key].append(len(subset))
        else:
            # add foo values to stop the bands if the GC_bin difference between bands is >= 2
            plsdb_hits_2[key].append([[ranges[i],ranges[i],"foo"], [ranges[i],ranges[i],"foo"]])
    
    ntotal = sum(plsdb_hits_2_mean[key])
    mean = np.mean(plsdb_hits_2_mean[key])
    std = np.std(plsdb_hits_2_mean[key])
    plsdb_hits_2_mean[key] = [ntotal, mean, std]
    
#get the mobilization band using the desired spreading_GC_value, std and GC bin
plsdb_hits_3 = {}
for key2 in plsdb_hits_2:
    spreading_GC_value_updated = spreading_GC_value+spreading_GC_dev
    mobilized_band, max_band = [], []
    min_band, line_count = 0, 0
    max_band_size = 0
    if key2 not in plsdb_hits_3.keys():
        plsdb_hits_3[key2] = []
    for upper_band in range(0,len(plsdb_hits_2[key2])):
        host_gc = [element[1] for element in plsdb_hits_2[key2][upper_band]]
        if (max(host_gc) - min(host_gc)) >= (spreading_GC_value_updated - spreading_GC_dev):
            max_band.append(upper_band)
            line_count = line_count+1
            spreading_GC_value_updated = max(host_gc) - min(host_gc)
            new_mobilized_band = [element1[2] for element1 in plsdb_hits_2[key2][upper_band]]
            mobilized_band = mobilized_band+new_mobilized_band
            if len(new_mobilized_band) >= max_band_size:
                max_band_size = len(new_mobilized_band)
        else:
            if len(mobilized_band) >= 1:
                max_min_band = range(min_band+1,max_band[0]-1)[::-1]
                for lower_band in max_min_band:
                    host_gc = [element2[1] for element2 in plsdb_hits_2[key2][lower_band]]
                    if (max(host_gc) - min(host_gc)) >= (spreading_GC_value_updated - spreading_GC_dev):
                        spreading_GC_value_updated = max(host_gc) - min(host_gc)
                        new_mobilized_band = [element3[2] for element3 in plsdb_hits_2[key2][lower_band]]
                        mobilized_band = mobilized_band+new_mobilized_band
                        line_count = line_count+1
                    else:
                        break
                #only accept bands with >= number of genes than expected &
                #bands >= 5 genes
                if max_band_size >= plsdb_hits_2_mean[key2][1] and len(mobilized_band) >= 5:                   
                    plsdb_hits_3[key2].append(mobilized_band)
                    min_band = upper_band
            line_count = 0
            spreading_GC_value_updated = spreading_GC_value+spreading_GC_dev
            mobilized_band, max_band = [], []
            
            host_gc = [element[1] for element in plsdb_hits_2[key2][upper_band]]
            if len(mobilized_band) == 0 and (max(host_gc) - min(host_gc)) >= (spreading_GC_value_updated - spreading_GC_dev):
                max_band.append(upper_band)
                line_count = line_count+1
                spreading_GC_value_updated = max(host_gc) - min(host_gc)
                new_mobilized_band = [element1[2] for element1 in plsdb_hits_2[key2][upper_band]]
                mobilized_band = mobilized_band+new_mobilized_band
                if len(new_mobilized_band) >= max_band_size:
                    max_band_size = len(new_mobilized_band)

#compute the whole mobilization index & the mobilization index for each item in
#the list items
mi_dict = {}
mi_dict["Global"] = []
out_handle = open(conf["output_tsv"], "w")
out_handle.write("ARO_name\tGlobal\t"+"\t".join(items)+"\t"+"\t".join(items)+"\n")
for COG in plsdb_hits_3.keys():
    total_COGs.remove(COG)
    count = 0.0
    for spreading_band in plsdb_hits_3[COG]:
        spreading_band_nr = list(set(spreading_band))
        count = count + len(spreading_band_nr)
    total_mi = round(count/len(list(db_dict[COG].keys())),2)
    mi_dict["Global"].append(total_mi)
    out_handle.write(COG+"\t"+str(total_mi))
    # compute the mobilization index per item
    spreading_bands_mechanism = []
    for item in items:
        if item not in mi_dict.keys():
            mi_dict[item] = []
        db_count, plsdb_count = 0.0,0.0
        current_spreading_band_mechanism = []
        for db_gene in db_dict[COG]:
            if db_dict[COG][db_gene] == item:
                db_count = db_count + 1
        for plsdb_band in plsdb_hits_3[COG]:
            plsdb_band_nr = list(set(plsdb_band))
            for plsdb_hit in plsdb_band_nr:
                if db_dict[COG][plsdb_hit] == item:
                    plsdb_count = plsdb_count +1
                    current_spreading_band_mechanism.append(plsdb_hit)
                    current_spreading_band_mechanism = list(set(current_spreading_band_mechanism))
        spreading_bands_mechanism.append(current_spreading_band_mechanism)
        try:
            item_mi = round((plsdb_count/db_count),2)
            out_handle.write("\t"+str(item_mi))
            mi_dict[item].append(item_mi)
        except:
            out_handle.write("\t")
    out_handle.write("\t")
    for print_spreading_band in spreading_bands_mechanism:
        out_handle.write(";".join(print_spreading_band)+"\t")
    out_handle.write("\n")

# print COGs that are not redundant
for absent_COG in total_COGs:
    out_handle.write(absent_COG+"\t0")
    mi_dict["Global"].append(0)
    for item in items:
        check = None
        for db_gene in db_dict[COG]:
            if db_dict[COG][db_gene] == item:
                check = True
                break
        if check == True:
            out_handle.write("\t0")
            mi_dict[item].append(0)
        elif check == None:
            out_handle.write("\t")
    out_handle.write("\n")
    
# print summary values
out_handle.write("Mean\t"+str(np.mean(mi_dict["Global"])))
for item in items:
    out_handle.write("\t"+str(np.mean(mi_dict[item])))
out_handle.write("\n")

stdev = []
out_handle.write("STDEV\t"+str(np.std(mi_dict["Global"])))
stdev.append(np.std(mi_dict["Global"]))
for item in items:
    out_handle.write("\t"+str(np.std(mi_dict[item])))
    stdev.append(np.std(mi_dict[item]))
out_handle.write("\n")
    
count_v = []
out_handle.write("Count\t"+str(len(mi_dict["Global"])))
count_v.append(len(mi_dict["Global"]))
for item in items:
    out_handle.write("\t"+str(len(mi_dict[item])))
    count_v.append(len(mi_dict[item]))
out_handle.write("\n")

out_handle.write("STDERRMEAN")
for j in range(0,len(stdev)):
    out_handle.write("\t"+str(stdev[j]/math.sqrt(count_v[j])))
out_handle.write("\n")


out_handle.close()
