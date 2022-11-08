# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules
import csv
import json

with open("test_input/4_dissemination_index.json") as json_conf : 
    conf = json.load(json_conf)

# translate CARD antibiotics names into antibiotic classes
ab_file = "utilities/antibiotic_ARO_translation.tsv"
ab_dict = {}
with open(ab_file, 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        ab_dict[row["CARD_Best_Hit_Drug_Class_UNIQUE"]] = row["Antibiotic_Class"]

# open the PLSDB working file 
file = conf["arg_db"]
db_dict_ab, db_dict_ab_split, db_dict_rm = {}, {}, {}
with open(file, 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        if row["Plot"] == "1":
           #### ANTIBIOTICS STUFF
           # translate CARD antibiotic list into antibiotic classes
           antibiotics = row["CARD_Best_Hit_Drug_Class_original"].split(";")
           antibiotics_f = []
           for ab in antibiotics:
               antibiotics_f.append(ab_dict[ab])
               antibiotics_f = list(set(antibiotics_f))
           # check if antibiotic class is in the dictionary
           for ab_f in antibiotics_f:
               if ab_f not in db_dict_ab.keys():
                   db_dict_ab[ab_f] = {"antibiotic efflux": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic inactivation": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target alteration": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target protection": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target replacement": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "reduced permeability to antibiotic": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "resistance by absence": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}}
                       
                   db_dict_ab_split[ab_f] = {}
                   db_dict_ab_split[ab_f]["Specific"] = {"antibiotic efflux": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic inactivation": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target alteration": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target protection": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target replacement": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "reduced permeability to antibiotic": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "resistance by absence": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}}
                       
                   db_dict_ab_split[ab_f]["Multiple"] = {"antibiotic efflux": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic inactivation": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target alteration": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target protection": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "antibiotic target replacement": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "reduced permeability to antibiotic": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}, \
                   "resistance by absence": {"Total": 0, "Spreading": 0, "Mutant_t": 0, "Regulator_t": 0,"Mutant_s": 0, "Regulator_s": 0}}
                   
               # count total ARG for each antibiotic class and % of spreading
               db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Total"] += 1 
               db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
               # save mutant/regulator info
               if row["Mutation_tag"] == "Y":
                   db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_t"] += 1
                   db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
               if row["Regulator_tag"] == "Y":
                   db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_t"] += 1
                   db_dict_ab[ab_f][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
                   
               # if resistance mechanism is specific for only one ab
               if len(antibiotics_f) == 1:
                   db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Total"] += 1 
                   db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
                   # save mutant/regulator info
                   if row["Mutation_tag"] == "Y":
                       db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_t"] += 1
                       db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
                   if row["Regulator_tag"] == "Y":
                       db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_t"] += 1
                       db_dict_ab_split[ab_f]["Specific"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
               if len(antibiotics_f) >= 2:
                   db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Total"] += 1 
                   db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
                   # save mutant/regulator info
                   if row["Mutation_tag"] == "Y":
                       db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_t"] += 1
                       db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Mutant_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
                   if row["Regulator_tag"] == "Y":
                       db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_t"] += 1
                       db_dict_ab_split[ab_f]["Multiple"][row["CARD_Best_Hit_Resistance_Mechanism"]]["Regulator_s"] += int(row["Gene_in_spreadingband_ARO_COG"])
            
           #### RESISTANCE MECHANISM STUFF
           rm = row["CARD_Best_Hit_Resistance_Mechanism"]
           if rm not in db_dict_rm.keys():
               db_dict_rm[rm] = {"Global": {"Total": 0, "Spreading": 0, "Total_specific": 0, "Spreading_specific": 0}, \
                                 "No_mut": {"Total": 0, "Spreading": 0, "Total_specific": 0, "Spreading_specific": 0}, \
                                 "No_reg": {"Total": 0, "Spreading": 0, "Total_specific": 0, "Spreading_specific": 0}, \
                                 "No_mut_reg": {"Total": 0, "Spreading": 0, "Total_specific": 0, "Spreading_specific": 0}}
           # count total ARG for each resistance mechanism and % of spreading
           db_dict_rm[rm]["Global"]["Total"] += 1
           db_dict_rm[rm]["Global"]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
           db_dict_rm[rm]["No_mut"]["Total"] += 1
           db_dict_rm[rm]["No_mut"]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
           db_dict_rm[rm]["No_reg"]["Total"] += 1
           db_dict_rm[rm]["No_reg"]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
           db_dict_rm[rm]["No_mut_reg"]["Total"] += 1
           db_dict_rm[rm]["No_mut_reg"]["Spreading"] += int(row["Gene_in_spreadingband_ARO_COG"])
           # if resistance mechanism is specific for only one ab
           if len(antibiotics_f) == 1:
               db_dict_rm[rm]["Global"]["Total_specific"] += 1
               db_dict_rm[rm]["Global"]["Spreading_specific"] += int(row["Gene_in_spreadingband_ARO_COG"])
               db_dict_rm[rm]["No_mut"]["Total_specific"] += 1
               db_dict_rm[rm]["No_mut"]["Spreading_specific"] += int(row["Gene_in_spreadingband_ARO_COG"])
               db_dict_rm[rm]["No_reg"]["Total_specific"] += 1
               db_dict_rm[rm]["No_reg"]["Spreading_specific"] += int(row["Gene_in_spreadingband_ARO_COG"])
               db_dict_rm[rm]["No_mut_reg"]["Total_specific"] += 1
               db_dict_rm[rm]["No_mut_reg"]["Spreading_specific"] += int(row["Gene_in_spreadingband_ARO_COG"])
           # remove if mutant
           if row["Mutation_tag"] == "Y":
               db_dict_rm[rm]["No_mut"]["Total"] -= 1
               db_dict_rm[rm]["No_mut"]["Spreading"] -= int(row["Gene_in_spreadingband_ARO_COG"])
               db_dict_rm[rm]["No_mut_reg"]["Total"] -= 1
               db_dict_rm[rm]["No_mut_reg"]["Spreading"] -= int(row["Gene_in_spreadingband_ARO_COG"])
               if len(antibiotics_f) == 1:
                   db_dict_rm[rm]["No_mut"]["Total_specific"] -= 1
                   db_dict_rm[rm]["No_mut"]["Spreading_specific"] -= int(row["Gene_in_spreadingband_ARO_COG"])
                   db_dict_rm[rm]["No_mut_reg"]["Total_specific"] -= 1
                   db_dict_rm[rm]["No_mut_reg"]["Spreading_specific"] -= int(row["Gene_in_spreadingband_ARO_COG"])
           # remove if regulator
           if row["Regulator_tag"] == "Y":
               db_dict_rm[rm]["No_reg"]["Total"] -= 1
               db_dict_rm[rm]["No_reg"]["Spreading"] -= int(row["Gene_in_spreadingband_ARO_COG"])
               db_dict_rm[rm]["No_mut_reg"]["Total"] -= 1
               db_dict_rm[rm]["No_mut_reg"]["Spreading"] -= int(row["Gene_in_spreadingband_ARO_COG"])
               if len(antibiotics_f) == 1:
                   db_dict_rm[rm]["No_reg"]["Total_specific"] -= 1
                   db_dict_rm[rm]["No_reg"]["Spreading_specific"] -= int(row["Gene_in_spreadingband_ARO_COG"])
                   db_dict_rm[rm]["No_mut_reg"]["Total_specific"] -= 1
                   db_dict_rm[rm]["No_mut_reg"]["Spreading_specific"] -= int(row["Gene_in_spreadingband_ARO_COG"])

# output saving
out_handle_ab = open(conf["output_ab_tsv"], "w")
out_handle_rm = open(conf["output_rm_tsv"], "w")

#### ANTIBIOTICS STUFF
out_handle_ab.write("Antibiotic\tMechanism\tHits\tMutants\tSpreader_Mutants\tRegulators\tSpreader_Regulator\tDissemination_index\tARO_COG_in_band\n")
for key_ab in db_dict_ab.keys():
    hits, mutants_t, regulators_t, mutants_s, regulators_s, band = 0,0,0,0,0,0
    for key_ab_rm in db_dict_ab[key_ab].keys():
        try:
            index = (float(db_dict_ab[key_ab][key_ab_rm]["Spreading"])/db_dict_ab[key_ab][key_ab_rm]["Total"])
        except ZeroDivisionError:
            index = 0
        out_handle_ab.write(key_ab+"\t"+key_ab_rm+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Total"])+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Mutant_t"])+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Mutant_s"])+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Regulator_t"])+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Regulator_s"])+"\t"+str(index)+"\t"+str(db_dict_ab[key_ab][key_ab_rm]["Spreading"])+"\n")
        hits += db_dict_ab[key_ab][key_ab_rm]["Total"]
        mutants_t += db_dict_ab[key_ab][key_ab_rm]["Mutant_t"]
        regulators_t += db_dict_ab[key_ab][key_ab_rm]["Regulator_t"]
        mutants_s += db_dict_ab[key_ab][key_ab_rm]["Mutant_s"]
        regulators_s += db_dict_ab[key_ab][key_ab_rm]["Regulator_s"]
        band += db_dict_ab[key_ab][key_ab_rm]["Spreading"]
    mi_global = float(band/hits)
    out_handle_ab.write(key_ab+"\tall mechanisms\t"+str(hits)+"\t"+str(mutants_t)+"\t"+str(mutants_s)+"\t"+str(regulators_t)+"\t"+str(regulators_s)+"\t"+str(mi_global)+"\t"+str(band)+"\n\n")

for specificity in ["Specific", "Multiple"]:
    for key_ab in db_dict_ab_split.keys():
        hits, mutants_t, regulators_t, mutants_s, regulators_s, band = 0,0,0,0,0,0
        for key_ab_rm in db_dict_ab_split[key_ab][specificity].keys():
            try:
                index = (float(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Spreading"])/db_dict_ab_split[key_ab][specificity][key_ab_rm]["Total"])
            except ZeroDivisionError:
                index = 0
            # out_handle_ab.write(key_ab+"\t"+key_ab_rm+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Total"])+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Mutant_t"])+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Mutant_s"])+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Regulator_t"])+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Regulator_s"])+"\t"+str(index)+"\t"+str(db_dict_ab_split[key_ab][specificity][key_ab_rm]["Spreading"])+"\n")
            hits += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Total"]
            mutants_t += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Mutant_t"]
            regulators_t += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Regulator_t"]
            mutants_s += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Mutant_s"]
            regulators_s += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Regulator_s"]
            band += db_dict_ab_split[key_ab][specificity][key_ab_rm]["Spreading"]
        try:
            mi_global = float(band/hits)
        except ZeroDivisionError:
            mi_global = 0
        out_handle_ab.write(key_ab+"\tall mechanisms\t"+str(hits)+"\t"+str(mutants_t)+"\t"+str(mutants_s)+"\t"+str(regulators_t)+"\t"+str(regulators_s)+"\t"+str(mi_global)+"\t"+str(band)+"\t"+specificity+"\n")

#### RESISTANCE MECHANISM STUFF
orders = ["Global", "No_mut", "No_reg", "No_mut_reg"]
for order in orders:
    out_handle_rm.write(order+"\n")
    out_handle_rm.write("Resistance_mechanism\tCount\tARO_COG_in_band\tDissemination_index\n")
    for key_rm in db_dict_rm.keys():
        try:
            index = float(db_dict_rm[key_rm][order]["Spreading"])/db_dict_rm[key_rm][order]["Total"]
        except ZeroDivisionError:
            index = 0
        out_handle_rm.write(key_rm+"\t"+str(db_dict_rm[key_rm][order]["Total"])+"\t"+str(db_dict_rm[key_rm][order]["Spreading"])+"\t"+str(index)+"\n")
        try:
            index_specific = float(db_dict_rm[key_rm][order]["Spreading_specific"])/db_dict_rm[key_rm][order]["Total_specific"]
        except ZeroDivisionError:
            index_specific = 0
        out_handle_rm.write(key_rm+"_specific\t"+str(db_dict_rm[key_rm][order]["Total_specific"])+"\t"+str(db_dict_rm[key_rm][order]["Spreading_specific"])+"\t"+str(index_specific)+"\n")
        try:
            index_multiple = float((db_dict_rm[key_rm][order]["Spreading"]-db_dict_rm[key_rm][order]["Spreading_specific"]))/(db_dict_rm[key_rm][order]["Total"]-db_dict_rm[key_rm][order]["Total_specific"])
        except ZeroDivisionError:
            index_multiple = 0
        out_handle_rm.write(key_rm+"_multiple\t"+str(db_dict_rm[key_rm][order]["Total"]-db_dict_rm[key_rm][order]["Total_specific"])+"\t"+str(db_dict_rm[key_rm][order]["Spreading"]-db_dict_rm[key_rm][order]["Spreading_specific"])+"\t"+str(index_multiple)+"\n")
    out_handle_rm.write("\n")

out_handle_ab.close()
out_handle_rm.close()
