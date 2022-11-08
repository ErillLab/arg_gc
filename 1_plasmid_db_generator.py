# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules
import csv
import json
import os
import subprocess
import time
from Bio import Entrez

def gb2faa(input_file, col_name, output_file):
    
    '''
    Download fasta_cds_aa using a list of tsv-formated GenBank accession numbers
    '''
    
    # save the GenBank accession numbers into a list
    db_ids = []
    # ISO-8859-1 for handling utf-8 problems
    with open(input_file, 'r', encoding='ISO-8859-1') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter = "\t")
        for row in csvreader:
            db_ids.append(row[col_name])
            db_ids = list(set(db_ids))
    
    # save fasta_cds_aa of each GenBank into a single FASTA file
    out_handle = open(output_file, "w")
    for nucleotide_id in db_ids:
        print ("Downloading "+nucleotide_id+"...")
        handle = Entrez.efetch(db="nucleotide", id = nucleotide_id, rettype="fasta_cds_aa")
        out_handle.write(handle.read())
        time.sleep(2)
        handle.close()
    out_handle.close()


def blastdb(input_file, output_file):
    
    '''
    Construct the blast_db using the bash script makeblastdb.sh
    '''
    
    bashCommand = "bash makeblastdb.sh "+input_file+" "+output_file
    subprocess.call(bashCommand.split())
 

############################################################################

# Load the configuration file
with open("test_input/1_plasmid_db_generator.json") as json_conf : 
    conf = json.load(json_conf)

# tell to who I am
Entrez.email = conf["email"]
Entrez.api_key = conf["api_key"]

# copy the makeblastdb.sh script to the current dir
bashCommand = "cp utilities/makeblastdb.sh ."
subprocess.call(bashCommand.split())

# get faa file
gb2faa(conf["plsdb_csv"], "ACC_NUCCORE", conf["db_fasta"])
# & construct the BLAST_DB
blastdb(conf["db_fasta"], conf["db_blastp_name"])
