# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules
import csv, json, os, time, subprocess
from Bio import  Entrez, SearchIO, SeqIO
from Bio.SeqUtils import GC
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO

def protein_retrieval(protein_acc, output_file):
    
    """ Download all protein sequences from NCBI RefSeq/GenBank database and 
        save them into a FASTA file.
        The headers are modified to obtain the following structure:
        >ProteinIdentifier__Species
    """
    
    out_handle = open(output_file, "w")
    handle = Entrez.efetch(db="protein", id=protein_acc, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    for record in records:
        out_handle.write(">"+record.id+"__"+record.description.split("[")[1].split("]")[0].replace(" ","_")+"\n")
        out_handle.write(str(record.seq)+"\n")
    out_handle.close()
    handle.close()
    time.sleep(2)


def blast_search(query_file, database, cutoff, nhits, coverage):
    
    """ Local BLASTp search to detect orthologues.
        Receives a query protein accession, an e-value cut off a coverage cut off
        and the maximum number of hits to be retrieved.

        Returns a list containing the protein accessions for the BLASTP hits.
    """

    #blastp against local database
    blastp_cline = NcbiblastpCommandline(query=query_file, db=database,\
                                        evalue=cutoff, outfmt='"6 std qcovs"',\
                                        num_alignments = nhits)
    stdout, sterr = blastp_cline()
    
    blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",\
                                   fields=['std','qcovs'])
        
    #return the hits passing retrieved thresholds & save them into a dict
    #with acc as keys and E-val as values
    blastp_results = {}
    for blast_record in blast_records:
        for alignment in blast_record:
            if alignment.query_coverage >= coverage:              
                blastp_hit = str(alignment.id_all[0])
                nucc = blastp_hit.split("lcl|")[-1].split(".")[0]
                prot = blastp_hit.split("lcl|")[-1].split(".")[1].split("_prot_")[-1]
                acc = prot+"__"+nucc
                for hsp in alignment.hsps:
                    evalue = hsp.evalue
                blastp_results[acc] = evalue
    
    #remove the query FASTA file
    os.remove(query_file)   
    
    return blastp_results


def HMM_annotation(protein_file, hmm_db, evalue):
    
    """ Annotate a given protein file with single FASTA seq & return an HMM annotation
    """
    
    bashCommand = "hmmscan --tblout cog_temp_arg.txt -E " +str(evalue)+" "+hmm_db+" "+protein_file
    subprocess.call(bashCommand.split())
    cog_file = "cog_temp_arg.txt"
    HMM_acc = ""
    with open(cog_file) as fp2:
        line = fp2.readline()
        cnt = 1
        while line:
            line = fp2.readline()
            cnt += 1
            if not line.startswith("#") and line != "":
                if "fam" in hmm_db:
                    HMM_acc = line.split(".")[0].split(" ")[-1]
                    break
                else:
                    HMM_acc = line.split(" ")[0]
                    break           
    os.remove("cog_temp_arg.txt")
            
    return HMM_acc


def genome_record_retrieval(ortholog_acc, nucleotide, sleepy):
    
    """ Takes a protein accession as an input. Retrieves its IPG record.

        The idea here is to obtain, prioritarily, data from complete genome 
        records if they exist, from RefSeq (AC_ and NC_ accessions) . If no 
        RefSeq is available, then select complete genome records from GenBank 
        (AE, CP, CY accessions). Otherwise, select contigs or WGS scaffolds from 
        RefSeq (NT_, NW_, NZ_). If that fails, get contigs or WGS scaffolds from 
        GenBank (AAAA-AZZZ). Only when nothing else is available, select direct 
        GenBank submissions (U, AF, AY, DQ).

        Prioritizes each type of accession and returns the best record.

        Priority indices range from 7 (best, for a complete RefSeq record) and 6
        (complete GenBank record), to 5 (for complete RefSeq WGS) and all the 
        way to 3 (undetermined GenBank records).
        
        It returns a composite record with the nucleotide accession number for
        the "best" coding region, and the position and orientation of the CDS 
        within that accession, as well as the prioritization score obtained.
    """

    #Download IPG record for the specific ortholog
    records = None
    while records == None:
        try:
            records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
                                        rettype='ipg', retmode='xml'))
        except:
            print ("IPG record retrieval error: "+ortholog_acc)
            time.sleep(sleepy)
            pass
    
    #from the IPG record, retrieve all the genome accessions from all CDS
    #keeping only accession, location of start and strand, as well as 
    #priority score
    if 'ProteinList' in records["IPGReport"].keys():
        for idprotein in records["IPGReport"]["ProteinList"]:
            if 'CDSList' in idprotein.keys():
                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver']
                    if cds_acc.split(".")[0] == nucleotide:
                        cds_start = cds.attributes['start']
                        cds_stop = cds.attributes['stop']
                        cds_strand = cds.attributes['strand']
                        cds_org = cds.attributes['org'].replace(" ","_")
                        cds_scr = 0
                        
                        #create and append record
                        cds_rec = {'acc':cds_acc, 'start':cds_start, \
                                   'stop':cds_stop, 'strand':cds_strand,\
                                   'p_score':cds_scr, 'org':cds_org}
                        
                        return cds_rec
            else:
                return (None)
    #GenBank record that has no proper IPG record (yes, they exist;
    #see for instance: https://www.ncbi.nlm.nih.gov/protein/RJR51       119.1)
    #in these cases, there should be a CDS within the protein record that 
    #contains the information we want; priority should be lowest
    else:
#        TO BE IMPLEMENTED
#        records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
#                                            rettype='genbank', retmode='xml'))
        return (None)

    time.sleep(sleepy)
    

def GC_retrieval(nuc_id, start, stop):
    
    """ Download nuc data & compute the GC content
        Save organism info
    """
    results = []
    handle = Entrez.efetch(db="nucleotide", id = nuc_id, seq_start=start, seq_stop=stop, rettype = "gb")
    records = SeqIO.parse(handle, "gb")
    for record in records:
        organism = record.annotations["organism"]
        GC_result = float(GC(record.seq))
        break
    handle.close()
    
    results.append(organism)
    results.append(GC_result)
    return results

def GC_host_retrieval(organism,NCBI_file):

    """ Get the host GC content matching the organism to the NCBI assembly info
    """
    
    # save NCBI_file info into a dictionary
    organism_gb_assembly = {}
    with open(NCBI_file, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            if row:
                organism_gb_assembly[row["Organism Name"]] = round(float(row["GC%"]),2)
    
    # define the level of the match
    priority_list = ["Strain", "Species", "Genus"]
    match_key = ""
    # doing the match
    results = [None,None]
    for priority in priority_list:
        for genome_key in organism_gb_assembly.keys():
            if priority == "Strain":
                if genome_key == organism:
                    match_key = genome_key
                    if len(organism.split(" ")) > 2:
                        results[1] = priority
                    else:
                        results[1] = "Species"
                    break
            if priority == "Species":
                if " ".join(genome_key.split(" ")[0:2]) == " ".join(organism.split(" ")[0:2]):
                    match_key = genome_key
                    results[1] = priority
                    break
            if priority == "Genus":
                if " ".join(genome_key.split(" ")[0]) == " ".join(organism.split(" ")[0]):
                    match_key = genome_key
                    results[1] =  priority
                    break
        if match_key != "":
            results[0] = organism_gb_assembly[match_key]
            break     
    
    return results

############################################################################

# Load the configuration file
with open("test_input/2_ARG_prediction.json") as json_conf : 
    conf = json.load(json_conf)

# tell to who I am
Entrez.email = conf["email"]
Entrez.api_key = conf["api_key"]


## Get BLASTP hits & save CARD associated info
plasmid_hits = {}
errors = []
# open CARD metadata
with open(conf["card_db_csv"], 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        if row:
            # blast_search against desired plasmid database using each CARD accession
            card_acc = row["Protein Accession"]
            blastp_result = {}
            if card_acc:
                try:
                    protein_retrieval(card_acc, "card_temp.faa")
                    print ("BLASTp using "+card_acc+"...")
                    blastp_result = blast_search("card_temp.faa", conf["plasmid_db_name"], conf["blastp_eval"],1000,conf["blastp_query_cov"])
                except:
                    errors.append(card_acc)
                    pass
                
                if blastp_result:
                    # save blastp_result into plasmid_hits dict
                    for blastp_key in blastp_result.keys():
                        if blastp_key not in plasmid_hits.keys():
                            plasmid_hits[blastp_key] = {"best_hit":card_acc,"evalue":blastp_result[blastp_key],
                                                        "queries": []}
                        # check the best query for each blastp_result and overwrite the results
                        # the idea is to save the best hit to match CARD info
                        if blastp_result[blastp_key] < plasmid_hits[blastp_key]["evalue"]:
                            plasmid_hits[blastp_key]["evalue"] = blastp_result[blastp_key]
                            plasmid_hits[blastp_key]["best_hit"] = card_acc
                        # save this query
                        plasmid_hits[blastp_key]["queries"].append(card_acc)
                        plasmid_hits[blastp_key]["queries"] = list(set(plasmid_hits[blastp_key]["queries"]))
                        plasmid_hits[blastp_key]["gene_family"] = row["AMR Gene Family"]
                        plasmid_hits[blastp_key]["resistance_mechanism"] = row["Resistance Mechanism"]
                        plasmid_hits[blastp_key]["drug"] = row["Drug Class"]


## Get BLASTP hits info (GC, host GC, Pfam&COG annotation...)
for arg_gene in plasmid_hits.keys():
    print (arg_gene)
    # Pfam&COG annotation 
    # avoid pseudogenes
    if not arg_gene.split("__")[0].isdigit():
    	# save a temp file with aa seq
        protein_retrieval(arg_gene.split("__")[0],"temp_arg_gene.faa")
        # perform the hmmscan
        COG_annotation = HMM_annotation("temp_arg_gene.faa",conf["COG_db"], conf["hmmer_eval"])
        Pfam_annotation = HMM_annotation("temp_arg_gene.faa",conf["Pfam_db"], conf["hmmer_eval"])
        # save into the plasmid_hits dictionary
        plasmid_hits[arg_gene]["COG"] = COG_annotation
        plasmid_hits[arg_gene]["Pfam"] = Pfam_annotation
    
        ## GC analysis
        # get the genetic location of each arg_gene
        arg_gene_location = {}
        arg_gene_location = genome_record_retrieval(arg_gene.split("__")[0], arg_gene.split("__")[1], 0)
        if arg_gene_location:
            # compute the GC content & save info about the host
            try:
                GC_info = GC_retrieval(arg_gene_location["acc"], int(arg_gene_location["start"]), int(arg_gene_location["stop"]))
                plasmid_hits[arg_gene]["gc"] = round(GC_info[1],2)
                plasmid_hits[arg_gene]["host"] = GC_info[0]
                # get the host GC content trying to match the name of the host to the whole NCBI assembly info
                GC_host_info = GC_host_retrieval(plasmid_hits[arg_gene]["host"],conf["ncbi_assembly_file"])
                plasmid_hits[arg_gene]["gc_host"] = GC_host_info[0]
                plasmid_hits[arg_gene]["gc_host_match"] = GC_host_info[1]
            except:
                pass


## Save final results into a json file
with open(conf["output_json"], "w") as outfile: 
    json.dump(plasmid_hits, outfile)


## Print download errors if exist
if len(errors) >= 1:
    out_handle_error = open("errors.txt", "w")
    for error in errors:
        out_handle_error.write("Error trying to download "+error+"\n")
    out_handle_error.close()
