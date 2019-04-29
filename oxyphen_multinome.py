from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.feature_selection import RFE

from Bio.ExPASy import Enzyme
import numpy as np
import pandas as pd
import os, glob

CONFIG_FILE = open("SETTINGS", "r").read().splitlines()

def read_config():
    multinome_folder=""
    for line in CONFIG_FILE:
        if line.startswith("BLAST_PATH"):
            blast_path = line.split("=")[1]
        if line.startswith("NUM_THREADS"):
            num_threads = float(line.split("=")[1])
        if line.startswith("PROTEOMES_FOLDER") and line.split("=")[1]:
            ### enter multinome mode
            multinome_folder = line.split("=")[1]

    return blast_path, num_threads, multinome_folder

def do_oxyphen(proteome, output_filename, ec_classes_file):

    '''
    Read and parse enzyme.dat file
    '''
    input_name = "DATA/enzyme.dat"
    output_name = 'DATA/ec_uniprot.tsv'

    ### program ###
    handle = open(input_name)
    records = Enzyme.parse(handle)

    out = dict() #dict of dicts, first key: EC number, second key: field
    transferred = dict() #dict of lists
    for record in records:
        if 'Transferred entry:' in record['DE']:
            record['DE'] = record['DE'].rstrip('.') #remove period
            record['DE'] = record['DE'].replace('Transferred entry:',' ') #remove title
            record['DE'] = record['DE'].replace(',',' ') #remove commas
            record['DE'] = record['DE'].replace('and',' ') #remove and
            point_to = record['DE'].split()
            transferred[record['ID']] = point_to
        else:
            out[record['ID']] = dict()
            out[record['ID']]['uniprot'] = ' '.join([x[0] for x in record['DR']])
            out[record['ID']]['description'] = record['DE']
            out[record['ID']]['transferred'] = False

    for id in transferred:
        out[id] = dict()
        out[id]['uniprot'] = ' '.join([out[x]['uniprot'] for x in transferred[id]])
        out[id]['description'] = 'Transferred entry: ' + ' '.join(transferred[id])
        out[id]['transferred'] = True
    df = pd.DataFrame.from_dict(out, orient='index')
    df.index.name = 'EC'
    df.to_csv(output_name, sep='\t')

    '''
    Take a subset of ecs of interest
    '''

    oxidases = tuple(open("DATA/oxygen_ecclasses", "r").read().splitlines())

    infile = open("DATA/ec_uniprot.tsv", "r").readlines()
    outfile = open("DATA/ec_uniprot_oxidases.tsv", "w")

    for line in infile:
        if line.startswith("EC"):
            outfile.write(line)
        elif line.startswith(oxidases):
            outfile.write(line)

    outfile.close()

    '''
    write a file with one uniprot ID per line, containing all of the
    uniprot IDs mentioned in uniprot column of the input file

    Ignore EC numbers that have been transferred
    '''

    input = "DATA/ec_uniprot_oxidases.tsv"
    output = "DATA/uniprot_ids.txt"

    df = pd.read_table(input)
    df.dropna(subset=['uniprot'], inplace=True) #ignore EC numbers with no uniprot ids associated

    df = df[df.transferred == False] #ignore EC numbers that are obsolete due to transfer

    unique_uniprot = set(" ".join(df.uniprot.values).split(" "))

    with open(output, "w") as outfile:
        for id in unique_uniprot:
            outfile.write(id + "\n")
    outfile.close()

    '''
    Make blastdb out of the swissprot subset
    '''

    blast_path, num_threads, multinome_folder = read_config()

    os.system("%s -in DATA/sprot_subset.fasta -dbtype prot -out DATA/sprot_subset -hash_index" % (os.path.join(blast_path, "makeblastdb")))

    '''
    Blast our pre-selected proteomes against the uniprot subset
    '''
    print "Performing Blast searches against oxygen-utilizing database..."
    os.system("%s -max_target_seqs 1 -outfmt '6 qseqid sseqid pident evalue qcovs' -query %s -db DATA/sprot_subset -out DATA/new_sequences_sprot_enzyme.tab -num_threads %d" % (os.path.join(blast_path, "blastp"), proteome, num_threads) )

    '''
    Filter Blast output.
    '''
    evalue = 10e-3
    identity = 40.0
    coverage = 40.0

    print "Filtering Blast output: evalue",evalue, " identity", identity, " coverage", coverage
    hits_table_file_name = "DATA/new_sequences_sprot_enzyme.tab"
    hits_table_file_name_filtered_out = open("DATA/new_sequences_sprot_enzyme_filtered.tab", "w")


    hits_table_file_name_filtered_out.write("\t".join(["hit","subject","id","len","eval","cov"])+"\n")


    for line in open(hits_table_file_name, "r").read().splitlines():
        if line.startswith("#"): continue

        query, target, ident, eval, cover = line.split("\t")
        eval = float(eval)
        ident = float(ident)
        cover = float(cover)

        if eval <= evalue and ident >= identity and cover >= coverage:
            hits_table_file_name_filtered_out.write(line+"\n")

    hits_table_file_name_filtered_out.close()

    hits_table_file_name_filtered = "DATA/new_sequences_sprot_enzyme_filtered.tab"
    enzyme_table_file_name = 'DATA/ec_uniprot_oxidases.tsv'

    hits = pd.read_csv(hits_table_file_name_filtered, sep="\t", header=0)
    enzyme = pd.read_csv(enzyme_table_file_name, sep="\t", header=0)

    hits.fillna('', inplace=True)  #replace empty values with blank spaces
    enzyme.fillna('', inplace=True)

    enzyme = enzyme[enzyme.transferred == False] #drop transferred EC numbers

    hits.subject = hits.subject.str[3:9] #take just the uniprot ID from the name

    def get_ecs(uniprot):
        if uniprot == '': #ignore invalid uniprot ids
            return ''
        else:
            return ' '.join(enzyme.EC[enzyme.uniprot.str.contains(uniprot)].values)

    hits['EC'] = hits.subject.apply(get_ecs)

    output_file_name = output_filename
    hits.to_csv(output_file_name, sep="\t", index=False)

    ### read final mapping output

    mapping_out = open(output_file_name, "r").read().splitlines()
    ecs_dict = {}

    for line in mapping_out[1:]:
        splitted = line.split("\t")
        ecs = splitted[-1]

        for ec in ecs.split():
            if ec not in ecs_dict:
                ecs_dict[ec] = []
            ecs_dict[ec].append(splitted[0])

    print "\n\n"
    print len(ecs_dict), "oxygen-utilizing enzymes were found from classes", ecs_dict.keys()

    ec_out = open(ec_classes_file, "w")
    ec_out.write("\t".join(ecs_dict.keys()))

    ec_out.close()
    #print "Detailed mapping can be found in OUTPUT/oxygen_utilizing_annot.tsv file"
    #print "Executing SVM classifier..."

    infile = open("DATA/model_svm", "r").read().splitlines()

    classifier_input = []
    classes = []
    ec_classes = []

    for line in infile:

    	if line.startswith("@attribute") and "class" not in line:
    		ec_classes.append(line.split()[1].replace("'",""))

def do_for_all():
    """
    Execute oxyphen for all proteomes in the input folder
    """
    blast_path, num_threads, multinome_folder = read_config()
    proteomes = glob.glob(os.path.join(multinome_folder,"*"))

    print "\n\nPROTEOMES IN YOUR PROTEOMES_FOLDER DIRECTORY:\n", "\n".join(proteomes)

    for proteome in proteomes:
        fname = os.path.splitext(os.path.basename(proteome))[0]
        output_filename = os.path.join("OUTPUT",fname+"_oxygen_utilizing_annot.tsv")

        ec_classes_file = os.path.join("OUTPUT",fname+"_EC_CLASSES.txt")

        print "\n\nRUNNING OXYPHEN FOR PROTEOME %s" % proteome
        do_oxyphen(proteome, output_filename, ec_classes_file)

        print "\n\nOUTPUT MAPPING FILE FOR THIS PROTEOME CAN BE FOUND IN %s" % output_filename
        print "LIST OF EC CLASSES FOR THIS PROTEOME CAN BE FOUND IN %s" % ec_classes_file

if __name__ == '__main__':
    do_for_all()

    print "YOUR OUTPUT FILES CAN BE FOUND IN OUTPUT/ FOLDER"
