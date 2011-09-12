from sys import argv, exit
from libs.common import ensure_dir
from libs.ngs_process import simple_q2a, trim_illumina
from libs.chop_block import chop_genome, chop_multifasta
from libs.blasting import mux_batch_blast
from libs.parsing import glomp_blast_out, glomp_good_reads, extract_read_sets
from libs.config import datasets, references, directories, bin_types
from libs.contigs_process import batch_contig_annot


if len(argv) > 1 | argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python main_script.py [step#]\n", \
          "For better results, run from within iPython."
    exit()

print "##################################################\n", \
      "### PhageOMatic v. 0.1                         ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) < 2:
    step = 0
else:
    step = int(argv[1])

if step is 0:
    ### STEP 0: Ensure that all base directories exist ###
    print "\n###", step, ". Setting up the work environment ###"
    for dir_name in directories.keys():
        ensure_dir(directories[dir_name])
    step +=1

if step is 1:
    ### STEP 1: Trim & bin reads (based on results of FastQC) ###
    print "\n###", step, ". Trim & bin, then split for batching ###"
    for dataset in datasets:
        dataset['trim_files'] = []
        bin_counts = trim_illumina(dataset)
        print bin_counts # TO LOG
        dataset['mft_files'] = []
        for trim_file in dataset['trim_files']:
            mft_count = simple_q2a(dataset, trim_file)
            #print mft_count # TO LOG
        for mft_file in dataset['mft_files']:
            print mft_file
            par_count = chop_multifasta(dataset, mft_file)
            #print par_count # TO LOG
    step +=1

if step is 2:
    ### STEP 2: Chop reference genomes and generate BLAST databases ###
    print "\n### ", step, ". Chop references and make Blast DBs ###"
    for ref in references:
        if ref['type'] == 'filter':
            rec_count = chop_genome(ref)
        elif ref['type'] == 'rescue':
            pass    # TODO: add chopNconcat function to pipeline
    step +=1

if step is 3:
    ### STEP 3: Blast read bins against reference databases ###
    print "\n### ", step, ". Blast read sets against reference databases ###"
    for dataset in datasets:
        for type in bin_types:
            status = mux_batch_blast(dataset, type)
    step +=1

if step is 4:
    ### STEP 4: Consolidate BLAST output files and output unique reads ###
    print "\n### ", step, ". Consolidate and analyze Blast output ###"
    for dataset in datasets:
        for ref_genome in dataset['ref_nicks']:
            # consolidate blast_out files and parse match positions
            status = glomp_blast_out(dataset, ref_genome)
        # combine match reads from all references
        for bin_type in bin_types:
            status = glomp_good_reads(dataset, bin_type)
    step +=1

if step is 5:
    ### STEP 5: Extract selected reads using boolean mask ###
    print "\n### ", step, ". Extract selected reads ###"
    for dataset in datasets:
        for type in bin_types:
            extract_read_sets(dataset, type)
    step +=1

if step is 6:
    ### STEP 6: Assemble with Velvet ###
    print "\n### ", step, ". Assemble with Velvet ###"
    for dataset in datasets:
        # prep data for assembly
        # run assembler
        pass # not yet ready for integration into the pipeline
    step +=1

if step is 7:
    ### STEP 7: Annotate size-appropriate contigs ###
    print "\n### ", step, ". Extract and annotate contigs ###"
    for dataset in datasets:
        wtf = batch_contig_annot(dataset)
    step +=1

if step > 7:
    print "Error: step number is too high!"