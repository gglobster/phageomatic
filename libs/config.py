import numpy

### Config file ###

## Data identification

# datasets
datasets = [{'source_file': 'reads.txt',
             'full_name': 'Phage from',
             'nickname': '',
             'f_nick': '',
             'ref_nicks': ('')}]
# references
references = [{'source_file': 'g.fas',
               'full_name': '',
               'nickname': '',
               'type': 'filter'},
              {'source_file': 'phages.fas', # already chopped
               'full_name': 'phages',
               'nickname': 'phages',
               'type': 'rescue'}]


## Directory structure

# project root directory
root_dir = 'data/run/'

directories = {
'ori_data_dir': root_dir+'original/',      # original data
'trim_dir': root_dir+'trimmed/',           # trimmed reads
'chopped_dir': root_dir+'chopped/',        # multifasta ref chop
'mft_dir': root_dir+'mf_trim/',            # multifasta trim reads
'parcels_dir': root_dir+'parcels/',        # parceled multifasta
'reference_dir': root_dir+'reference/',    # reference genomes
'blast_db_dir': root_dir+'blast_db/',      # blast databases
'blast_out_dir': root_dir+'blast_out/',    # blast outputs
'match_dir': root_dir+'matches/',          # match sets
'assembly_dir': root_dir+'assembly/',      # assembly
'annot_dir': root_dir+'annotation/',       # annotation
'reports_dir': root_dir+'reports/',        # reports
'select_dir': root_dir+'selection/'        # selection
}

## Parameters

# binning parameters
bin_types = ('_L', '_S')#, '_R')
# trimming parameters
trim_param = {'full': 101,
              'min_threshold': 50,
              'clip_start': 5,
              'clip_end': 21}
# reference chopping parameters
chop_param = {'size': 100, 'parceln': 10000}
# Blast parameters
blast_prefs = {'evalue': 0.001,
               'outfmt_pref': 6,
               'score': 100,
               'length': 60}
# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint8'),
                           ('q_end', 'uint8'),
                           ('r_start', 'uint8'),
                           ('r_end', 'uint8'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])
##