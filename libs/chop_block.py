from config import directories as dirs, chop_param
from datetime import datetime
from common import ensure_dir
from loaders import load_agnostic, batch_iterator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from blasting import make_blastDB

def chop_genome(genome):
    """Chop genome to size.

    Outputs a single multifasta file with all the sequence records resulting
    from the chopping process. Also generates a Blast database.

    """
    # Identify the genome
    nickname = genome['nickname']
    # Unpack chopping parameters
    seg_size = chop_param['size']
    # Prep data source file
    source_file = dirs['reference_dir']+genome['source_file']
    # Prep output file
    out_file_root = dirs['chopped_dir']+nickname
    output_handle = out_file_root+'_chop_'+str(seg_size)+'.fas'
    # Save filename for later reference
    genome['chopped_file'] = output_handle
    # Signal the process start
    print "-- Chopping", nickname, "to size", seg_size, "--"
    print datetime.now()
    # Load the genome into memory
    g_data = load_agnostic(source_file)
    genome_length = len(g_data[0].seq)
    # Create list of segment coordinates
    seg_coords = []
    index = 0
    while index < genome_length:
        incr_index = index+seg_size
        pair = (index, incr_index)
        seg_coords.append(pair)
        index = incr_index
    print "\t"+str(len(seg_coords)), "segments created"
    # Extract the segment sequences, writing to file as we go
    records = []
    for coord_pair in seg_coords:
        start = coord_pair[0]
        stop = coord_pair[1]
        seg_seq = g_data[0].seq[start:stop]
        seg_record = SeqRecord(seg_seq,
                               id=nickname+"_"+str(start)+"_"+str(stop))
        records.append(seg_record)
    SeqIO.write(records, output_handle, 'fasta')
    print "\twritten to", output_handle
    # Make a BLAST DB from the output file
    db_name = dirs['blast_db_dir']+nickname#+str(seg_size)
    make_blastDB(db_name, output_handle, "nucl")
    print "\tDB", db_name, "created"
    print "-- Finished --"
    print datetime.now()
    return len(seg_coords)

def chop_multifasta(dataset, mft_file):
    """Split a master file into smaller multifasta files.

    This is useful for making reasonably-sized batch BLAST jobs.
    Iterator function adapted from http://biopython.org/wiki/Split_large_file

    """
    # Identify the genome
    nickname = dataset['nickname']
    # Identify the trim file type
    ttype = mft_file['type']
    # Prep output file
    dir_root = dirs['parcels_dir']+nickname+"/"
    out_file_root = dir_root+nickname+ttype
    ensure_dir(dir_root)
    # Save filenames for later reference
    dataset['parcel_files'] = {'root': out_file_root, 'suffix': '.fas'}
    # Unpack chopping parameters
    parceln = chop_param['parceln']
    # Signal the process start
    print "-- Splitting ", nickname+ttype, "into batches of", parceln, "reads --"
    print datetime.now()
    # Set up iterator function
    record_iter = SeqIO.parse(open(mft_file['name']),"fasta")
    parcel_files = []
    for i, batch in enumerate(batch_iterator(record_iter, parceln)) :
        filename = out_file_root+"_%i.fas" % (i+1)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        parcel_files.append(filename)
        print "\twrote %i records to %s" % (count, filename)
    print "-- Finished --"
    print datetime.now()
    return len(parcel_files)
