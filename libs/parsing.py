import os, re, numpy
from loaders import read_array
from config import directories as dirs, bin_types, blast_dtypes, \
    chop_param as cpm, references
import matplotlib.pyplot as pylot
from common import ensure_dir
from datetime import datetime
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Blast import NCBIXML

def glomp_blast_out(dataset, ref_nick):
    """Consolidate Blast output files."""
    # Identify the genome
    nickname = dataset['nickname']
    # Determine the input file root
    root_dir = dirs['blast_out_dir']+nickname+"/"+ref_nick+"/"
    file_root = root_dir+nickname
    # Signal process start
    print "-- Consolidating B_out for", nickname, "against", ref_nick, "--"
    print datetime.now()
    # Cycle through bin types
    series_index = 0
    averages = []       # for comparing series later
    binned_pos = []
    for bin_type in bin_types:
        index = 1
        bin_arrays =[]
        while os.path.isfile(file_root+bin_type+"_"+str(index)+"_blast.out"):
            infile = file_root+bin_type+"_"+str(index)+"_blast.out"
            rec_array = read_array(infile, blast_dtypes)
            if len(rec_array) > 0:
                bin_arrays.append(rec_array)
            index +=1
        print "\t\t"+str(len(bin_arrays)), "arrays for", \
        nickname+bin_type, "series"
        if len(bin_arrays) > 0:
            series = numpy.hstack(bin_arrays)
        else:
            series = []
        print "\t\t"+str(len(series)), "total records in", \
        nickname+bin_type, "series"
        # Save to file
        cons_outfile = file_root+bin_type+"_cons_out.npy"
        numpy.save(cons_outfile, series)
        # Evaluate match positions on reference
        positions = []
        match_read = []
        for row in series:
            # collect match read info while we're at it
            # use regex to extract query index
            query_pattern = re.compile(r'\w*_(\d*)')
            query_match = query_pattern.match(row[0])
            query_index = int(query_match.group(1))
            match_read.append(query_index)
            # use regex to extract ref coords
            ref_pattern = re.compile(r'\w*_\d*_(\d*)')
            ref_match = ref_pattern.match(row[1])
            ref_pos = int(ref_match.group(1))
            pos_scaled = ref_pos/cpm['size']   # adjust to db segment length
            positions.append(pos_scaled)
        # uniquify the match read array
        unique_matches = numpy.unique(match_read)
        print "\t"+str(len(unique_matches)), "unique matches for", bin_type
        # write to file for future use
        match_dir_root = dirs['match_dir']+nickname+"/"+ref_nick+"/"
        ensure_dir(match_dir_root)
        match_outfile = match_dir_root+nickname+bin_type+"_match.npy"
        numpy.save(match_outfile, unique_matches)
        # now count ocurrences per position
        pos_np = numpy.array(positions)
        binned = numpy.bincount(pos_np)
        binned_pos.append(binned)
        pos_count_average = numpy.average(binned)
        averages.append((pos_count_average, series_index))
        series_index +=1
    # compare series
    averages.sort()
    averages.reverse()
    order_indices = []
    for pair in averages:
        order_indices.append(pair[1])
    # identify reference
    ref_name = [reference['full_name'] for reference in references if
                reference['nickname'] is ref_nick]
    # prep directory & file
    fig_root = dirs['reports_dir']+"match_figs/"
    fig_file = fig_root+nickname+"_"+ref_nick+".png"
    ensure_dir(fig_root)
    # generate a figure
    pylot.autoscale(enable=True, axis='both', tight=True)
    pylot.xlabel('Position on the chromosome (/'+str(cpm['size'])+')')
    pylot.ylabel('Number of matches (includes multiples)')
    pylot.title(nickname+' matches to '+ref_name)
    pylot.grid(True)
    for index in order_indices:
        label_root = nickname+bin_types[index]
        label_str = label_root+" ("+str(numpy.sum(binned_pos[index]))+")"
        pylot.plot(binned_pos[index], label=label_str)
    pylot.legend(loc=1)
    pylot.savefig(fig_file, dpi=None, facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format=None)
    pylot.clf()
    print "\t"+str(series_index), "series consolidated and parsed"
    print "-- Done, see plot --"
    print datetime.now()
    return "OK"

def glomp_good_reads(dataset, bin_type):
    """Use matching reads lists to glomp good reads."""
    # Identify dataset
    nickname = dataset['nickname']
    # Identify input directory and filenames root
    match_dir_root = dirs['match_dir']+nickname+"/"
    # Signal process start
    print "-- Glomping matches against all references for", nickname, "--"
    print datetime.now()
    filter_files = []
    rescue_files = []
    for ref_nick in dataset['ref_nicks']:
        infile = match_dir_root+ref_nick+"/"+nickname+bin_type+"_match.npy"
        ref_type = [reference['type'] for reference in references if
                    reference['nickname'] is ref_nick]
        if ref_type[0] is 'filter':
            filter_files.append({'ref_nick': ref_nick, 'matches': infile})
        elif ref_type[0] is 'rescue':
            rescue_files.append({'ref_nick': ref_nick, 'matches': infile})
    # process filter files
    print "\tprocessing filter references", [filter_file['ref_nick'] for
                                             filter_file in filter_files]
    filter_arrays = []
    for file in filter_files:
        data_array = numpy.load(file['matches'])
        filter_arrays.append(data_array)
    array_index = 0
    filter_IRA = filter_arrays[array_index]
    filter_URA = filter_arrays[array_index]
    while array_index < len(filter_arrays)-1:
        array_index +=1
        filter_IRA = numpy.intersect1d(filter_IRA, filter_arrays[array_index])
        filter_URA = numpy.union1d(filter_URA, filter_arrays[array_index])
    print "\t\t"+str(len(filter_IRA)), "present in all filter references"
    print "\t\t"+str(len(filter_URA)), "matching reads all together (union)"
    # process rescue files
    print "\tprocessing rescue references", [rescue_file['ref_nick'] for
                                             rescue_file in rescue_files]
    rescue_arrays = []
    for file in rescue_files:
        data_array = numpy.load(file['matches'])
        rescue_arrays.append(data_array)
    array_index = 0
    rescue_IRA = rescue_arrays[array_index]
    rescue_URA = rescue_arrays[array_index]
    while array_index < len(rescue_arrays)-1:
        array_index +=1
        rescue_IRA = numpy.intersect1d(rescue_IRA, rescue_arrays[array_index])
        rescue_URA = numpy.union1d(rescue_URA, rescue_arrays[array_index])
    print "\t\t"+str(len(rescue_IRA)), "present in all rescue references"
    print "\t\t"+str(len(rescue_URA)), "matching reads all together (union)"
    # prepare for masking
    print "\tpreparing selection masks"
    q2a_file = dirs['mft_dir']+nickname+"/"+nickname+bin_type+"_track.txt"
    dtype = numpy.dtype([('title', 'S50'), ('bincode', 'S15')])
    pair_array = read_array(q2a_file, dtype, separator='\t')
    # create masking arrays
    mask = numpy.zeros(len(pair_array), bool)
    mask = numpy.invert(mask)
    # filter out baddies - flip to False
    for item in filter_URA:
        mask[item-1] = False    # False means reject
        if item%2==0:
            mask[item-2] = False    # even numbers are /2, flip previous
        else:
            mask[item] = False      # odd numbers are /1, flip next

    # rescue goodies - flip to True
    for item in rescue_URA:
        mask[item-1] = True     # True means accept
        if item%2==0:
            mask[item-2] = True     # even numbers are /2, flip previous
        else:
            mask[item] = True       # odd numbers are /1, flip next
    # save mask to file (where True means keep, False means reject)
    mask_file = dirs['trim_dir']+nickname+bin_type+"_mask.npy"
    numpy.save(mask_file, mask)
    # separate the two sets and write to file
    bin_accept_titles = pair_array[mask]
    accept_file = dirs['select_dir']+nickname+bin_type+"_accept.npy"
    numpy.save(accept_file, bin_accept_titles)
    inv_mask = numpy.invert(mask)
    bin_reject_titles = pair_array[inv_mask]
    reject_file = dirs['select_dir']+nickname+bin_type+"_reject.npy"
    numpy.save(reject_file, bin_reject_titles)
    print "\t\t"+str(len(bin_accept_titles)), "reads to accept"
    print "\t\t"+str(len(bin_reject_titles)), "reads to reject"
    print "-- Done! --"
    print datetime.now()

def extract_read_sets(dataset, bin_type):
    """Use masks to extract read sets."""
    # Identify dataset
    nickname = dataset['nickname']
    # Identify input & output files
    source_file = dirs['trim_dir']+nickname+bin_type+".txt"
    mask_file = dirs['trim_dir']+nickname+bin_type+"_mask.npy"
    out_file = dirs['assembly_dir']+nickname+bin_type+"_ok.txt"
    # Signal process start
    print "-- Extracting reads for", nickname+bin_type, "--"
    print datetime.now()
    out_handle = open(out_file, 'w')
    counter = 0
    tru_count = 0
    mask = numpy.load(mask_file)
    for title, seq, qual in FastqGeneralIterator(open(source_file)) :
        if mask[counter] == True:
            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            tru_count +=1
        else:
            pass
        counter +=1
        if counter%100000==0:
            print "\t"+str(tru_count), "reads selected of", str(counter), \
            "reads processed"
    out_handle.close()
    print "-- Done! --"
    print datetime.now()

def collect_cogs(blast_out):
    results = {}
    blast_records = NCBIXML.parse(open(blast_out))
    for record in blast_records:
        if record.alignments: # ignores searches with no hits
            top_hit = record.alignments[0]
            hit_name = top_hit.title
            results[record.query_id] = hit_name
        else:
            results[record.query_id] = 'no match'
    return results
