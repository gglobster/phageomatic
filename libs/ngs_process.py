from datetime import datetime
from common import ensure_dir
from config import directories as dirs, trim_param as tpm
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def trim_illumina(dataset):
    """Trim paired-end reads based on Illumina's PHRED quality scores.

    Trimmed reads are binned by size. Reads are loaded and processed as
    plain strings to avoid extra overhead from SeqIO. Modified from:
    http://news.open-bio.org/news/2010/04/illumina-q2-trim-fastq/

    """
    # Prep data source file
    source_file = dirs['ori_data_dir']+dataset['source_file']
    # Prep bin directories & files
    fn_root = dirs['trim_dir']+dataset['nickname']
    long_fname = fn_root+'_L.txt'
    short_fname = fn_root+'_S.txt'
    reject_fname = fn_root+'_R.txt'
    # Save filenames for next step (not saving rejected reads)
    dataset['trim_files'].append({'type': '_L', 'name': long_fname})
    dataset['trim_files'].append({'type': '_S', 'name': short_fname})
    # Unpack clipping sizes
    clip_start = tpm['clip_start']
    clip_end = tpm['clip_end']
    # Signal the process start
    print "-- Trimming", dataset['nickname'], "--"
    print datetime.now()
    # Open files for writing
    L_trim_out = open(long_fname, "w")
    S_trim_out = open(short_fname, "w")
    R_out = open(reject_fname, "w")
    # Initialize counters
    total_count = 0
    long_bin_count = 0
    short_bin_count = 0
    rejected_count = 0
    # Cycle through reads
    for titles, seqs, quals in FastqGGIterator(open(source_file)) :
        F_title, R_title = titles
        F_seq, R_seq = seqs
        F_qual, R_qual = quals
        #Find the location of the first "B" (PHRED quality 2)
        F_trim = F_qual.find("B")
        R_trim = R_qual.find("B")
        if F_trim == -1:
            F_trim = R_trim
        if R_trim == -1:
            R_trim = F_trim
        trim = min(F_trim, R_trim)
        # no B found, no need to trim, just clip
        if trim == -1:
            trim = tpm['full']-clip_end
            L_trim_out.write("@%s\n%s\n+\n%s\n" % (F_title,
                                                   F_seq[clip_start:trim],
                                                   F_qual[clip_start:trim]))
            L_trim_out.write("@%s\n%s\n+\n%s\n" % (R_title,
                                                   R_seq[clip_start:trim],
                                                   R_qual[clip_start:trim]))
            long_bin_count +=1
        # B found below min size, no point keeping the read
        elif trim < tpm['min_threshold']:
            R_out.write("@%s\n%s\n+\n%s\n" % (F_title, F_seq, F_qual))
            R_out.write("@%s\n%s\n+\n%s\n" % (R_title, R_seq, R_qual))
            rejected_count +=1
        # B found after min size, trim to minimum size
        elif trim < tpm['full']:
            trim = tpm['min_threshold']
            S_trim_out.write("@%s\n%s\n+\n%s\n" % (F_title,
                                                   F_seq[clip_start:trim],
                                                   F_qual[clip_start:trim]))
            S_trim_out.write("@%s\n%s\n+\n%s\n" % (R_title,
                                                   R_seq[clip_start:trim],
                                                   R_qual[clip_start:trim]))
            short_bin_count +=1
        # increment counter
        total_count +=1
        # report on the progress
        if total_count%200000==0:
            print "\t", total_count, "read pairs processed"
    # Close bin files
    L_trim_out.close()
    S_trim_out.close()
    R_out.close()
    # Summarize count results
    counts_dict = {'long_trim': long_bin_count,
                   'rejected': rejected_count,
                   'short_trim': short_bin_count,
                   'total': total_count
                   }
    print "-- Finished --"
    print datetime.now()
    return counts_dict

def simple_q2a(dataset, trim_file):
    """."""
    # Identify the genome
    nickname = dataset['nickname']
    # Identify the trim file type
    ttype = trim_file['type']
    # Identify the source file
    source_file = trim_file['name']
    # Prep output files
    dir_root = dirs['mft_dir']+nickname+"/"
    ensure_dir(dir_root)
    out_file = dir_root+nickname+ttype+'.fas'
    track_file = dir_root+nickname+ttype+'_track.txt'
    # Save filenames for later reference
    dataset['mft_files'].append({'type': ttype,
                                 'name': out_file,
                                 'track': track_file})
    # Signal the process start
    print "-- Converting ", nickname+ttype, "to multifasta --"
    print datetime.now()
    # Set up iterator
    multifasta = open(out_file, 'w')
    tracker = open(track_file, 'w')
    read_count = 0
    for title, seq, qual in FastqGeneralIterator(open(source_file)) :
        read_count +=1
        id_string = nickname+"_"+str(read_count)
        multifasta.write(">"+id_string+"\n"+seq+"\n")
        tracker.write(title+"\t"+id_string+"\n")
        if read_count%100000==0:
            print "\t"+str(read_count), "reads processed"
    multifasta.close()
    tracker.close()
    return read_count

def FastqGGIterator(handle):
    """Iterate over 2 Fastq records as string tuples.

    Modified from Bio.SeqIO.QualityIO import FastqGeneralIterator to return
    two reads at a time.

    """
    handle_readline = handle.readline
    while True:
        line = handle_readline()
        if line == "" : return
        if line[0] == "@":
            break
    while True:
        title_lines = []
        seq_strings = []
        quality_strings = []
        count = 0
        while count < 2 :
            if line[0] != "@":
                raise ValueError("Bad formatting of record start lines.")
            title_line = line[1:].rstrip()
            title_lines.append(title_line)
            seq_string = handle_readline().rstrip()
            seq_strings.append(seq_string)
            while True:
                line = handle_readline()
                if not line:
                    raise ValueError("End of file without quality info.")
                if line[0] == "+":
                    second_title = line[1:].rstrip()
                    if second_title and second_title != title_line:
                        raise ValueError("Seq and qual captions differ.")
                    break
                seq_string += line.rstrip() #removes trailing newlines
            if " " in seq_string or "\t" in seq_string:
                raise ValueError("Whitespace not allowed in the sequence.")
            seq_len = len(seq_string)
            quality_string = handle_readline().rstrip()
            quality_strings.append(quality_string)
            while True:
                line = handle_readline()
                if not line : break #end of file
                if line[0] == "@":
                    if len(quality_string) >= seq_len:
                        break
                quality_string += line.rstrip()
            if seq_len != len(quality_string):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i)." \
                                 % (title_line, seq_len, len(quality_string)))
            count +=1
        yield (title_lines, seq_strings, quality_strings)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"