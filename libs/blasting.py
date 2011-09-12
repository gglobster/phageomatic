import os, subprocess
from config import blast_prefs, directories as dirs
from common import ensure_dir
from datetime import datetime
from Bio.Blast.Applications import NcbiblastnCommandline, \
    NcbirpsblastCommandline
from Bio.Blast import NCBIWWW

def make_blastDB(name, infile, db_type):
    """Make BLAST database from FASTA input file."""
    cline = "makeblastdb -in "+ infile +" -dbtype "+ db_type +" -title "+  \
            infile +" -out "+ name +" -parse_seqids"
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def local_blastn_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastn against local database."""
    cline = NcbiblastnCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_rpsblast_2file(query_file, dbfile_path, outfile, prefs):
    """Perform RPS Blast against local database."""
    cline = NcbirpsblastCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=5) # must output XML!
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def remote_blastp_2file(query_string, database, outfile, evalue):
    """Perform blastp against remote database."""
    result_handle = NCBIWWW.qblast('blastp',
                                   database,
                                   query_string,
                                   expect=evalue)
    save_file = open(outfile, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

def mux_batch_blast(dataset, bin_type):
    """Send batch jobs to Blast. Muxes to multiple reference DBs."""
    # Identify the genome
    nickname = dataset['nickname']
    # Determine the input file root
    root_dir = dirs['parcels_dir']+nickname+"/"
    file_root = root_dir+nickname+bin_type
    # Identify the references to blast against
    ref_nicks = dataset['ref_nicks']
    for ref_nick in ref_nicks:
        # Identify Blast DB
        db_path = dirs['blast_db_dir']+ref_nick
        # Prep output directory
        out_dir = dirs['blast_out_dir']+nickname+"/"+ref_nick+"/"
        ensure_dir(out_dir)
        # Signal process start
        print "--- Blasting", nickname+bin_type, "against", ref_nick, "---"
        print datetime.now()
        index = 1
        while os.path.isfile(file_root+"_"+str(index)+".fas"):
            query_file = file_root+"_"+str(index)+".fas"
            outfile = out_dir+nickname+bin_type+"_"+str(index)+"_blast.out"
            print "\tblasting", query_file
            local_blastn_2file(query_file, db_path, outfile, blast_prefs)
            index +=1
        print "--- Finished BLAST run ---"
        print datetime.now()
        print index, "parcel files blasted"
    return "OK"


