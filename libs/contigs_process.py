import re, subprocess
from loaders import load_multifasta, load_fasta, load_genbank
from writers import write_fasta, write_genbank
from common import ensure_dir
from os import path
from blasting import local_rpsblast_2file
from parsing import collect_cogs
from config import directories as dirs, blast_prefs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from graphics import ContigDraw

def train_prodigal(seq_file, training_file):
    """Train Prodigal on the entire dataset."""
    cline = "prodigal -i "+seq_file+" -t "+training_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def run_prodigal(in_file, an_gbk, an_aa, trn_file):
    """Annotate sequence records individually using Prodigal."""
    cline = "prodigal -i "+in_file+" -o "+an_gbk+" -a "+an_aa+" -t "+trn_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def batch_contig_annot(dataset):
    """Extract and annotate contigs."""
    # identify dataset contig file
    contigs_file = dirs['assembly_dir']+dataset['f_nick']+'/'+'contigs.fa'
    # locate the COG database
    cog_db = dirs['blast_db_dir']+'Cog_LE/Cog'
    # make the training file
    training_file = dirs['annot_dir']+dataset['f_nick']+'/'+'contigs.trn'
    #train_prodigal(contigs_file, training_file)
    # set output dirs
    fas_out_dir = dirs['annot_dir']+dataset['f_nick']+'/fasta/'
    gbk_out_dir = dirs['annot_dir']+dataset['f_nick']+'/predict/'
    aa_out_dir = dirs['annot_dir']+dataset['f_nick']+'/aa/'
    blast_out_dir = dirs['annot_dir']+dataset['f_nick']+'/rpsblast/'
    solid_out_dir = dirs['annot_dir']+dataset['f_nick']+'/genbank/'
    maps_out_dir = dirs['annot_dir']+dataset['f_nick']+'/maps/'
    ensure_dir(fas_out_dir)
    ensure_dir(gbk_out_dir)
    ensure_dir(aa_out_dir)
    ensure_dir(blast_out_dir)
    ensure_dir(solid_out_dir)
    # set phage hit collector
    contig_hits = {}
    sp_hit_list = dirs['annot_dir']+dataset['f_nick']+'/'\
                                   +dataset['f_nick']+'_kw_hits.html'
    all_hit_list = dirs['annot_dir']+dataset['f_nick']+'/'\
                                    +dataset['f_nick']+'_all_hits.html'
    sp_hit_list_handle = open(sp_hit_list, 'w')
    all_hit_list_handle = open(all_hit_list, 'w')
    sp_hit_list_handle.write("<ul>")
    all_hit_list_handle.write("<ul>")
    # load all contigs
    contigs_list = load_multifasta(contigs_file)
    # cycle through contigs
    ctg_count = 0
    gene_count = 0
    for contig in contigs_list:
        ctg_count +=1
        # use regex to acquire relevant record ID info
        pattern = re.compile(r'NODE_(\d*)_length_(\d*)_cov_(\d*)')
        match = pattern.match(contig.id)
        nick = match.group(1)+'_'+match.group(2)+'_'+match.group(3)
        contig.id = nick
        fasta_out = fas_out_dir+nick+'.fas'
        # write record to file
        write_fasta(fasta_out, contig)
        # create contig entry in dict
        contig_hits[nick] = []
        # run the annotation
        annot_gbk = gbk_out_dir+nick+'.gbk'
        annot_aa = aa_out_dir+nick+'.fas'
        #run_prodigal(fasta_out, annot_gbk, annot_aa, training_file)
        # blast the amino acids against COG
        print '\tblasting', dataset['f_nick'], nick
        blast_out = blast_out_dir+nick+'.xml'
        if path.isfile(blast_out):
            print "\t\talready blasted"
        else:
            local_rpsblast_2file(annot_aa, cog_db, blast_out, blast_prefs)
        # collect best hits
        rec_cogs = collect_cogs(blast_out)
        map_file = maps_out_dir+nick+'.pdf'
        # consolidate annotated genbank file
        record = load_fasta(fasta_out)
        aa_defs = load_multifasta(annot_aa)
        features = []
        counter = 1
        ctg_flag_1 = 0
        ctg_flag_2 = 0
        for protein in aa_defs:
            gene_count +=1
            # get feature details from description line
            # necessary because the prodigal output is not parser-friendly
            pattern = re.compile(r'\d+_\d+_\d+_\d+_\d+\s+\S+\s+(\d+)\s+\S+\s+(\d+)\s+\S+\s+(\S*\d)')
            match = pattern.match(protein.description)
            start_pos = int(match.group(1))
            end_pos = int(match.group(2))
            strand_pos = int(match.group(3))
            feat_loc = FeatureLocation(start_pos, end_pos)
            annotation = rec_cogs['Query_'+str(counter)]
            if ctg_flag_1 is 0:
                all_hit_list_handle.write("</ul><br><a href='"
                                          +"../../../../"
                                          +map_file
                                          +"'>Contig "
                                          +nick+"</a><ul>")
                ctg_flag_1 = 1
            all_hit_list_handle.write("<li>"+str(counter)
                                            +'. '+annotation+"</li>")
            # detect phage content in annotation
            phi_pattern = re.compile(r".+(COG\d+).+"
                                      "(phage|capsid|muramidase|tail|"
                                      "replication|helicase|polymerase|"
                                      "integrase|recombinase"
                                      "suppressor|hydrolase|transposase).+",
                                     re.IGNORECASE)
            phi_match = phi_pattern.match(annotation)
            if phi_match:
                hit_flag = 'on'
                hit_dict = {'CDS': counter,
                            'annot': annotation,
                            'COGs': phi_match.group}
                contig_hits[nick].append(hit_dict)
                # write out to summary file
                if ctg_flag_2 is 0:
                    sp_hit_list_handle.write("</ul><br><a href='"
                                             +"../../../../"
                                             +map_file
                                             +"'>Contig "
                                             +nick+"</a><ul>")
                    ctg_flag_2 = 1
                sp_hit_list_handle.write("<li>"+str(counter)
                                          +'. '+annotation+"</li>")
            else:
                hit_flag = 'off'
            # consolidation feature annotations
            quals = {'note': protein.description,
                     'fct': annotation,
                     'flag': hit_flag}
            feature = SeqFeature(location=feat_loc,
                                 strand=strand_pos,
                                 id=protein.id,
                                 type='CDS',
                                 qualifiers=quals)
            features.append(feature)
            counter +=1
        record.features = features
        record.description = dataset['f_nick']+'_contig_'+nick
        record.name = nick
        record.dbxrefs = ['Project:np1']
        record.seq.alphabet = generic_dna
        gbk_out = solid_out_dir+nick+'.gbk'
        write_genbank(gbk_out, record)
        # generate graphical map
        ContigDraw(nick, gbk_out, map_file)
    sp_hit_list_handle.write("</ul>")
    all_hit_list_handle.write("</ul>")
    sp_hit_list_handle.close()
    all_hit_list_handle.close()
    print "\t", gene_count, "predicted genes in", ctg_count, "contigs"

