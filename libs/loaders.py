from Bio import SeqIO, GenBank
import numpy

def load_agnostic(seqfile):
    """Load single-record file of unspecified format."""
    # present version is single-record fasta or genbank only
    while True:
        try: seq_record = load_fasta(seqfile)
        except: print "\tnot single-record fasta"
        else:
            rec_type = 'fasta'
            print "\tfound a fasta file"
            return seq_record, rec_type
            break
        try: seq_record = load_genbank(seqfile)
        except: print "not genbank"
        else:
            rec_type = 'genbank'
            print "\tfound a genbank file"
            return seq_record, rec_type
            break
        raise Exception("Cannot open query file!")

def load_multifasta(seqfile):
    """Load multiple records from Fasta file."""
    from Bio import SeqIO
    input_handle = open(seqfile, 'rU')
    multifasta = SeqIO.parse(input_handle, 'fasta')
    fasta_list = list(multifasta)
    for record in fasta_list:
        assert record.id
    return fasta_list

def load_fasta(seqfile):
    """Load single-record Fasta file."""
    input_handle = open(seqfile, 'rU')
    fasta_record = SeqIO.read(input_handle, 'fasta')
    assert fasta_record.id
    input_handle.close()
    return fasta_record

def load_genbank(seqfile):
    """Load single-record GenBank file."""
    parser = GenBank.FeatureParser()
    input_handle = open(seqfile, 'rU')
    gb_record = parser.parse(input_handle)
    input_handle.close()
    return gb_record

def batch_iterator(iterator, batch_size) :
    # copied from http://biopython.org/wiki/Split_large_file
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def read_array(filename, dtype, separator='\t'):
    # From Numpy cookbook
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    cast = numpy.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return numpy.rec.array(data, dtype=dtype)