#!/usr/bin/env python
"""Cut up fasta file in non-overlapping or overlapping parts of equal length.
"""
import argparse
from Bio import SeqIO


def cut_up_fasta(fastfiles, chunk_size, overlap, merge_last):
    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            if (not merge_last and len(record.seq) > chunk_size) or (merge_last and len(record.seq) >= 2 * chunk_size):
                i = 0
                for split_seq in chunks(record.seq, chunk_size, overlap, merge_last):
                    print ">%s.%i\n%s" % (record.id, i, split_seq)
                    i = i + 1
            else:
                print ">%s\n%s" % (record.id, record.seq)


def chunks(l, n, o, merge_last):
    """ Yield successive n-sized chunks from l with given overlap o between the
    chunks.
    """
    assert n > o

    if not merge_last:
        for i in xrange(0, len(l), n - o):
            yield l[i:i + n]
    else:
        for i in xrange(0, len(l) - n + 1, n - o):
            yield l[i:i + n] if i + n + n - o <= len(l) else l[i:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-c", "--chunk_size", default=1999, type=int, help="Chunk size\n")
    parser.add_argument("-o", "--overlap_size", default=1900, type=int, help="Overlap size\n")
    parser.add_argument("-m", "--merge_last", default=False, action="store_true", help="Concatenate final part to last contig\n")
    args = parser.parse_args()
    cut_up_fasta(args.contigs, args.chunk_size, args.overlap_size, args.merge_last)
