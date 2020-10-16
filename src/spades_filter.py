#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Filter a sapdes assembly result on cover and lenght"""

import sys
import argparse
from Bio import SeqIO
import os

desc = "Filter a spades assembly results"
command = argparse.ArgumentParser(prog='spades_filter.py', \
    description=desc, usage='%(prog)s [options] contigs')
command.add_argument('-c','--cover', type=int, \
    nargs='?', default=5, \
    help='Minimun coverage of the contig, default:5')
command.add_argument('-l','--lenght', type=int, \
    nargs='?', default=300, \
    help='Minimun lenght of the contig, default:300')
command.add_argument('-o','--output', type=argparse.FileType("w"), \
    nargs='?', default=sys.stdout, \
    help='Output filter assembly fasta, default:stdout')
command.add_argument('contigs', type=argparse.FileType("r"), \
    help='Input assembly fasta file from Spades')

if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    if args.lenght < 0 :
        raise Exception("Min lenght must be superior to 0")
    output = args.output
    
    for seq in SeqIO.parse(args.contigs, "fasta"):
        h = seq.id.split("_")
        if len(h) != 6 or h[0] != "NODE" or h[2] != "length" or h[4] != "cov":
            raise Exception("This seems not a assmbly fasta file from Spades\n"+\
                            seq.id + " is not in format : NODE_1_length_324709_cov_39.9142\n")
        if int(h[3]) >= args.lenght and float(h[5]) >= args.cover:
            SeqIO.write(seq, output, "fasta")
            
