#! /usr/bin/python
# -*- coding: utf-8 -*-

"""Compar les metriques de contigs assemblé  
"""

from Bio import SeqIO
import os
import sys
import argparse

desc = "Evaluation of contigs assembly."
command = argparse.ArgumentParser(prog='compar_assembly_contigs.py', \
                                      description=desc, usage='%(prog)s [options] contigs files')
command.add_argument('-n', '--nvalue', nargs="?", \
                         type=int, default=50, \
                         help='Defined N value to computed, default:50')
command.add_argument('-m', '--min', nargs="?", \
                         type=int, default=150, \
                         help='Defined minimun contigs size to be considered, default:150')
command.add_argument('-o', "--output", \
                         type=argparse.FileType("w"), nargs='?', default=sys.stdout, \
                         help='Export comparaison result to : STDOUT')
command.add_argument('files', type=argparse.FileType("r"), nargs='+', \
                         help='Fasta contigs file')

def read_fasta_file(fasta_f, min_len):
    contig = []
    for seq in SeqIO.parse(fasta_f, 'fasta'):
        if len(seq) > min_len:
            contig.append(len(seq))
    contig.sort(reverse=True)
    return contig

def get_N(contig, n):
    """return the len of Nth contigs were n% of total len were attempted
    contigs are order by low to high"""
    seuil = sum(contig)*n/100
    value = 0
    for i in contig:
        value += i
        if value > seuil:
            return i

def print_lines(name, contigs, out, nvalue):
    sep="\t"
    towrite = [name]
    towrite.append(len(contigs))
    towrite.append(sum(contigs))
    towrite.append(max(contigs))
    towrite.append(get_N(contigs, nvalue))
    out.write(sep.join(map(str, towrite)) + "\n")

def print_header(out, nvalue):
    sep="\t"
    out.write(sep.join(["Name", "Nbr Contigs", "Total Lenght", "Best Contigs", "N"+str(nvalue)]) + "\n")


if __name__=='__main__':
    # parse command line
    args = command.parse_args()
    out = args.output
    nvalue = args.nvalue
    if nvalue <1 or nvalue >100:
        raise Exception("N value must be comprise between ")
    print_header(out, nvalue)
    for fasta_f in args.files:
        contigs = read_fasta_file(fasta_f, args.min)
        print_lines(fasta_f.name, contigs, out, nvalue)
    
    # #print calculation
    # print_lines("Files",contigs.keys(), out)
    # print_lines("Nbr contigs",map(len,contigs.values()), out)    
    # print_lines("Total length",map(sum,contigs.values()), out)   
    # print_lines("Best contig",map(max,contigs.values()), out)
    # print_lines("N"+str(nvalue),map(get_N,contigs.values(), [nvalue]*len(contigs)), out) 

    
