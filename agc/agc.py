#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "STERLIN Aude"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["STERLIN Aude"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "STERLIN Aude"
__email__ = "sterlin.aude@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()
#==============================================================
# 1 : Dé-duplication en séquence "complète"

def read_fasta(amplicon_file, minseqlen):
    '''
    Sequences generator.
    '''
    with gzip.open(amplicon_file, 'rb') as file:
        lines=[line.decode().strip() for line in file]

    names_index=[]
    sequences=[]
    
    for i in range(len(lines)):
        line=lines[i]

        if line[0] == '>':
            names_index.append(i)

    for i in range (len(names_index)):

        start_index=names_index[i]

        if i!=len(names_index)-1:

            end_index=names_index[i+1]
            sequence=''.join([lines[j] for j in range(start_index,end_index)])

        else :

            end_index=len(lines)
            sequence=''.join([lines[j] for j in range(start_index,end_index)])

        sequences.append(sequence)


    end_of_name_index=[sequences[i].find('fastq') for i in range(len(sequences))]
    names=[sequences[i][:end_of_name_index[i]+5] for i in range(len(sequences))]

    clean_sequences=[sequences[i][end_of_name_index[i]+5:] for i in range(len(sequences))]

    for sequence in clean_sequences:
        if len(sequence)>=minseqlen:
            yield sequence

def dereplication_fulllength(amplicon_file, minseqlen, mincount, max_studied=1000):
    unique_sequences=[]
    occurences=[]
#     i=0

    sequences=read_fasta(amplicon_file, minseqlen)
    for sequence in sequences:

        if sequence not in unique_sequences:

            unique_sequences.append(sequence)
            occurences.append(1)
        else : 
            index=unique_sequences.index(sequence)
            occurences[index]=occurences[index]+1
#         i+=1
#         print(len(occurences))
#         if len(occurences)>max_studied:
#             break

    zipped=sorted(zip(occurences,unique_sequences),reverse=True)
    unique_sorted=[seq for _,seq in zipped]
    occurences_sorted=[occ for occ,_ in zipped]
    print(occurences_sorted)
    for i in range(len(occurences_sorted)):
        if occurences_sorted[i]>mincount:
            yield [unique_sequences[i], occurences_sorted[i]]

#==============================================================
# Part 2

def get_chunks(sequence, chunk_size):
    
    N_sous_sequences=len(sequence)//chunk_size
    
    sous_sequences=[sequence[i*chunk_size:(i+1)*chunk_size] 
                    for i in range(N_sous_sequences)]
    
    return sous_sequences
    
# print(get_chunks('AAATTTCCCGGGAAATTTGGGCCC',5))

def cut_kmer(sequence, kmer_size):
    unique_kmers=[]
    for i in range(len(sequence)):
        kmer=sequence[i:i+kmer_size]
        if kmer not in unique_kmers and len(kmer)==kmer_size:
            yield kmer

# resu=cut_kmer('AAATTTGCCC',5)
# for res in resu:
#     print(res)

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    cut_kmer_result=cut_kmer(sequence, kmer_size)
    for kmer in cut_kmer_result:
        if kmer not in list(kmer_dict.keys()):
            kmer_dict[kmer]=[id_seq]
        else: 
            if id_seq not in kmer_dict[kmer]:
                kmer_dict[kmer]=kmer_dict[kmer].append(id_seq)
            
    return kmer_dict

# kmer_dict={}
# sequence='AAATTTTTTTCGCGGGAAAA'
# id_seq=0
# kmer_size=3
# print(get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size))

def search_mates(kmer_dict, sequence, kmer_size):
    
    return [i[0] for i in 
            Counter([ids for kmer in cut_kmer(sequence, kmer_size) 
                     if kmer in kmer_dict 
                     for ids in kmer_dict[kmer]]).most_common(8)]

def get_identity(alignment_list):
    N_identical=0
    seq0=alignment_list[0]
    seq1=alignment_list[1]
    for i in range(len(seq0)):
        if seq0[i]==seq1[i]:
            N_identical+=1
            
    return N_identical/len(seq0)

# print(get_identity(['AAATGCGTAA','AAA___GTAA']))
def get_unique(ids):
    return {}.fromkeys(ids).keys()

def detect_chimera(perc_identity_matrix):
    

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    output_file = args.output_file
    
    # Part 1:
#     sequences=read_fasta(amplicon_file, minseqlen)
#     i=0
#     for sequence in sequences :
#         print(sequence[0:100])
#         i=i+1
# #         if i > 10:
# #             break
#     result=dereplication_fulllength(amplicon_file, minseqlen,mincount)
# #     for res in result:
# #         print(res)


if __name__ == '__main__':
    main()