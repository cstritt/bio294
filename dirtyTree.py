#!/usr/bin/env python3
# -*- coding: utf-8 -*-



""" create fast & dirty gene tree

Created on Fri Nov  8 10:17:51 2019
@author: cristobal
"""


import sys
import strumenti

from Bio import AlignIO

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo


def fastagrep(gene, path_to_fasta, outname):
    """ Fetch fasta entry/ies from a multifasta file.
    """
    
    switch = 0
    gene_seq = ''
    
    with open(path_to_fasta) as f:
        for line in f:
            if line.startswith('>'):
                
                if isinstance(gene, str):
                    switch = 1 if gene in line else 0
                    
                elif isinstance(gene, list):
                    switch = 1 if any([x in line for x in gene]) else 0
                    
            if switch == 1:
                gene_seq += line
            
    with open(outname, 'w') as g:
        g.write(gene_seq)
        
        
        
def local_blast(gene_fasta, subject):
    """ Creates a file in the working directory (.blastout) and additionally
    returns a list with the blast results.
    """
    
    blast_cline = NcbiblastpCommandline(query = gene_fasta, 
                                         subject = subject,
                                         outfmt = '6 pident length evalue sseqid', 
                                         out = gene_fasta + ".blastout")
    blast_cline()
    
    blast_res = []
    with open(gene_fasta + '.blastout') as f:
        for line in f:
            blast_res.append(line.strip().split())
            
    return blast_res


def filter_blast_results(path_to_blast_results, 
                         min_perc_id = 30, 
                         min_aln_len = 150, 
                         max_eval = 0.0001):
    
    """ Returns names of genes which have passed the filtering
    """
    homologs = set()
    
    with open(path_to_blast_results) as f:
        for line in f:
            res = line.strip().split()

            pident, l, evalu, seqid = res
            if int(l) > min_aln_len and float(pident) > min_perc_id and float(evalu) < max_eval:
                homologs.add(seqid)
    
    return list(homologs)
    


def muscle_align(fasta_in, outname):
    """
    """
    
    cline = MuscleCommandline(input=fasta_in, 
                          out=outname)

    cline()
    aln = AlignIO.read(outname, 'fasta')
    return aln
        


def main():
    
    try:
        gene = sys.argv[1]
        gene_space = sys.argv[2]
        
    except IndexError:
        sys.exit("Usage: dirtyTree.py <gene name> <gene space>")
        
    
    # Fetch fasta entry of the focal gene
    fastagrep(gene, gene_space, gene + ".fasta")
    
    
    # Blast and filter results
    blast_res = local_blast(gene + ".fasta", gene_space)

    homologs = filter_blast_results(gene + ".blastout")
    if len(homologs)  == 0:
        sys.exit("No homologs found for %s, exiting ..." % (gene))
    
    
    # Fetch sequences of homologs, align them, and trim the alignment    
    fastagrep(homologs, gene_space, gene + ".homologs.fasta")
    
    aln = muscle_align(gene + ".homologs.fasta", gene + ".homologs.aligned.fasta")
    
    aln_trimmed = strumenti.trim_alignment(aln, max_prop_missing = 0.2)
    
    
    # Estimate a gene tree using the neighbor-joining algorithm        
    tree = strumenti.neighbor_joining_tree(aln_trimmed, 'blosum62')
    
    sys.stdout.write("\nHomolog tree for %s\n" % (gene))
    Phylo.draw_ascii(tree, file=sys.stdout, column_width=100)
    sys.stderr.write("\n%i positions were used to estimate the tree.\n" % (aln_trimmed.get_alignment_length()))
    
if __name__ == '__main__':
    main()
    