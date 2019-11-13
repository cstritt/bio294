#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import Bio
#import sys
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo import TreeConstruction

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def trim_alignment(aln, max_prop_missing = 0):

    """ Given a multisequence alignment, remove all positionswith
    a proportion of missing bases larger than max_prop_missing. 
    """
    
    remove = []
    for position in range(aln.get_alignment_length()):
        letters = [aln[x][position] for x in range(len(aln))]
        prop_missing = letters.count('-') / float(len(letters))
    
        if prop_missing > max_prop_missing:
            remove.append(position)

    #sys.stderr.write("%i out of %i positions will be trimmed\n" % (len(remove), aln.get_alignment_length()))

    aln_trimmed = Bio.Align.MultipleSeqAlignment([])
    
    for entry in aln:
    
        seq_trimmed = [entry.seq[i] for i in range(aln.get_alignment_length()) if i not in remove]
        seq_trimmed = "".join(seq_trimmed)
        seq_rec = SeqRecord(Seq(seq_trimmed), id=entry.id, name=entry.id, description="")
        aln_trimmed.append(seq_rec)
        
    return aln_trimmed
        

def neighbor_joining_tree(aln, prot_model='blosum62'):
    
    """ Estimate a tree from an alignment using neighbor-joining
    """

    calculator = TreeConstruction.DistanceCalculator(prot_model)
    constructor = TreeConstruction.DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    return(tree)

    

def online_blast(fasta_path, 
                 min_aln_len=150, 
                 max_evalue=0.001):
    
    """ fragment: implement filtering...
    """
    
    fasta_string = open(fasta_path).read()
    result_handle = NCBIWWW.qblast("blastp", "refseq_protein", fasta_string)
    
    blast_records = NCBIXML.parse(result_handle)
    blast_records = list(blast_records)
    return blast_records
    
    
def blastp(query, subject):
    
    cmd = ['blastp',
           '-query', query,
           '-subject', subject,
           '-outfmt', '6 pident length evalue sseqid']

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = proc.stdout.read()
    blast_out = output.splitlines()

    return blast_out