#!/usr/bin/env python

from argparse import ArgumentParser
from bs4 import BeautifulSoup
import csv

# argparser    
parser = ArgumentParser(description = 'Parses MEME-MAST xml output.',
    usage = 'python mast2tsv.py mast.output.xml mast.output.tsv')
# positional arguments
parser.add_argument('infile',
    help = 'XML file from MAST.')
parser.add_argument('outfile',
    help = 'Output filename.')
# parse
args = parser.parse_args()

try:
	with open(args.infile, 'r') as f:
		xml = BeautifulSoup(f)

		# parse motifs
		m = xml.mast.motifs.findAll("motif")
		motifs = {}
		for motif in m:
			motifs[str(int(motif['alt'][5:6]) - 1)] = str(motif['id'])

		# parse sequences
		seqs = xml.mast.sequences.findAll("sequence")
		output = [['peakName', 'motifPvalue', 'motifPosition','motif']]
		for seq in seqs:
			if seq.hit:
				seq_hits = seq.find_all('hit')
				for hit in seq_hits:
					output.append([str(seq['name']), float(hit['pvalue']), str(hit['pos']), motifs[hit['idx']]])
			else:
				output.append([str(seq['name']), float(seq.score['combined_pvalue']), str('NA'), str('NA')])
except IOError:
    print(" '%s' file not openable." % args.infile)
    sys.exit(0)

try:
    with open(args.outfile, 'w') as f:
            wr = csv.writer(f, delimiter = '\t', lineterminator='\n')
            for line in output:
                wr.writerow(line)
except IOError:
    print(" '%s' file not writable." % args.outfile)
    sys.exit(0)
