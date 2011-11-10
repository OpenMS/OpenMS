# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow $
# --------------------------------------------------------------------------

import re
import random
import math
import sys
import argparse

## holds FASTA header + sequence
class FASTAEntry:
    pass
		

## grab entries from FASTA file
def nextEntry(fileobj):
				entry = FASTAEntry()
				entry.header = fileobj.readline()	
				entry.sequence = ""
				for line in fileobj:	
								if '>' == line[0]:	
												yield entry	
												entry.header = line	
												entry.sequence = ""
								else:	
												entry.sequence += line	
				yield entry	


## sample abundance from Gaussian in log space
def sampleAbundance(mu=3, sigma=1):
		return math.exp(random.gauss(mu, sigma))


def main(argv):
		## we use ArgumentParser, which requires 2.7
		if sys.version_info < (2, 7):
			raise "This script requires python 2.7 or greater"

		## add weight filtering functionality if BioPython is available
		try:
			from Bio.SeqUtils.ProtParam import ProteinAnalysis	
			has_biopython = 1
		except :
			has_biopython = 0
			
  		
		parser = argparse.ArgumentParser(description='Add abundance to FASTA files.')
		parser.add_argument('infile', type=argparse.FileType('r'), help='Input FASTA file')
		parser.add_argument('outfile', type=argparse.FileType('w'), help='Output FASTA file')
		
		parser.add_argument('--mu', dest='mu', action='store', default=3, help='mean of gaussian in log space')
		parser.add_argument('--sigma', dest='sigma', action='store', default=1, help='sd of gaussian in log space')
		parser.add_argument('--sample', dest='sample', action='store', default=0, help='Number of entries to keep (for sampling a bigger FASTA file)')
		parser.add_argument('--random', dest='random', action='store_true', help='Randomly shuffle entries before sampling (only if --sample is given). If not given, the first \'X\' samples are used.')
		if (has_biopython):
			parser.add_argument('--weight_low', dest='weight_low', action='store', default=0, help='minimum molecular weight of protein')
			parser.add_argument('--weight_up', dest='weight_up', action='store', default=0, help='Maximum molecular weight of protein (use 0 for unlimited)')
		else:
			print "Warning: protein weight filtering not supported, as BioPython module is not installed."
			
		## argument parsing
		args = parser.parse_args()
		fileobj = args.infile
		fileoutobj = args.outfile
		sample_size = int(args.sample)
		sample_random = bool(args.random)
		if (has_biopython):
			weight_low = float(args.weight_low)
			weight_up = float(args.weight_up)
			if (weight_up <= 0): weight_up = sys.float_info.max
			
		
		## list of final entries
		fasta_entries = []
		
		for entry in nextEntry(fileobj):
				header = entry.header
				## check if it contains 'intensity'?
				rep = re.compile(r"\[# *(.*) *#\]")
				m = rep.search(header)
				header_new = ""
				other = []
				if (m):
					header_new = header.replace(m.group(0), "") ## delete meta
					for element in m.group(1).split(','):
							#print "element:", element
							if (element.find("intensity") == -1):
									other.append(element)
				else:
					header_new = header	## nothing to replace

				## create new metainfo array
				i = "intensity=" + str(sampleAbundance(float(args.mu), float(args.sigma)))
				other.append(i)

				entry.header = header_new.rstrip() + "[# " + (", ").join(other) + " #]"
				
				if (has_biopython):
					sequence = "".join(entry.sequence.split("\n"))
					##
					## BioPython does not like some AA letters - they need replacement
					##
					## replace "U" (Selenocystein) with "C" (Cystein)
					sequence = sequence.replace("U","C")
					## replace "X" (unknown) with "P" (Proline) [arbitrary choice - but weight of 115 is very close to averagine]
					sequence = sequence.replace("X","P")
					## replace "B" (Asparagine or aspartic acid) with "N" (Asparagine) 
					sequence = sequence.replace("B","N")
					## replace "Z" (Glutamine or glutamic acid) with "Q" (Glutamine) 
					sequence = sequence.replace("Z","Q")
					## replace "Z" (Glutamine or glutamic acid) with "Q" (Glutamine) 
					sequence = sequence.replace("Z","Q")
					## replace "J" (Leucine or Isoleucine) with "L" (Leucine) 
					sequence = sequence.replace("J","L")
					analysed_seq = ProteinAnalysis(sequence)
					weight = analysed_seq.molecular_weight()
					if (not(weight_low <= weight and weight <= weight_up)):
						continue
				
				
				fasta_entries.append(entry.header + "\n" + entry.sequence)
				 
				## only read to sample size (the rest is thrown away anyways)
				if (sample_size > 0 and not(sample_random)):
					if (len(fasta_entries) >= sample_size):
						break
					
				
		## select subset (if required)		
		if (sample_size > 0):
			indices = range(0,len(fasta_entries))
			## random sampling only makes sense if we take a subset
			if (sample_random and sample_size < len(fasta_entries)):
				random.shuffle(indices)
			indices = [indices[i] for i in range(0,sample_size)]
			fasta_entries = [fasta_entries[i] for i in indices]
			
		## write to file
		for entry in fasta_entries:
			fileoutobj.write(entry)		

if __name__ == "__main__":
		main(sys.argv)
