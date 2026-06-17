#!/usr/bin/env python
import sys

def main(options):

    # test parameter handling
    print(options.infile, options.outfile, options.enzyme, options.min_length, options.max_length, options.missed_cleavages)

def handle_args():
    import argparse

    usage = ""
    usage += "\nProteinDigestor −− In silico digestion of proteins."

    parser = argparse.ArgumentParser(description = usage)
    parser.add_argument('-in', dest='infile', help='An input file containing amino acid sequences [fasta]')
    parser.add_argument('-out', dest='outfile', help='Output digested sequences in idXML format [idXML]')
    parser.add_argument('-enzyme', dest='enzyme', help='Enzyme used for digestion')
    parser.add_argument('-min_length', type=int, dest='min_length', help ='Minimum length of peptide')
    parser.add_argument('-max_length', type=int, dest='max_length', help='Maximum length of peptide')
    parser.add_argument('-missed_cleavages', type=int, dest='missed_cleavages', help='The number of allowed missed cleavages')

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)

