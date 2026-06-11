#!/usr/bin/env python
import sys
from pyopenms import *

def main(options):
    
    # read fasta file  
    ff = FASTAFile()
    ff.readStart(options.infile)
    fe = FASTAEntry()

    # use ProteaseDigestion class 
    digestor = ProteaseDigestion()
    digestor.setEnzyme(options.enzyme)
    digestor.setMissedCleavages(options.missed_cleavages)

    # protein and peptide datastructure
    protein_identifications = []
    peptide_identifications = []
    protein_identification = ProteinIdentification()
    protein_identifications.append(protein_identification)
    temp_pe = PeptideEvidence()

    # number of dropped peptides due to length restriction
    dropped_by_length = 0

    while(ff.readNext(fe)):  
        # construct ProteinHit and fill it with sequence information
        temp_protein_hit = ProteinHit()
        temp_protein_hit.setSequence(fe.sequence)
        temp_protein_hit.setAccession(fe.identifier)

        # save the ProteinHit in a ProteinIdentification Object 
        protein_identification.insertHit(temp_protein_hit)
       
        # construct a PeptideHit and save the ProteinEvidence (Mapping) for the specific 
        # current protein
        temp_peptide_hit = PeptideHit()
        temp_pe.setProteinAccession(fe.identifier);
        temp_peptide_hit.setPeptideEvidences([temp_pe])

        # digestion 
        current_digest = []
        aaseq = AASequence()
        if (options.enzyme == "none"):
            current_digest.append(aaseq.fromString(fe.sequence))
        else:
            dropped_by_length += digestor.digest(aaseq.fromString(fe.sequence), current_digest, options.min_length, options.max_length)

        for seq in current_digest:
            # fill the PeptideHit and PeptideIdentification datastructure
            peptide_identification = PeptideIdentification() 
            temp_peptide_hit.setSequence(seq)
            peptide_identification.insertHit(temp_peptide_hit)
            peptide_identifications.append(peptide_identification)

    print(str(dropped_by_length) + " peptides have been dropped due to the length restriction.")
    idxml = IdXMLFile()
    idxml.store(options.outfile, protein_identifications, peptide_identifications)
        
def handle_args():
    import argparse 
    
    # get available enzymes from ProteaseDB (in OpenMS)
    all_enzymes = []
    p_db=ProteaseDB().getAllNames(all_enzymes)

    usage = ""
    usage += "\nProteinDigestor −− In silico digestion of proteins."

    parser = argparse.ArgumentParser(description = usage)
    parser.add_argument('-in', dest='infile', help='An input file containing amino acid sequences [fasta]')
    parser.add_argument('-out', dest='outfile', help='Output digested sequences in idXML format [idXML]')
    parser.add_argument('-enzyme', dest='enzyme', help='Enzymes which can be used for digestion: '+ ', '.join(map(bytes.decode, all_enzymes)))
    parser.add_argument('-min_length', type=int, dest='min_length', help ='Minimum length of peptide')
    parser.add_argument('-max_length', type=int, dest='max_length', help='Maximum length of peptide')
    parser.add_argument('-missed_cleavages', type=int, dest='missed_cleavages', help='The number of allowed missed cleavages')

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)

