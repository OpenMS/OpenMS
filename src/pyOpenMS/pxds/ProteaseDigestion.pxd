from Types cimport *
from String cimport *
from AASequence cimport *
from EnzymaticDigestion cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDigestion.h>" namespace "OpenMS":

    cdef cppclass ProteaseDigestion(EnzymaticDigestion):
        # wrap-inherits:
        #    EnzymaticDigestion
        #
        # wrap-doc:
        #     Class for the enzymatic digestion of proteins
        #     -----
        #     Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
        #     Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
        #     due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
        #     peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
        #     Thus no random selection of just n specific missed cleavage sites is performed.
        #     -----
        #     Usage:
        #         from pyopenms import *
        #         from urllib.request import urlretrieve
        #         #
        #         urlretrieve ("http://www.uniprot.org/uniprot/P02769.fasta", "bsa.fasta")
        #         #
        #         dig = ProteaseDigestion()
        #         dig.setEnzyme('Lys-C')
        #         bsa_string = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
        #         bsa_oms_string = String(bsa_string) # convert python string to OpenMS::String for further processing
        #         #
        #         minlen = 6
        #         maxlen = 30
        #         #
        #         # Using AASequence and digest
        #         result_digest = []
        #         result_digest_min_max = []
        #         bsa_aaseq = AASequence.fromString(bsa_oms_string)
        #         dig.digest(bsa_aaseq, result_digest)
        #         dig.digest(bsa_aaseq, result_digest_min_max, minlen, maxlen)
        #         print(result_digest[4].toString()) # GLVLIAFSQYLQQCPFDEHVK
        #         print(len(result_digest)) # 57 peptides
        #         print(result_digest_min_max[4].toString()) # LVNELTEFAK
        #         print(len(result_digest_min_max)) # 42 peptides
        #         #
        #         # Using digestUnmodified without the need for AASequence from the EnzymaticDigestion base class
        #         result_digest_unmodified = []
        #         dig.digestUnmodified(StringView(bsa_oms_string), result_digest_unmodified, minlen, maxlen)
        #         print(result_digest_unmodified[4].getString()) # LVNELTEFAK
        #         print(len(result_digest_unmodified)) # 42 peptides

      ProteaseDigestion() nogil except +

      ProteaseDigestion(ProteaseDigestion &) nogil except + 

      void setEnzyme(String name) nogil except + # wrap-doc:Sets the enzyme for the digestion (by name)

      Size digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +

      Size digest(AASequence & protein, libcpp_vector[AASequence] & output, Size min_length, Size max_length) nogil except +
          # wrap-doc:
          #     Performs the enzymatic digestion of a protein.
          #     -----
          #     :param protein: Sequence to digest
          #     :param output: Digestion products (peptides)
          #     :param min_length: Minimal length of reported products
          #     :param max_length: Maximal length of reported products (0 = no restriction)
          #     :returns: Number of discarded digestion products (which are not matching length restrictions)

      Size peptideCount(AASequence & protein) nogil except + # wrap-doc:Returns the number of peptides a digestion of protein would yield under the current enzyme and missed cleavage settings

      bool isValidProduct(AASequence protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except +
          # wrap-doc:
          #     Variant of EnzymaticDigestion::isValidProduct() with support for n-term protein cleavage and random D|P cleavage
          #     -----
          #     Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the flags provided here
          #     -----
          #     :param protein: Protein sequence
          #     :param pep_pos: Starting index of potential peptide
          #     :param pep_length: Length of potential peptide
          #     :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
          #     :param allow_nterm_protein_cleavage: Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 and protein starts with 'M'
          #     :param allow_random_asp_pro_cleavage: Allow cleavage at D|P sites to count as n/c-terminal
          #     :returns: True if peptide has correct n/c terminals (according to enzyme, specificity and above flags)

      bool isValidProduct(String protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except + # wrap-doc:Forwards to isValidProduct using protein.toUnmodifiedString()


