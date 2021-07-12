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
        #
        #     Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
        #     Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
        #     due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
        #     peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
        #     Thus no random selection of just n specific missed cleavage sites is performed.

      ProteaseDigestion() nogil except +
      ProteaseDigestion(ProteaseDigestion) nogil except +

      void setEnzyme(String name) nogil except + # wrap-doc:Sets the enzyme for the digestion (by name)

      Size digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +

      Size digest(AASequence & protein, libcpp_vector[AASequence] & output, Size min_length, Size max_length) nogil except +
          # wrap-doc:
          #     Performs the enzymatic digestion of a protein.
          #
          #     :param protein: Sequence to digest
          #     :param output: Digestion products (peptides)
          #     :param min_length: Minimal length of reported products
          #     :param max_length: Maximal length of reported products (0 = no restriction)
          #     :return: Number of discarded digestion products (which are not matching length restrictions)

      Size peptideCount(AASequence & protein) nogil except + # wrap-doc:Returns the number of peptides a digestion of protein would yield under the current enzyme and missed cleavage settings

      bool isValidProduct(AASequence protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except +
          # wrap-doc:
          #     Variant of EnzymaticDigestion::isValidProduct() with support for n-term protein cleavage and random D|P cleavage
          #
          #     Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the flags provided here
          #
          #     :param protein: Protein sequence
          #     :param pep_pos: Starting index of potential peptide
          #     :param pep_length: Length of potential peptide
          #     :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
          #     :param allow_nterm_protein_cleavage: Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 and protein starts with 'M'
          #     :param allow_random_asp_pro_cleavage: Allow cleavage at D|P sites to count as n/c-terminal
          #     :return: True if peptide has correct n/c terminals (according to enzyme, specificity and above flags)

      bool isValidProduct(String protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except + # wrap-doc:Forwards to isValidProduct using protein.toUnmodifiedString()


