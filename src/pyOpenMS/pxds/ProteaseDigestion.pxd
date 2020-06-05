from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDigestion.h>" namespace "OpenMS":

    cdef cppclass ProteaseDigestion:

      ProteaseDigestion() nogil except +
      ProteaseDigestion(ProteaseDigestion) nogil except + #wrap-ignore

      SignedSize getMissedCleavages() nogil except +
      void setMissedCleavages(SignedSize missed_cleavages) nogil except +

      # not wrapped due to name clash with DigestionEnzyme.h (???)

      # void setEnzyme(DigestionEnzyme* enzyme) nogil except +

      String getEnzymeName() nogil except +
      void setEnzyme(String name) nogil except +

      Size digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +
      Size digest(AASequence & protein, libcpp_vector[AASequence] & output, Size min_length, Size max_length) nogil except +
      Size peptideCount(AASequence & protein) nogil except +

      # Returns the specificity for the digestion
      Specificity getSpecificity() nogil except +
      # Sets the specificity for the digestion (default is SPEC_FULL).
      void setSpecificity(Specificity spec) nogil except +
      # convert spec string name to enum
      Specificity getSpecificityByName(String name) nogil except +

      bool isValidProduct(AASequence protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except +
      bool isValidProduct(String protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages, bool methionine_cleavage) nogil except +

      # void digestUnmodified(StringView sequence, libcpp_vector[ StringView ] & output, Size min_length, Size max_length) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/ProteaseDigestion.h>" namespace "OpenMS::ProteaseDigestion":

    # protected
    # cdef cppclass BindingSite:
    #   BindingSite()
    #   Size position
    #   String AAname

    # cdef cppclass CleavageModel:
    #   CleavageModel()
    #   double p_cleave
    #   double p_miss

    cdef enum Specificity:
        # wrap-attach:
        #    ProteaseDigestion
        SPEC_NONE = 0, # no requirements on start / end
        SPEC_SEMI = 1, # semi specific, i.e., one of the two cleavage sites must fulfill requirements
        SPEC_FULL = 2, # fully enzyme specific, e.g., tryptic (ends with KR, AA-before is KR), or peptide is at protein terminal ends
        SPEC_UNKNOWN = 3
        SPEC_NOCTERM = 8, # no requirements on CTerm (currently not supported in the class)
        SPEC_NONTERM = 9, # no requirements on NTerm (currently not supported in the class)
        SIZE_OF_SPECIFICITY = 10
