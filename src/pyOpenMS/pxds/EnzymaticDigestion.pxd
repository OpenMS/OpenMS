from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS":

    cdef cppclass EnzymaticDigestion:

      EnzymaticDigestion() nogil except +
      EnzymaticDigestion(EnzymaticDigestion) nogil except + #wrap-ignore

      SignedSize getMissedCleavages() nogil except +
      void setMissedCleavages(SignedSize missed_cleavages) nogil except +

      # not wrapped due to name clash with Enzyme.h

      # Enzyme getEnzyme()  nogil except +
      # void setEnzyme(Enzyme enzyme) nogil except +
      # Enzyme getEnzymeByName(String & name) nogil except +

      String getEnzymeName() nogil except +
      void setEnzyme(String name) nogil except +

      void digest(AASequence & protein, libcpp_vector[AASequence] & output) nogil except +
      Size peptideCount(AASequence & protein) nogil except +

      # Returns the specificity for the digestion
      Specificity getSpecificity() nogil except +
      # Sets the specificity for the digestion (default is SPEC_FULL).
      void setSpecificity(Specificity spec) nogil except +
      # convert spec string name to enum
      Specificity getSpecificityByName(String name)nogil except +

      bool isValidProduct(AASequence protein, Size pep_pos, Size pep_length, 
                          bool methionine_cleavage, bool ignore_missed_cleavages) nogil except +
      bool isValidProduct(String protein, Size pep_pos, Size pep_length, 
                          bool methionine_cleavage, bool ignore_missed_cleavages) nogil except +
						  
      # void digestUnmodifiedString(StringView sequence, libcpp_vector[ StringView ] & output, Size min_length, Size max_length) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/EnzymaticDigestion.h>" namespace "OpenMS::EnzymaticDigestion":

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
        #    EnzymaticDigestion
      SPEC_FULL,    # fully enzyme specific, e.g., tryptic (ends with KR, AA-before is KR), or peptide is at protein terminal ends
      SPEC_SEMI,    # semi specific, i.e., one of the two cleavage sites must fulfill requirements
      SPEC_NONE,    # no requirements on start / end
      SIZE_OF_SPECIFICITY

