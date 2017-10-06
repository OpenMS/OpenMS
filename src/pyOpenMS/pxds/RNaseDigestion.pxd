from Types cimport *
from String cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/RNaseDigestion.h>" namespace "OpenMS":

    cdef cppclass RNaseDigestion:

      RNaseDigestion() nogil except + #wrap-ignore
      RNaseDigestion(RNaseDigestion) nogil except + #wrap-ignore

      SignedSize getMissedCleavages() nogil except +
      void setMissedCleavages(SignedSize missed_cleavages) nogil except +

      String getEnzymeName() nogil except +
      void setEnzyme(String name) nogil except +

      void digest(String& rna, libcpp_vector[String]& output, Size min_length, Size max_length) nogil except +

      # Returns the specificity for the digestion
      Specificity getSpecificity() nogil except +
      # Sets the specificity for the digestion (default is SPEC_FULL).
      void setSpecificity(Specificity spec) nogil except +
      # convert spec string name to enum
      Specificity getSpecificityByName(String name) nogil except +

      bool isValidProduct(String protein, Size pep_pos, Size pep_length,
                          bool ignore_missed_cleavages) nogil except +

      # void digestUnmodified(StringView sequence, libcpp_vector[ StringView ] & output, Size min_length, Size max_length) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/RNaseDigestion.h>" namespace "OpenMS::RNaseDigestion":

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
        #    RNaseDigestion
        SPEC_FULL,    # fully enzyme specific, e.g., tryptic (ends with KR, AA-before is KR), or peptide is at protein terminal ends
        SPEC_SEMI,    # semi specific, i.e., one of the two cleavage sites must fulfill requirements
        SPEC_NONE,    # no requirements on start / end
        SIZE_OF_SPECIFICITY
