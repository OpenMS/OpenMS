from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from PILISModel cimport *
from PeptideHit cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISCrossValidation.h>" namespace "OpenMS":
    
    cdef cppclass PILISCrossValidation(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISCrossValidation() nogil except +
        PILISCrossValidation(PILISCrossValidation) nogil except +
        # void setOptions(Map[ String, Option ] & rhs) nogil except +
        # void setOption(String & name, Option & option) nogil except +
        # TODO apply does not compile
        # void apply(Param & PILIS_param, PILISModel & base_model, libcpp_vector[ PILIS_Peptide ] & peptides) nogil except +
        # NESTED STL # DoubleReal scoreHits(libcpp_vector[ libcpp_vector[ libcpp_vector[ MSSpectrum[RichPeak1D] ] ] ] & sim_spectra, libcpp_vector[ libcpp_vector[ MSSpectrum[RichPeak1D] ] ] & exp_spectra) nogil except +


cdef extern from "<OpenMS/ANALYSIS/ID/PILISCrossValidation.h>" namespace "OpenMS":
    
    cdef cppclass PILIS_Peptide "OpenMS::PILISCrossValidation::Peptide":

      # TODO does not compile as attributes -> getter/setters necessary?
      # AASequence sequence
      Int charge
      # MSSpectrum[RichPeak1D] spec
      # libcpp_vector[PeptideHit] hits

cdef extern from "<OpenMS/ANALYSIS/ID/PILISCrossValidation.h>" namespace "OpenMS::PILISCrossValidation":
    
    cdef cppclass PILIS_Option "OpenMS::PILISCrossValidation::Option":
        PILIS_Option() nogil except +
        PILIS_Option(PILIS_Option) nogil except +
        # PILIS_Option_Type type
        Int int_min
        Int int_max
        Int int_stepsize
        DoubleReal dbl_min
        DoubleReal dbl_max
        DoubleReal dbl_stepsize
        PILIS_Option(PILIS_Option_Type t, DoubleReal min_, DoubleReal max_, DoubleReal stepsize) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/PILISCrossValidation.h>" namespace "OpenMS::PILISCrossValidation::Option":
    cdef enum PILIS_Option_Type "OpenMS::PILISCrossValidation::Option::Type":
        #wrap-attach:
        #    PILIS_Option
        INT
        DOUBLE
        BOOL
        STRINGLIST

