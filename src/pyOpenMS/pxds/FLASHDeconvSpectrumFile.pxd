from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from String cimport *
from MSExperiment cimport *
from DeconvolvedSpectrum cimport *
from FLASHHelperClasses cimport *
from Param cimport *

cdef extern from "<OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>" namespace "OpenMS":

    cdef cppclass FLASHDeconvSpectrumFile:
        # wrap-doc:
        #  FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats.
        #  This class provides static methods for writing deconvolved spectrum data.
        #  Note: Methods taking std::ostream are not directly exposed. Use file-based workflows.

        FLASHDeconvSpectrumFile() except + nogil
        FLASHDeconvSpectrumFile(FLASHDeconvSpectrumFile &) except + nogil

        # Static method for writing mzML files - this takes file paths which works well in Python
        @staticmethod
        void writeMzML(MSExperiment & map_in,
                       libcpp_vector[DeconvolvedSpectrum] & deconvolved_spectra,
                       String & deconvolved_mzML_file,
                       String & annotated_mzML_file,
                       int mzml_charge,
                       libcpp_vector[double] tols) except + nogil
        # wrap-doc:
        #  Write deconvolved and annotated mzML files.
        #  :param map_in: The original MSExperiment
        #  :param deconvolved_spectra: The deconvolved spectra
        #  :param deconvolved_mzML_file: Output path for deconvolved mzML
        #  :param annotated_mzML_file: Output path for annotated mzML
        #  :param mzml_charge: Charge to use in output mzML
        #  :param tols: Mass tolerances per MS level
