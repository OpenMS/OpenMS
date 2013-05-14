from smart_ptr cimport shared_ptr
from TargetedExperiment cimport *
from OpenSwathDataStructures cimport *
from LightTargetedExperiment cimport LightTargetedExperiment, LightPeptide

from MSSpectrum cimport *
from Peak1D cimport *
from MSChromatogram cimport *
from ChromatogramPeak cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>" namespace "OpenMS":

    cdef cppclass OpenSwathDataAccessHelper:

        # OpenSwathDataAccessHelper() nogil except +
        # OpenSwathDataAccessHelper(OpenSwathDataAccessHelper) nogil except + # wrap-ignore

        # Convert a SpectrumPtr to an OpenMS Spectrum
        void convertToOpenMSSpectrum(shared_ptr[Spectrum] sptr, MSSpectrum[Peak1D] & spectrum) nogil except +

        # Convert an OpenMS Spectrum to an SpectrumPtr
        shared_ptr[Spectrum] convertToSpectrumPtr(MSSpectrum[Peak1D] spectrum) nogil except + # wrap-ignore

        # Convert a ChromatogramPtr to an OpenMS Chromatogram
        void convertToOpenMSChromatogram(MSChromatogram[ChromatogramPeak] & chromatogram, shared_ptr[Chromatogram] cptr) nogil except +

        # convert from the OpenMS Targeted experiment to the light Targeted Experiment
        void convertTargetedExp(TargetedExperiment & transition_exp_, LightTargetedExperiment & transition_exp)

        # convert from the LightPeptide to an OpenMS AASequence (with correct modifications)
        void convertPeptideToAASequence(LightPeptide & peptide, AASequence & aa_sequence) nogil except +

