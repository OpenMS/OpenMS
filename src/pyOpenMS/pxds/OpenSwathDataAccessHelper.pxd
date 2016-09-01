from smart_ptr cimport shared_ptr
from TargetedExperiment cimport *
from TargetedExperimentHelper cimport *
from OpenSwathDataStructures cimport *
from LightTargetedExperiment cimport LightTargetedExperiment, LightCompound

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
        void convertToOpenMSSpectrum(shared_ptr[OSSpectrum] sptr, MSSpectrum[Peak1D] & spectrum) nogil except +

        # Convert an OpenMS Spectrum to an SpectrumPtr
        shared_ptr[OSSpectrum] convertToSpectrumPtr(MSSpectrum[Peak1D] spectrum) nogil except + # wrap-ignore

        # Convert a ChromatogramPtr to an OpenMS Chromatogram
        void convertToOpenMSChromatogram(MSChromatogram[ChromatogramPeak] & chromatogram, shared_ptr[OSChromatogram] cptr) nogil except +

        # convert from the OpenMS Targeted experiment to the light Targeted Experiment
        void convertTargetedExp(TargetedExperiment & transition_exp_, LightTargetedExperiment & transition_exp) nogil except +

        # convert from the LightCompound to an OpenMS AASequence (with correct modifications)
        void convertPeptideToAASequence(LightCompound & peptide, AASequence & aa_sequence) nogil except +
        
        # convert from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide
        void convertTargetedCompound(Peptide pep, LightCompound & p) nogil except +

