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

# this class has addons, see the ./addons folder (../addons/OpenSwathDataAccessHelper.pyx)

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>" namespace "OpenMS":

    cdef cppclass OpenSwathDataAccessHelper:

        OpenSwathDataAccessHelper() nogil except + # compiler
        OpenSwathDataAccessHelper(OpenSwathDataAccessHelper &) nogil except + # compiler

        void convertToOpenMSSpectrum(shared_ptr[OSSpectrum] sptr, MSSpectrum & spectrum) nogil except + # wrap-doc:Converts a SpectrumPtr to an OpenMS Spectrum

        shared_ptr[OSSpectrum] convertToSpectrumPtr(MSSpectrum spectrum) nogil except + # wrap-ignore

        shared_ptr[OSChromatogram] convertToChromatogramPtr(MSChromatogram chromatogram) nogil except + # wrap-ignore

        void convertToOpenMSChromatogram(shared_ptr[OSChromatogram] cptr, MSChromatogram & chromatogram) nogil except + # wrap-doc:Converts a ChromatogramPtr to an OpenMS Chromatogram
        void convertToOpenMSChromatogramFilter(MSChromatogram & chromatogram, shared_ptr[OSChromatogram] cptr,
                                               double rt_min, double rt_max) nogil except +

        void convertTargetedExp(TargetedExperiment & transition_exp_, LightTargetedExperiment & transition_exp) nogil except + # wrap-doc:Converts from the OpenMS TargetedExperiment to the OpenMs LightTargetedExperiment

        void convertPeptideToAASequence(LightCompound & peptide, AASequence & aa_sequence) nogil except + # wrap-doc:Converts from the LightCompound to an OpenMS AASequence (with correct modifications)
        
        void convertTargetedCompound(Peptide pep, LightCompound & p) nogil except + # wrap-doc:Converts from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide

