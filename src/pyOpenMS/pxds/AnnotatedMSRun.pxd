from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from Types cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/METADATA/AnnotatedMSRun.h>" namespace "OpenMS":
    
    cdef cppclass AnnotatedMSRun:
        # wrap-doc:
        #  Class for storing MS run data with peptide and protein identifications
        #
        #  This class stores an MSExperiment (containing spectra) along with peptide and protein
        #  identifications. Each spectrum in the MSExperiment is associated with a single
        #  PeptideIdentification object. Object gets typically not manually created but generated
        #  by the IDMapper class.
        #
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    run = AnnotatedMSRun()
        #    exp = MSExperiment()
        #    MzMLFile().load(path_to_file, exp)
        #    run.setMSExperiment(exp)
        #    run.setPeptideIdentifications(my_peptide_ids)
        
        AnnotatedMSRun() nogil except +
        AnnotatedMSRun(MSExperiment) nogil except +
        AnnotatedMSRun(AnnotatedMSRun) nogil except +
        
        # Protein identification methods
        libcpp_vector[ProteinIdentification]& getProteinIdentifications() nogil except +
        const libcpp_vector[ProteinIdentification]& getProteinIdentifications() nogil except + # wrap-ignore
        
        # Peptide identification methods
        libcpp_vector[PeptideIdentification]& getPeptideIdentifications() nogil except +
        const libcpp_vector[PeptideIdentification]& getPeptideIdentifications() nogil except + # wrap-ignore
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification]& ids) nogil except +
        
        # MSExperiment methods
        MSExperiment& getMSExperiment() nogil except +
        const MSExperiment& getMSExperiment() nogil except + # wrap-ignore
        void setMSExperiment(MSExperiment& experiment) nogil except +
        
        # Access methods
        libcpp_pair[MSSpectrum&, PeptideIdentification&] operator[](size_t idx) nogil except + # wrap-ignore
        
