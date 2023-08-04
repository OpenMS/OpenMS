from Types cimport *
from ConsensusMap cimport *
from Peak2D cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>" namespace "OpenMS":
    
    cdef cppclass ProteinInference "OpenMS::ProteinInference":
        # wrap-doc:
            #  [experimental class] given a peptide quantitation, infer corresponding protein quantities
            #  
            #  Infers protein ratios from peptide ratios (currently using unique peptides only).
            #  Use the IDMapper class to add protein and peptide information to a
            #  quantitative ConsensusMap prior to this step

        ProteinInference() except + nogil 
        ProteinInference(ProteinInference &) except + nogil 

        void infer(ConsensusMap & consensus_map, UInt reference_map) except + nogil 
            # wrap-doc:
            #  Given a peptide quantitation, infer corresponding protein quantities
            #  
            #  Infers protein ratios from peptide ratios (currently using unique peptides only).
            #  Use the IDMapper class to add protein and peptide information to a
            #  quantitative ConsensusMap prior to this step
            #  
            #  
            #  :param consensus_map: Peptide quantitation with ProteinIdentifications attached, where protein quantitation will be attached
            #  :param reference_map: Index of (iTRAQ) reference channel within the consensus map
