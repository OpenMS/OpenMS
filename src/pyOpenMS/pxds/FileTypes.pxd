
cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS::FileTypes":

    cdef enum Type:
    
            UNKNOWN,
            DTA,
            DTA2D,
            MZDATA,
            MZXML,
            FEATUREXML,
            IDXML,
            CONSENSUSXML,
            MGF,
            INI,
            TOPPAS,
            TRANSFORMATIONXML,
            MZML,
            MS2,
            PEPXML,
            PROTXML,
            MZIDENTML,
            GELML,
            TRAML,
            MSP,
            OMSSAXML,
            MASCOTXML,
            PNG,
            XMASS,
            TSV,
            PEPLIST,
            HARDKLOER,
            KROENIK,
            FASTA,
            EDTA
            SIZE_OF_TYPE
