from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS":

    cdef cppclass FileTypes:

        # compiler
        FileTypes() except + nogil  # wrap-doc:Centralizes the file types recognized by FileHandler
        FileTypes(FileTypes &) except + nogil  # compiler

        @staticmethod
        String typeToName(FileType t) except + nogil  # wrap-doc:Returns the name/extension of the type

        @staticmethod
        String typeToDescription(FileType t) except + nogil  # wrap-doc:Returns the human-readable explanation of the type

        @staticmethod
        String typeToMZML(FileType t) except + nogil  # wrap-doc:Returns the mzML name

        @staticmethod
        FileType nameToType(String name) except + nogil  
            # wrap-doc:
                #  Converts a file type name into a Type 
                #  
                #  
                #  :param name: A case-insensitive name (e.g. FASTA or Fasta, etc.)

cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS::FileTypes":

    cdef enum class FileType "OpenMS::FileTypes::Type":
        # wrap-attach:
        #    FileTypes
    
          UNKNOWN,            # < Unknown file extension
          DTA,                # < DTA file (.dta)
          DTA2D,              # < DTA2D file (.dta2d)
          MZDATA,             # < MzData file (.mzData)
          MZXML,              # < MzXML file (.mzXML)
          FEATUREXML,         # < %OpenMS feature file (.featureXML)
          IDXML,              # < %OpenMS identification format (.idXML)
          CONSENSUSXML,       # < %OpenMS consensus map format (.consensusXML)
          MGF,                # < Mascot Generic Format (.mgf)
          INI,                # < %OpenMS parameters file (.ini)
          TOPPAS,             # < %OpenMS parameters file with workflow information (.toppas)
          TRANSFORMATIONXML,  # < Transformation description file (.trafoXML)
          MZML,               # < MzML file (.mzML)
          CACHEDMZML,         # < CachedMzML file (.cachedmzML)
          MS2,                # < MS2 file (.ms2)
          PEPXML,             # < TPP pepXML file (.pepXML)
          PROTXML,            # < TPP protXML file (.protXML)
          MZIDENTML,          # < mzIdentML (HUPO PSI AnalysisXML followup format) (.mzid)
          QCML,               # < qcML (will undergo standardisation maybe) (.qcml)
          MZQC,               # < mzQC (HUPO PSI format) (.mzQC)
          GELML,              # < GelML (HUPO PSI format) (.gelML)
          TRAML,              # < TraML (HUPO PSI format) for transitions (.traML)
          MSP,                # < NIST spectra library file format (.msp)
          OMSSAXML,           # < OMSSA XML file format for peptide identifications (.xml)
          MASCOTXML,          # < Mascot XML file format for peptide identifications (.xml)
          PNG,                # < Portable Network Graphics (.png)
          XMASS,              # < XMass Analysis file (fid)
          TSV,                # < msInspect file (.tsv)
          MZTAB,              # < mzTab file (.mzTab)
          PEPLIST,            # < specArray file (.peplist)
          HARDKLOER,          # < hardkloer file (.hardkloer)
          KROENIK,            # < kroenik file (.kroenik)
          FASTA,              # < FASTA file (.fasta)
          EDTA,               # < enhanced comma separated files (RT, m/z, Intensity, [meta])
          CSV,                # < general comma separated files format (might also be tab or space separated!!!), data should be regular, i.e. matrix form
          TXT,                # < any text format, which has only loose definition of what it actually contains -- thus it is usually hard to say where the file actually came from (e.g. PepNovo).
          OBO,                # < Controlled Vocabulary format
          HTML,               # < any HTML format
          ANALYSISXML,        # < analysisXML format
          XSD,                # < XSD schema format
          PSQ,                # < NCBI binary blast db
          MRM,                # < SpectraST MRM List
          SQMASS,             # < SqLite format for mass and chromatograms
          PQP,                # < OpenSWATH Peptide Query Parameter (PQP) SQLite DB
          MS,                 # < SIRIUS file (.ms)
          OSW,                # < OpenSWATH OpenSWATH report (OSW) SQLite DB
          PSMS,               # < Percolator tab-delimited output (PSM level)
          PIN,                # < Percolator tab-delimited input (PSM level)
          PARAMXML,           # < internal format for writing and reading parameters (also used as part of CTD)
          SPLIB,              # < SpectraST binary spectral library file
          NOVOR,              # < Novor custom parameter file
          XQUESTXML,          # < xQuest XML file format for protein-protein cross-link identifications (.xquest.xml)
          SPECXML,            # < xQuest XML file format for matched spectra (.spec.xml)
          JSON,               # < JavaScript Object Notation file (.json)
          RAW,                # < Thermo Raw File (.raw)
          OMS,                # < OpenMS database file
          EXE,                # < Executable (.exe)
          XML,                # < any XML format
          BZ2,                # < any BZ2 compressed file
          GZ,                 # < any Gzipped file
          PARQUET,            # < Apache Parquet file format (.parquet, .pqt)
          SIZE_OF_TYPE        # < No file type. Simply stores the number of types

    cdef enum class FileProperties "OpenMS::FileTypes::FileProperties":
        # wrap-attach:
        #    FileTypes
        READABLE,                     # < SOMETHING in OpenMS can read this
        WRITEABLE,                    # < SOMETHING in OpenMS can write this
        PROVIDES_SPECTRUM,            # < Provides spectrum data
        PROVIDES_EXPERIMENT,          # < Provides experiment data
        PROVIDES_FEATURES,            # < Provides features
        PROVIDES_CONSENSUSFEATURES,   # < Provides consensus features
        PROVIDES_IDENTIFICATIONS,     # < Provides identifications
        PROVIDES_TRANSITIONS,         # < Provides transitions
        PROVIDES_QUANTIFICATIONS,     # < Provides quantifications
        PROVIDES_TRANSFORMATIONS,     # < Provides transformations
        PROVIDES_QC,                  # < Provides quality control data
        SIZE_OF_FILEPROPERTIES        # < Not a property, just the number of 'em

cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS":

    cdef enum class FilterLayout "OpenMS::FilterLayout":
        COMPACT,    # < make a single item, e.g. 'all readable files (*.mzML *.mzXML);;'
        ONE_BY_ONE, # < list all types individually, e.g. 'mzML files (*.mzML);;mzXML files (*.mzXML);;'
        BOTH        # < combine COMPACT and ONE_BY_ONE
