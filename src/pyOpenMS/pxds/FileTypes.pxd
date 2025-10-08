from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS":

    cdef cppclass FileTypes:

        # compiler
        FileTypes() except + nogil  # wrap-doc:Centralizes the file types recognized by FileHandler
        FileTypes(FileTypes &) except + nogil  # compiler

        String typeToName(FileType t) except + nogil  # wrap-doc:Returns the name/extension of the type

        String typeToMZML(FileType t) except + nogil  # wrap-doc:Returns the mzML name

        FileType nameToType(String name) except + nogil  
            # wrap-doc:
                #  Converts a file type name into a Type 
                #  
                #  
                #  :param name: A case-insensitive name (e.g. FASTA or Fasta, etc.)

cdef extern from "<OpenMS/FORMAT/FileTypes.h>" namespace "OpenMS::FileTypes":

    cdef enum FileType "OpenMS::FileTypes::Type":
    
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
          GELML,              # < GelML (HUPO PSI format) (.gelML)
          TRAML,              # < TraML (HUPO PSI format) for transitions (.traML)
          MSP,                # < NIST spectra library file format (.msp)
          OMSSAXML,           # < OMSSA XML file format for peptide identifications (.xml)
          MASCOTXML,          # < Mascot XML file format for peptide identifications (.xml)
          PNG,                # < Portable Network Graphics (.png)
          XMASS,              # < XMass Analysis file (fid)
          TSV,                # < msInspect file (.tsv)
          PEPLIST,            # < specArray file (.peplist)
          HARDKLOER,          # < hardkloer file (.hardkloer)
          KROENIK,            # < kroenik file (.kroenik)
          FASTA,              # < FASTA file (.fasta)
          EDTA,               # < enhanced comma separated files (RT, m/z, Intensity, [meta])
          CSV,                # < general comma separated files format (might also be tab or space separated!!!), data should be regular, i.e. matrix form
          TXT,                # < any text format, which has only loose definition of what it actually contains -- thus it is usually hard to say where the file actually came from (e.g. PepNovo).
          OBO,                # < Controlled Vocabulary format
          HTML,               # < any HTML format
          XML,                # < any XML format
          ANALYSISXML,        # < analysisXML format
          XSD,                # < XSD schema format
          PSQ,                # < NCBI binary blast db
          MRM,                # < SpectraST MRM List
          SQMASS,             # < SqLite format for mass and chromatograms
          PQP,                # < OpenSWATH Peptide Query Parameter (PQP) SQLite DB
          OSW,                # < OpenSWATH OpenSWATH report (OSW) SQLite DB
          PSMS,               # < Percolator tab-delimited output (PSM level)
          PARAMXML,           # < internal format for writing and reading parameters (also used as part of CTD)
          SIZE_OF_TYPE        # < No file type. Simply stores the number of types
