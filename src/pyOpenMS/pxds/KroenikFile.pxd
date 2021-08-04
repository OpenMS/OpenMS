from String cimport *

from MSSpectrum cimport *
from Peak1D cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/FORMAT/KroenikFile.h>" namespace "OpenMS":

    cdef cppclass KroenikFile:
        # wrap-doc:
                #   File adapter for Kroenik (HardKloer sibling) files
                #   -----
                #   The first line is the header and contains the column names:
                #   File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
                #   -----
                #   Every subsequent line is a feature
                #   -----
                #   All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
                #   metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications"
                #   -----
                #   The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file

        KroenikFile() nogil except +

        void store(String filename, MSSpectrum & spectrum)  nogil except + # wrap-doc:Stores a MSExperiment into a Kroenik file
        void load(String filename, FeatureMap & feature_map) nogil except +
            # wrap-doc:
                #   Loads a Kroenik file into a featureXML
                #   -----
                #   The content of the file is stored in `features`
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be opened
                #   :raises:
                #     Exception: ParseError is thrown if an error occurs during parsing
