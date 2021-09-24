from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from FeatureMap cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/MzQCFile.h>" namespace "OpenMS":

    cdef cppclass MzQCFile "OpenMS::MzQCFile":
        # wrap-doc:
        #   File adapter for mzQC files used to load and store mzQC files
        #   -----
        #   This class collects the data for the mzQC File

        MzQCFile() nogil except +
        
        MzQCFile(MzQCFile) nogil except + #wrap-ignore

        void store(String input_file,
                   String output_file,
                   MSExperiment & exp,
                   String contact_name,
                   String contact_address,
                   String description,
                   String label,
                   FeatureMap & feature_map,
                   libcpp_vector[ ProteinIdentification ] & prot_ids, 
                   libcpp_vector[ PeptideIdentification ] & pep_ids) nogil except +
                   # wrap-doc:
                   #   Stores QC data in mzQC file with JSON format
                   #   :param input_file: MzML input file name
                   #   :param output_file: MzQC output file name
                   #   :param exp: MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
                   #   :param contact_name: Name of the person creating the mzQC file
                   #   :param contact_address: Contact address (mail/e-mail or phone) of the person creating the mzQC file
                   #   :param description: Description and comments about the mzQC file contents
                   #   :param label: Qnique and informative label for the run
                   #   :param feature_map: FeatureMap from feature file (featureXML)
                   #   :param prot_ids: Protein identifications from ID file (idXML)
                   #   :param pep_ids: Protein identifications from ID file (idXML)
