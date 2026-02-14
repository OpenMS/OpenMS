from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from String cimport *
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from FileTypes cimport *
from Types cimport *
from PeakFileOptions cimport *
from FeatureFileOptions cimport *
from ConsensusMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from ProteinIdentification cimport *
from PeptideIdentificationList cimport *
from MSSpectrum cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS":
        # wrap-doc:
        #  Facilitates file handling by file type recognition
        #  This class provides file type recognition from the file name and
        #  for some types from the file content
        #  It offers a common interface to load MSExperiment data
        #  and allows querying for supported file types
        #  
        #  Usage:
        #
        #  .. code-block:: python
        #  
        #    MSExperiment exp;
        #    FileHandler().loadExperiment("test.mzXML", exp)
        #    FileHandler().loadExperiment("test.mzML", exp)
        #  

    cdef cppclass FileHandler:  # wrap=True
        FileHandler() except + nogil 
        FileHandler(FileHandler) except + nogil  # wrap-ignore

        void loadExperiment(String, MSExperiment &) except + nogil
            # wrap-doc:
            #  Loads a file into an MSExperiment
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param exp: The experiment to load the data into
            #  :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
            #  :param log: Progress logging mode
            #  :param rewrite_source_file: Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from
            #  :param compute_hash: If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile
            #  :return: true if the file could be loaded, false otherwise
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeExperiment(String, MSExperiment) except + nogil
            # wrap-doc:
            #  Stores an MSExperiment to a file\n
            #  
            #  The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used
            #  
            #  
            #  :param filename: The name of the file to store the data in
            #  :param exp: The experiment to store
            #  :param log: Progress logging mode
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadFeatures(String, FeatureMap &) except + nogil 
            # wrap-doc:
            #  Loads a file into a FeatureMap
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param map: The FeatureMap to load the data into
            #  :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
            #  :return: true if the file could be loaded, false otherwise
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeFeatures(String, FeatureMap) except + nogil
            # wrap-doc:
            #  Stores a FeatureMap to a file
            #  
            #  The file type to store the data in is determined by the file name.
            #  
            #  :param filename: The name of the file to store the data in
            #  :param map: The FeatureMap to store
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadConsensusFeatures(String, ConsensusMap &) except + nogil
            # wrap-doc:
            #  Loads a file into a ConsensusMap
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param map: The ConsensusMap to load the data into
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeConsensusFeatures(String, ConsensusMap) except + nogil
            # wrap-doc:
            #  Stores a ConsensusMap to a file
            #  
            #  The file type to store the data in is determined by the file name.
            #  
            #  :param filename: The name of the file to store the data in
            #  :param map: The ConsensusMap to store
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadIdentifications(String, libcpp_vector[ProteinIdentification] &, PeptideIdentificationList &) except + nogil
            # wrap-doc:
            #  Loads an identification file into proteinIdentifications and peptideIdentifications
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param protein_ids: The proteinIdentification vector to load the data into
            #  :param peptide_ids: The peptideIdentification list to load the data into
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeIdentifications(String, libcpp_vector[ProteinIdentification], PeptideIdentificationList) except + nogil
            # wrap-doc:
            #  Stores proteins and peptides into an Identification File
            #  
            #  
            #  :param filename: The file name of the file to write to
            #  :param protein_ids: The proteinIdentification vector to store
            #  :param peptide_ids: The peptideIdentification list to store
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadTransitions(String, TargetedExperiment &) except + nogil
            # wrap-doc:
            #  Loads transitions of a spectral library
            #  
            #  
            #  :param filename: The file name of the file to read
            #  :param library: The TargetedExperiment to load
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeTransitions(String, TargetedExperiment) except + nogil
            # wrap-doc:
            #  Stores transitions of a spectral library
            #  
            #  
            #  :param filename: The file name of the file to write
            #  :param library: The TargetedExperiment to store
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadTransformations(String, TransformationDescription &, bool) except + nogil
            # wrap-doc:
            #  Loads a file into Transformations
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param map: The TransformationDescription to load the data into
            #  :param fit_model: Call fitModel() on the map before returning
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeTransformations(String, TransformationDescription) except + nogil
            # wrap-doc:
            #  Stores Transformations to a file
            #  
            #  
            #  :param filename: The file name of the file to write
            #  :param map: The TransformationDescription to store
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        void loadSpectrum(String, MSSpectrum &) except + nogil
            # wrap-doc:
            #  Loads a single MSSpectrum from a file
            #  
            #  
            #  :param filename: The file name of the file to load
            #  :param spec: The spectrum to load the data into
            #  :raises:
            #    Exception: FileNotFound is thrown if the file could not be opened
            #  :raises:
            #    Exception: ParseError is thrown if an error occurs during parsing

        void storeSpectrum(String, MSSpectrum) except + nogil
            # wrap-doc:
            #  Stores a single MSSpectrum to a file
            #  
            #  
            #  :param filename: The file name of the file to store
            #  :param spec: The spectrum to store the data from
            #  :raises:
            #    Exception: UnableToCreateFile is thrown if the file could not be written

        PeakFileOptions  getOptions() except + nogil  # wrap-doc:Access to the options for loading/storing
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Sets options for loading/storing

        FeatureFileOptions getFeatOptions() except + nogil  # wrap-doc:Access to the feature file options for loading/storing
        void setFeatOptions(FeatureFileOptions) except + nogil  # wrap-doc:Sets feature file options for loading/storing

        @staticmethod
        int getType(const String& filename) except + nogil
            # wrap-doc:
            #  Determines the file type based on the file name and/or content
            #
            #  :param filename: Path to the file
            #  :returns: Integer representation of the file type

        @staticmethod
        FileType getTypeByFileName(const String & filename) except + nogil
            # wrap-doc:
            #  Determines the file type based on the file extension
            #
            #  :param filename: Path to the file
            #  :returns: The file type based on the extension

        @staticmethod
        FileType getTypeByContent(const String & filename) except + nogil
            # wrap-doc:
            #  Determines the file type based on the file content
            #
            #  :param filename: Path to the file
            #  :returns: The file type based on file content analysis

        @staticmethod
        FileType getConsistentOutputfileType(const String & output_filename, const String & requested_type) except + nogil
            # wrap-doc:
            #  Useful function for TOPP tools which have an 'out_type' parameter
            #
            #  Makes sure that the type derived from output_filename and requested_type are consistent
            #  :param output_filename: A full filename (with path) whose type is determined
            #  :param requested_type: A type as string, usually obtained from '-out_type' parameter
            #  :returns: A consistent file type or UNKNOWN upon conflict

        @staticmethod
        String computeFileHash(const String & filename) except + nogil
            # wrap-doc:
            #  Computes a SHA-1 hash of the file content
            #
            #  :param filename: Path to the file
            #  :returns: SHA-1 hash string of the file content

        @staticmethod
        bool isSupported(FileType type_) except + nogil
            # wrap-doc:
            #  Checks whether the given file type is supported
            #
            #  :param type_: The file type to check
            #  :returns: True if the file type is supported

        @staticmethod
        bool hasValidExtension(const String & filename, FileType type_) except + nogil
            # wrap-doc:
            #  Checks whether the file has a valid extension for the given type
            #
            #  :param filename: Path to the file
            #  :param type_: The expected file type
            #  :returns: True if the extension matches the file type

        @staticmethod
        String stripExtension(String file) except + nogil
            # wrap-doc:
            #  Returns the file name without the extension
            #
            #  :param file: Path to the file
            #  :returns: File path without extension

        @staticmethod
        String swapExtension(String filename, FileType new_type) except + nogil
            # wrap-doc:
            #  Removes the current extension (if any) and adds a new one
            #
            #  :param filename: Path to the file
            #  :param new_type: The new file type whose extension should be used
            #  :returns: File path with the new extension
