from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from String cimport *
from libcpp.string cimport string as libcpp_string
from FileTypes cimport *
from Types cimport *
from PeakFileOptions cimport *

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

        PeakFileOptions  getOptions() except + nogil  # wrap-doc:Access to the options for loading/storing
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Sets options for loading/storing

#
# wrap static method:
#
cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS::FileHandler":

    int getType(const String& filename) except + nogil  # wrap-attach:FileHandler
    FileType getTypeByFileName(const String & filename) except + nogil  # wrap-attach:FileHandler 
    FileType getTypeByContent(const String & filename) except + nogil  # wrap-attach:FileHandler 
    String computeFileHash(const String & filename) except + nogil  # wrap-attach:FileHandler 
    bool isSupported(FileType type_) except + nogil  # wrap-attach:FileHandler 
    bool hasValidExtension(const String & filename, FileType type_) except + nogil  # wrap-attach:FileHandler 

    # Returns the file name without the extension
    String stripExtension(String file) except + nogil  # wrap-attach:FileHandler
    # Removes the current extension (if any) and adds a new one
    String swapExtension(String filename, FileType new_type) except + nogil  # wrap-attach:FileHandler 
