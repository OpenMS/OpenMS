from String cimport *
from Param cimport *
from StringList cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/SYSTEM/File.h>" namespace "OpenMS":

    cdef cppclass File:
        pass

# File has only methods, which we wrap as declared below:
cdef extern from "<OpenMS/SYSTEM/File.h>" namespace "OpenMS::File":

    String getExecutablePath()  # wrap-attach:File

    # Method used to test if a @p file exists.
    bool exists(String file_) nogil except + # wrap-attach:File

    # Return true if the file does not exist or the file is empty
    bool empty(String file_) nogil except + # wrap-attach:File

    # Removes a file (if it exists).
    bool remove(String file_) nogil except + # wrap-attach:File

    # Removes the specified directory (absolute path). Returns true if successful.
    bool removeDirRecursively(String dir_name) nogil except + # wrap-attach:File

    # Replaces the relative path in the argument with the absolute path.
    String absolutePath(String file) nogil except + # wrap-attach:File

    # Returns the basename of the file (without the path).
    String basename(String file) nogil except + # wrap-attach:File

    # Returns the path of the file (without the file name).
    String path(String file) nogil except + # wrap-attach:File

    # Returns the file name without the extension
    String removeExtension(String file) nogil except + # wrap-attach:File

    # Return true if the file exists and is readable
    bool readable(String file) nogil except + # wrap-attach:File

    # Return true if the file is writable
    bool writable(String file) nogil except + # wrap-attach:File

    # Return true if the given path specifies a directory
    bool isDirectory(String path) nogil except + # wrap-attach:File

    # @brief Looks up the location of the file @p filename
    String find(String filename, StringList directories) nogil except + # wrap-attach:File

    # @brief Retrieves a list of files matching @p file_pattern in directory
    #        @p dir (returns filenames without paths unless @p full_path is true)
    bool fileList(String dir, String file_pattern, StringList output, bool full_path) nogil except + # wrap-attach:File

    # Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
    String getUniqueName() nogil except + # wrap-attach:File

    # Returns the OpenMS data path (environment variable overwrites the default installation path)
    String getOpenMSDataPath() nogil except + # wrap-attach:File

    # The current OpenMS temporary data path (for temporary files)
    String getTempDirectory() nogil except + # wrap-attach:File

    # The current OpenMS user data path (for result files)
    String getUserDirectory() nogil except + # wrap-attach:File

    # get the system's default OpenMS.ini file in the users home directory (&lt nogil except +home&gt nogil except +/OpenMS/OpenMS.ini)
    # or create/repair it if required
    Param getSystemParameters() nogil except + # wrap-attach:File

    # uses File::find() to search for a file names @p db_name
    # in the 'id_db_dir' param of the OpenMS system parameters
    # @exception FileNotFound is thrown, if the file is not found
    String findDatabase(String db_name) nogil except + # wrap-attach:File

    # Searchs for an executable with the given name.
    String findExecutable(String toolName) nogil except + # wrap-attach:File


