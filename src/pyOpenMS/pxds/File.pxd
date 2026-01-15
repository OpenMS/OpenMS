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
    bool exists(String file) except + nogil  # wrap-attach:File

    # Return true if the file does not exist or the file is empty
    bool empty(String file) except + nogil  # wrap-attach:File

    # Removes a file (if it exists).
    bool remove(String file) except + nogil  # wrap-attach:File

    # Removes the specified directory (absolute path). Returns true if successful.
    bool removeDirRecursively(String dir_name) except + nogil  # wrap-attach:File

    # Replaces the relative path in the argument with the absolute path.
    String absolutePath(String file) except + nogil  # wrap-attach:File

    # Returns the basename of the file (without the path).
    String basename(String file) except + nogil  # wrap-attach:File

    # Returns the path of the file (without the file name).
    String path(String file) except + nogil  # wrap-attach:File

    # Return true if the file exists and is readable
    bool readable(String file) except + nogil  # wrap-attach:File

    # Return true if the file is writable
    bool writable(String file) except + nogil  # wrap-attach:File

    # Return true if the given path specifies a directory
    bool isDirectory(String path) except + nogil  # wrap-attach:File

    # @brief Looks up the location of the file @p filename
    String find(String filename, StringList directories) except + nogil  # wrap-attach:File

    # @brief Retrieves a list of files matching @p file_pattern in directory
    #       @p dir (returns filenames without paths unless @p full_path is true)
    bool fileList(String dir, String file_pattern, StringList output, bool full_path) except + nogil  # wrap-attach:File

    # Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
    String getUniqueName() except + nogil  # wrap-attach:File

    # Returns the OpenMS data path (environment variable overwrites the default installation path)
    String getOpenMSDataPath() except + nogil  # wrap-attach:File

    String getOpenMSHomePath() except + nogil  # wrap-attach:File

    # The current OpenMS temporary data path (for temporary files)
    String getTempDirectory() except + nogil  # wrap-attach:File

    # The current OpenMS user data path (for result files)
    String getUserDirectory() except + nogil  # wrap-attach:File

    # get the system's default OpenMS.ini file in the users home directory (&lt except + nogil home&gt except + nogil /OpenMS/OpenMS.ini)
    # or create/repair it if required
    Param getSystemParameters() except + nogil  # wrap-attach:File

    # uses File::find() to search for a file names @p db_name
    # in the 'id_db_dir' param of the OpenMS system parameters
    # @exception FileNotFound is thrown, if the file is not found
    String findDatabase(String db_name) except + nogil  # wrap-attach:File

    # Searchs for an executable with the given name.
    String findExecutable(String toolName) except + nogil  # wrap-attach:File

    String getTemporaryFile(const String & alternative_file) except + nogil  # wrap-attach:File

    # Resolves a partial file name to a documentation file in the doc-folder.
    String findDoc(String filename) except + nogil  # wrap-attach:File

    bool rename(const String & old_filename, const String & new_filename, bool overwrite_existing, bool verbose) except + nogil  # wrap-attach:File

    # bool removeDir(const QString & dir_name) except + nogil 

