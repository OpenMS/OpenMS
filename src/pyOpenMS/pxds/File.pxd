from String cimport *
from Param cimport *
from StringList cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/SYSTEM/File.h>" namespace "OpenMS":

    cdef cppclass File:

        @staticmethod
        String getExecutablePath()

        # Method used to test if a @p file exists.
        @staticmethod
        bool exists(String file) except + nogil

        # Return true if the file does not exist or the file is empty
        @staticmethod
        bool empty(String file) except + nogil

        # Removes a file (if it exists).
        @staticmethod
        bool remove(String file) except + nogil

        # Removes the specified directory (absolute path). Returns true if successful.
        @staticmethod
        bool removeDirRecursively(String dir_name) except + nogil

        # Replaces the relative path in the argument with the absolute path.
        @staticmethod
        String absolutePath(String file) except + nogil

        # Returns the basename of the file (without the path).
        @staticmethod
        String basename(String file) except + nogil

        # Returns the path of the file (without the file name).
        @staticmethod
        String path(String file) except + nogil

        # Return true if the file exists and is readable
        @staticmethod
        bool readable(String file) except + nogil

        # Return true if the file is writable
        @staticmethod
        bool writable(String file) except + nogil

        # Return true if the given path specifies a directory
        @staticmethod
        bool isDirectory(String path) except + nogil

        # @brief Looks up the location of the file @p filename
        @staticmethod
        String find(String filename, StringList directories) except + nogil

        # @brief Retrieves a list of files matching @p file_pattern in directory
        #       @p dir (returns filenames without paths unless @p full_path is true)
        @staticmethod
        bool fileList(String dir, String file_pattern, StringList output, bool full_path) except + nogil

        # Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
        @staticmethod
        String getUniqueName() except + nogil

        # Returns the OpenMS data path (environment variable overwrites the default installation path)
        @staticmethod
        String getOpenMSDataPath() except + nogil

        @staticmethod
        String getOpenMSHomePath() except + nogil

        # The current OpenMS temporary data path (for temporary files)
        @staticmethod
        String getTempDirectory() except + nogil

        # The current OpenMS user data path (for result files)
        @staticmethod
        String getUserDirectory() except + nogil

        # get the system's default OpenMS.ini file in the users home directory (&lt except + nogil home&gt except + nogil /OpenMS/OpenMS.ini)
        # or create/repair it if required
        @staticmethod
        Param getSystemParameters() except + nogil

        # uses File::find() to search for a file names @p db_name
        # in the 'id_db_dir' param of the OpenMS system parameters
        # @exception FileNotFound is thrown, if the file is not found
        @staticmethod
        String findDatabase(String db_name) except + nogil

        # Searchs for an executable with the given name.
        @staticmethod
        String findExecutable(String toolName) except + nogil

        @staticmethod
        String getTemporaryFile(const String & alternative_file) except + nogil

        # Resolves a partial file name to a documentation file in the doc-folder.
        @staticmethod
        String findDoc(String filename) except + nogil

        @staticmethod
        bool rename(const String & old_filename, const String & new_filename, bool overwrite_existing, bool verbose) except + nogil

        # bool removeDir(const QString & dir_name) except + nogil
