// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/config.h>
#include <cstdlib>


namespace OpenMS
{
  class Param;
  class TOPPBase;

  /**
    @brief Basic file handling operations.

    @ingroup System
  */
  class OPENMS_DLLAPI File
  {
public:

    friend class TOPPBase;

    /**
      @brief Class representing a temporary directory
    
    */

    class OPENMS_DLLAPI TempDir
    {
    public:
      
      /// Construct temporary folder
      /// If keep_dir is set to true, the folder will not be deleted on destruction of the object.
      TempDir(bool keep_dir = false);

      /// Destroy temporary folder (can be prohibited in Constructor)
      ~TempDir();

      /// delete all means to copy or move a TempDir
      TempDir(const TempDir&) = delete;
      TempDir& operator=(const TempDir&) = delete;
      TempDir(TempDir&&) = delete;
      TempDir& operator=(TempDir&&) = delete;

      /// Return path to temporary folder
      const String& getPath() const;

    private:
      String temp_dir_;
      bool keep_dir_;
    };

    /// Retrieve path of current executable (useful to find other TOPP tools)
    /// The returned path is either just an EMPTY string if the call to system subroutines failed
    /// or the complete path including a trailing "/", to enable usage of this function as
    /// File::getExecutablePath() + "mytool"
    static String getExecutablePath();

    /// Method used to test if a @p file exists.
    static bool exists(const String& file);

    /// Return true if the file does not exist or the file is empty
    static bool empty(const String& file);

    /**
       @brief Rename a file
       
       If @p from and @p to point to the same file (symlinks are resolved),
       no action will be taken and true is returned.
       If the target already exists (and is not identical to the source),
       this function will fail unless @p overwrite_existing is true.
       
       @param from Source filename
       @param to Target filename
       @param overwrite_existing Delete already existing target, before renaming
       @param verbose Print message to OPENMS_LOG_ERROR if something goes wrong.
       @return True on success
    */
    static bool rename(const String& from, const String& to, bool overwrite_existing = true, bool verbose = true);

    /**
       @brief Copy directory recursively
       
       Copies a source directory to a new target directory (recursive).
       If the target directory already exists, files will be added.
       If files from the source already exist in the target, @p option allows for the following behaviour:
       
       OVERWRITE: Overwrite the file in the target directory if it already exists.
       SKIP: Skip the file in the target directory if it already exists.
       CANCEL: Cancel the copy process if file already exists in target directory - return false.

       @param from_dir Source directory
       @param to_dir Target directory
       @param option Specify the copy option (OVERWRITE, SKIP, CANCEL)
       @return True on success
    */
    enum class CopyOptions {OVERWRITE,SKIP,CANCEL};
    static bool copyDirRecursively(const QString &from_dir, const QString &to_dir, File::CopyOptions option = CopyOptions::OVERWRITE);

    /**
      @brief Removes a file (if it exists).

      @return Returns true if the file was successfully deleted (or if it did not exist).
    */
    static bool remove(const String& file);

    /// Removes the subdirectories of the specified directory (absolute path). Returns true if successful.
    static bool removeDirRecursively(const String& dir_name);

    /// Removes the directory and all subdirectories (absolute path).
    static bool removeDir(const QString& dir_name);

    /// Replaces the relative path in the argument with the absolute path.
    static String absolutePath(const String& file);

    /// Returns the basename of the file (without the path).
    static String basename(const String& file);

    /// Returns the path of the file (without the file name).
    static String path(const String& file);

    /// Return true if the file exists and is readable
    static bool readable(const String& file);

    /// Return true if the file is writable
    static bool writable(const String& file);

    /// Return true if the given path specifies a directory
    static bool isDirectory(const String& path);

    /**
      @brief Looks up the location of the file @p filename

      The following locations are checked in this order:
      - the directories in @p directories
      - the directory contained in the environment variable $OPENMS_DATA_PATH
      - the 'share/OpenMS/' directory of the OpenMS install directory

      @exception FileNotFound is thrown, if the file is not found
    */
    static String find(const String& filename, StringList directories = StringList());

    /**
      @brief Retrieves a list of files matching @p file_pattern in directory
             @p dir (returns filenames without paths unless @p full_path is true)

      @return true => there are matching files
    */
    static bool fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path = false);

    /**
      @brief Resolves a partial file name to a documentation file in the doc-folder.

      Using find() to locate the documentation file under OPENMS_DATA_PATH,
      OPENMS_SOURCE_PATH, OPENMS_BINARY_PATH + "/../../doc" (or a variation for
      MacOS packages)

      Will return the String with the full path to the local documentation. If
      this call fails, try the web documentation
      (http://www.openms.de/current_doxygen/) instead.
     
      @param String The doc file name to find.
      @return The full path to the requested file.

      @exception FileNotFound is thrown, if the file is not found
    */
    static String findDoc(const String& filename);
    
    /**
      @brief Returns a string, consisting of date, time, hostname, process id, and a incrementing number. This can be used for temporary files.

      @param include_hostname add hostname into result - potentially a long string
      @return a unique name
    */
    static String getUniqueName(bool include_hostname = true);

    /// Returns the OpenMS data path (environment variable overwrites the default installation path)
    static String getOpenMSDataPath();

    /// Returns the OpenMS home path (environment variable overwrites the default home path)
    static String getOpenMSHomePath();

    /// The current OpenMS temporary data path (for temporary files)
    static String getTempDirectory();

    /// The current OpenMS user data path (for result files)
    /// Tries to set the user directory in following order:
    ///   1. OPENMS_HOME_DIR if environmental variable set
    ///   2. "home_dir" entry in OpenMS.ini
    ///   3. user home directory
    static String getUserDirectory();

    /// get the system's default OpenMS.ini file in the users home directory (&lt;home&gt;/OpenMS/OpenMS.ini)
    /// or create/repair it if required
    /// order:
    ///   1. &lt;OPENMS_HOME_DIR&gt;/OpenMS/OpenMS.ini if environmental variable set
    ///   2. user home directory &lt;home&gt;/OpenMS/OpenMS.ini
    static Param getSystemParameters();

    /// uses File::find() to search for a file names @p db_name
    /// in the 'id_db_dir' param of the OpenMS system parameters
    /// @exception FileNotFound is thrown, if the file is not found
    static String findDatabase(const String& db_name);

    /**
      @brief Extract list of directories from a concatenated string (usually $PATH).

      Depending on platform, the components are split based on ":" (Linux/Mac) or ";" (Windows).
      All paths use the '/' as separator and end in '/'.
      E.g. for 'PATH=/usr/bin:/home/unicorn' the result is {"/usr/bin/", "/home/unicorn/"}
            or 'PATH=c:\temp;c:\Windows' the result is {"c:/temp/", "c:/Windows/"}

      Note: the environment variable is passed as input to enable proper testing (env vars are usually read-only).  
    */
    static StringList getPathLocations(const String& path = std::getenv("PATH"));

    /**
      @brief Searches for an executable with the given name (similar to @em where (Windows) or @em which (Linux/MacOS)

      This function can be used to find the full path+filename to an executable in
      the PATH environment. Only the @em first hit (by order in PATH) is returned.
      If the @p exe_filename has a relative or full path which points to an existing file, PATH information will not be used.
      The function returns true if the filename was found (exists) and false otherwise.
      Note: this does not require the file to have executable permission set (this is not tested)
      The returned content of @p exe_filename is only valid if true is returned.

      @param[in,out] exe_filename The executable to search for.
      @return true if @p exe_filename could be resolved to a full path and it exists
    */
    static bool findExecutable(OpenMS::String& exe_filename);

    /**
      @brief Searches for an executable with the given name.

      @param toolName The executable to search for.
      @exception FileNotFound is thrown, if the tool executable was not found.
    */
    static String findSiblingTOPPExecutable(const String& toolName);

    /**
      @brief Obtain a temporary filename, ensuring automatic deletion upon exit

      The file is not actually created and only deleted at exit if it exists.
      
      However, if 'alternative_file' is given and not empty, no temporary filename
      is created and 'alternative_file' is returned (and not destroyed upon exit).
      This is useful if you have an optional
      output file, which may, or may not be requested, but you need its content regardless,
      e.g. for intermediate plotting with R.
      Thus you can just call this function to get a file which can be used and gets automatically
      destroyed if needed.

      @param alternative_file If this string is not empty, no action is taken and it is used as return value
      @return Full path to a temporary file
    */
    static const String& getTemporaryFile(const String& alternative_file = "");

    /**
      @brief Helper function to test if filenames provided in two StringLists match.

      Passing several InputFilesLists is error-prone as users may provide files in a different order.
      To check for common mistakes this helper function checks:
      - if both file lists have the same length (returns false otherwise)
      - if the content is the same and provided in exactly the same order (returns false otherwise)

      Note: Because workflow systems may assign file names randomly a non-strict comparison mode is enabled by default.      
      Instead of the strict comparison (which returns false if there is a single mismatch), the non-strict comparison mode 
      only returns false if the unique set of filenames match but some positions differ, i.e., only the order has been mixed up.

      @param sl1 First StringList with filenames
      @param sl2 Second StringList with filenames
      @param basename If set to true, only basenames are compared
      @param ignore_extension If set to true, extensions are ignored (e.g., useful to compare spectra filenames to ID filenames)
      @param strict If set to true, no mismatches (respecting basename and ignore_extension parameter) are allowed. 
                    If set to false, only the order is compared if both share the same filenames.
      @return False, if both StringLists are different (respecting the parameters)
    */
    static bool validateMatchingFileNames(const StringList& sl1, const StringList& sl2, bool basename = true, bool ignore_extension = true, bool strict = false);
private:

    /// get defaults for the system's Temp-path, user home directory etc.
    static Param getSystemParameterDefaults_();

    /// Check if the given path is a valid OPENMS_DATA_PATH
    static bool isOpenMSDataPath_(const String& path);

#ifdef OPENMS_WINDOWSPLATFORM
    /**
      @brief Get list of file suffices to try during search on PATH (usually .exe, .bat etc)

      Input could be ".COM;.EXE;.BAT;.CMD;.VBS".
      If the result does not contain at least ".exe", then we assume the environment variable is broken and return a
      fallback, i.e. {".exe", ".bat"}.

      Note: the environment variable is passed as input to enable proper testing (env vars are usually read-only).

    */
    static StringList executableExtensions_(const String& ext = std::getenv("PATHEXT"));
#endif

    /**
      @brief Internal helper class, which holds temporary filenames and deletes these files at program exit
    */
    class TemporaryFiles_
    {
      public:
        TemporaryFiles_();
        /// create a new filename and queue internally for deletion
        const String& newFile();

        ~TemporaryFiles_();
      private:
        TemporaryFiles_(const TemporaryFiles_&) = delete; // copy is forbidden
        TemporaryFiles_& operator=(const TemporaryFiles_&) = delete;
        StringList filenames_;
    };


    /// private list of temporary filenames, which are deleted upon program exit
    static TemporaryFiles_ temporary_files_;
  };
}

