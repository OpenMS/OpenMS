// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_SYSTEM_FILE_H
#define OPENMS_SYSTEM_FILE_H

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/config.h>


namespace OpenMS
{

  class TOPPBase;

  /**
    @brief Basic file handling operations.

    @ingroup System
  */
  class OPENMS_DLLAPI File
  {
public:

    friend class TOPPBase;

    /// Retrieve path of current executable (useful to find other TOPP tools)
    /// The returned path is either just an EMPTY string if the call to system subroutines failed
    /// or the complete path including a trailing "/", to enable usage of this function as
    ///  File::getExecutablePath() + "mytool"
    static String getExecutablePath();

    /// Method used to test if a @p file exists.
    static bool exists(const String& file);

    /// Return true if the file does not exist or the file is empty
    static bool empty(const String& file);

    /**
      @brief Removes a file (if it exists).

      @return Returns true if the file was successfully deleted (or if it did not exist).
    */
    static bool remove(const String& file);

    /// Removes the specified directory (absolute path). Returns true if successful.
    static bool removeDirRecursively(const String& dir_name);

    /// Replaces the relative path in the argument with the absolute path.
    static String absolutePath(const String& file);

    /// Returns the basename of the file (without the path).
    static String basename(const String& file);

    /// Returns the path of the file (without the file name).
    static String path(const String& file);

    /**
      Returns the file name without the extension

      The extension is the suffix of the string upto and including the last dot.

      If no extension is found, the whole file name is returned
    */
    static String removeExtension(const String& file);

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

    /// Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
    static String getUniqueName();

    /// Returns the OpenMS data path (environment variable overwrites the default installation path)
    static String getOpenMSDataPath();

    /// The current OpenMS temporary data path (for temporary files)
    static String getTempDirectory();

    /// The current OpenMS user data path (for result files)
    static String getUserDirectory();

    /// get the system's default OpenMS.ini file in the users home directory (&lt;home&gt;/OpenMS/OpenMS.ini)
    /// or create/repair it if required
    static Param getSystemParameters();

    /// uses File::find() to search for a file names @p db_name
    /// in the 'id_db_dir' param of the OpenMS system parameters
    /// @exception FileNotFound is thrown, if the file is not found
    static String findDatabase(const String& db_name);

    /**
      @brief Searchs for an executable with the given name.

      @param toolName The executable to search for.
      @exception FileNotFound is thrown, if the tool executable was not found.
    */
    static String findExecutable(const String& toolName);

private:

    /// get defaults for the system's Temp-path, user home directory etc
    static Param getSystemParameterDefaults_();

    /// Check if the given path is a valid OPENMS_DATA_PATH
    static bool isOpenMSDataPath_(const String& path);
  };

}

#endif // OPENMS_SYSTEM_FILE_H
