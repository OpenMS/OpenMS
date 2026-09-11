// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <mutex>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief A uniquely named temporary directory, removed when the object is destroyed.

    The directory is created below SystemSettings::getTempDirectory() unless an
    explicit base directory is given. Destruction removes it recursively unless
    @p keep_dir was set.

    @ingroup System
  */
  class OPENMS_DLLAPI TempDir
  {
  public:
    /// Construct temporary folder under system temp directory
    /// If keep_dir is set to true, the folder will not be deleted on destruction of the object.
    TempDir(bool keep_dir = false);

    /// Construct temporary folder under a custom base directory
    /// Creates a unique subdirectory with a generated name under base_dir.
    /// If keep_dir is set to true, the folder will not be deleted on destruction of the object.
    /// @param base_dir The base directory under which to create the temp folder (e.g., user-specified temp path)
    /// @param keep_dir If true, the folder will not be deleted on destruction
    TempDir(const std::string& base_dir, bool keep_dir = false);

    /// Destroy temporary folder (can be prohibited in Constructor)
    ~TempDir();

    /// delete all means to copy or move a TempDir
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&&) = delete;
    TempDir& operator=(TempDir&&) = delete;

    /// Return path to temporary folder
    const std::string& getPath() const;

  private:
    std::string temp_dir_;
    bool keep_dir_;
  };

  /**
    @brief Temporary files that live until the process exits.

    Names are allocated below SystemSettings::getTempDirectory() and every file
    that still exists at exit is removed.

    @ingroup System
  */
  class OPENMS_DLLAPI TempFiles
  {
  public:
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

      @param[in] alternative_file If this string is not empty, no action is taken and it is used as return value
      @return Full path to a temporary file
    */
    static std::string getTemporaryFile(const std::string& alternative_file = "");

  private:
    /// Holds the allocated filenames and deletes the files at program exit
    class Registry_
    {
    public:
      Registry_(const Registry_&) = delete; // copy is forbidden
      Registry_& operator=(const Registry_&) = delete;
      Registry_();
      /// create a new filename and queue internally for deletion
      std::string newFile();

      ~Registry_();
    private:
      std::vector<std::string> filenames_;
      std::mutex mtx_;
    };

    /// filenames handed out so far, deleted upon program exit
    static Registry_ registry_;
  };
}
