// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QObject>

#include <functional> // for std::function

class QProcess; // forward declare to avoid header include
class QString;
#include <QtCore/qcontainerfwd.h> // for QStringList

namespace OpenMS
{

  /**
    @brief A wrapper around QProcess to conveniently start an external program and forward its outputs

    Use the custom Ctor to provide callback functions for stdout/stderr output or set them via setCallbacks().

    Running an external program blocks the caller, so do not use this in a main GUI thread
    (unless you have some other means to tell the user that no interaction is possible at the moment).

    If you want QMessageBoxes to be shown if something went wrong, use ExternalProcessMBox as a convenient wrapper instead.

  */
  class OPENMS_DLLAPI ExternalProcess
    : public QObject
  {
    Q_OBJECT

  public:
    /// result of calling an external executable
    enum class RETURNSTATE
    {
      SUCCESS,  ///< everything went smoothly (exit-code = 0)
      NONZERO_EXIT, /// finished, but returned with an exit-code other than 0
      CRASH, ///< ran, but crashed (segfault etc)
      FAILED_TO_START ///< executable not found or not enough access rights for user
    };

    /// Open mode for the process.
    enum class IO_MODE
    {
        NO_IO, ///< No read nor write access
        READ_ONLY,
        WRITE_ONLY,
        READ_WRITE
    };

    /// default Ctor; callbacks for stdout/stderr are empty
    ExternalProcess();

    /// set the callback functions to process stdout and stderr output when the external process generates it
    ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /// D'tor
    ~ExternalProcess() override ;

    /// re-wire the callbacks used during run()
    void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /**
      @brief Runs a program and calls the callback functions from time to time if output from the external program is available.

      @param exe The program to call (can contain spaces in path, no problem)
      @param args A list of extra arguments (can be empty)
      @param working_dir Execute the external process in the given directory (relevant when relative input/output paths are given). Leave empty to use the current working directory.
      @param verbose Report the call command and errors via the callbacks (default: false)
      @param[out] error_msg Message to display to the user if something went wrong (if return != SUCCESS)
      @param io_mode Open mode for the process (read access, write access, ...)
      @return Did the external program succeed (SUCCESS) or did something go wrong?
    */
    RETURNSTATE run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg, IO_MODE io_mode = IO_MODE::READ_WRITE);
    
    /**
      @brief Same as other overload, just without a returned error message
     */
    ExternalProcess::RETURNSTATE run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, IO_MODE io_mode = IO_MODE::READ_WRITE);

  private slots:
    void processStdOut_();
    void processStdErr_();

  private:
    QProcess* qp_; ///< pointer to avoid including the QProcess header here (it's huge)
    std::function<void(const String&)> callbackStdOut_;
    std::function<void(const String&)> callbackStdErr_;
  };
} // ns OpenMS
