// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/SYSTEM/ExternalProcess.h>

namespace OpenMS
{

  /**
    @brief A wrapper around ExternalProcess to conveniently show a MessageBox when an error occurs

    Use the custom Ctor to provide callback functions for stdout/stderr output or set them via setCallbacks().

    Running an external program blocks the caller, so do not use this in a main GUI thread
    (unless you have some other means to tell the user that no interaction is possible at the moment).

  */
  class OPENMS_GUI_DLLAPI ExternalProcessMBox
  {
  public:
    /// default Ctor; callbacks for stdout/stderr are empty
    ExternalProcessMBox();

    /// set the callback functions to process stdout and stderr output when the external process generates it
    ExternalProcessMBox(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /// D'tor
    ~ExternalProcessMBox();

    /// re-wire the callbacks used using run()
    void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /**
      @brief Runs a program by calling ExternalProcess::run and shows any error reported in @p error_msg as a MessageBox before this function returns

      @param parent Optional parent widget, used to position QMesssageBoxes above the parent
      @param exe The program to call (can contain spaces in path, no problem)
      @param args A list of extra arguments (can be empty)
      @param working_dir Execute the external process in the given directory (relevant when relative input/output paths are given). Leave empty to use the current working directory.
      @param verbose Report the call command and errors via the callbacks (default: false)
      @param[out] error_msg Message to display to the user or log somewhere if something went wrong (if return != SUCCESS)
      @return Did the external program succeed (SUCCESS) or did something go wrong?
    */
    ExternalProcess::RETURNSTATE run(QWidget* parent, const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg);

    /**
      @brief Same as other overload, just without a returned error message
    */
    ExternalProcess::RETURNSTATE run(QWidget* parent, const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose);


  private:
    ExternalProcess ep_;
  };
} // ns OpenMS
