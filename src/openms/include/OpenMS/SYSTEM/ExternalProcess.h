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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QObject>

#include <functional> // for std::function

class QProcess; // forward declare to avoid header include
class QString;
class QStringList;

namespace OpenMS
{

  /**
    @brief A wrapper around QProcess to conveniently start an external program and forward its outputs

    Use the custom Ctor to provide callback functions for stdout/stderr output or set them via setCallbacks().

    Running an external program blocks the caller, so do not use this in a main GUI thread
    (unless you have some other means to tell the user that no interaction is possible at the moment).

    @Note If you want QMessageboxes to be shown if something went wrong, use ExternalProcessMBox as a convenient wrapper instead.

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

    /// default Ctor; callbacks for stdout/stderr are empty
    ExternalProcess();

    /// set the callback functions to process stdout and stderr output when the external process generates it
    ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /// D'tor
    ~ExternalProcess();

    /// re-wire the callbacks used using run()
    void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr);

    /**
      @brief Runs a program and calls the callback functions from time to time if output from the external program is available.

      @param exe The program to call (can contain spaces in path, no problem)
      @param args A list of extra arguments (can be empty)
      @param verbose Report the call command and errors via the callbacks (default: false)
      @param working_dir Execute the external process in the given directory (relevant when relative input/output paths are given). Leave empty to use the current working directory.
      @param[out] error_string Message to display to the user if something went wrong (if return != SUCCESS)
      @return Did the external program succeed (SUCCESS) or did something go wrong?
    */
    RETURNSTATE run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg);
    
    /**
      @brief Same as other overload, just without a returned error message
     */
    ExternalProcess::RETURNSTATE run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose);

  private slots:
    void processStdOut_();
    void processStdErr_();

  private:
    QProcess* qp_; ///< pointer to avoid including the QProcess header here (it's huge)
    std::function<void(const String&)> callbackStdOut_;
    std::function<void(const String&)> callbackStdErr_;
  };
} // ns OpenMS
