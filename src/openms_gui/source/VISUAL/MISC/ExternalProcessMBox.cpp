// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MISC/ExternalProcessMBox.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
#include <QMessageBox>
#include <QCoreApplication>
#include <utility>

namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcessMBox::ExternalProcessMBox() = default;

  ExternalProcessMBox::ExternalProcessMBox(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
    : ep_(std::move(callbackStdOut), std::move(callbackStdErr))
  {
  }

  ExternalProcessMBox::~ExternalProcessMBox() = default;

  /// re-wire the callbacks used using run()
  void ExternalProcessMBox::setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
  {
    ep_.setCallbacks(std::move(callbackStdOut), std::move(callbackStdErr));
  }

  ExternalProcess::RETURNSTATE ExternalProcessMBox::run(QWidget* parent, const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose, String& error_msg)
  {
    // Pass idle_callback that pumps the GUI event loop to prevent freezing
    auto rs = ep_.run(exe, args, working_dir, verbose, error_msg, ExternalProcess::IO_MODE::READ_WRITE, {},
                      []() { QCoreApplication::processEvents(); });

    if (!error_msg.empty())
    {
      QMessageBox::critical(parent, "Error", toQString(error_msg));
    }

    return rs;
  }

  ExternalProcess::RETURNSTATE ExternalProcessMBox::run(QWidget* parent, const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose)
  {
    String error_msg;
    auto rs = ep_.run(exe, args, working_dir, verbose, error_msg, ExternalProcess::IO_MODE::READ_WRITE, {},
                      []() { QCoreApplication::processEvents(); });

    if (!error_msg.empty()) QMessageBox::critical(parent, "Error", toQString(error_msg));

    return rs;
  }

} // ns OpenMS
