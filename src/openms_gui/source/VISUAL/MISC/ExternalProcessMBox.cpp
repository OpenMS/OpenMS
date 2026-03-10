// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MISC/ExternalProcessMBox.h>
#include <QMessageBox>
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

  ExternalProcess::RETURNSTATE ExternalProcessMBox::run(QWidget* parent, const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg)
  {
    const std::string exe_std = exe.toStdString();
    std::vector<std::string> args_std;
    args_std.reserve(static_cast<size_t>(args.size()));
    for (const auto& a : args) args_std.emplace_back(a.toStdString());
    const std::string wd_std = working_dir.toStdString();

    auto rs = ep_.run(exe_std, args_std, wd_std, verbose, error_msg);

    if (!error_msg.empty())
    {
      QMessageBox::critical(parent, "Error", QString::fromStdString(static_cast<const std::string&>(error_msg)));
    }

    return rs;
  }

  ExternalProcess::RETURNSTATE ExternalProcessMBox::run(QWidget* parent, const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose)
  {
    String error_msg;
    const std::string exe_std = exe.toStdString();
    std::vector<std::string> args_std;
    args_std.reserve(static_cast<size_t>(args.size()));
    for (const auto& a : args) args_std.emplace_back(a.toStdString());
    const std::string wd_std = working_dir.toStdString();

    auto rs = ep_.run(exe_std, args_std, wd_std, verbose, error_msg);

    if (!error_msg.empty())
    {
      QMessageBox::critical(parent, "Error", QString::fromStdString(static_cast<const std::string&>(error_msg)));
    }

    return rs;
  }

} // ns OpenMS
