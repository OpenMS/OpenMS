// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QApplication>
#include <QProgressDialog>
#include <iostream>

namespace OpenMS
{
  GUIProgressLoggerImpl::GUIProgressLoggerImpl() :
    dlg_(nullptr),
    begin_(0),
    end_(0),
    current_(0)
  {
  }

  void GUIProgressLoggerImpl::startProgress(const SignedSize begin, const SignedSize end, const String& label, const int /* current_recursion_depth */) const
  {
    begin_ = begin;
    current_ = begin_;
    end_ = end;
    delete dlg_; // delete old dialog, if present
    dlg_ = new QProgressDialog(label.c_str(), QString(), int(begin), int(end));
    dlg_->setWindowTitle(label.c_str());
    dlg_->setWindowModality(Qt::WindowModal);
    dlg_->show();
    QApplication::processEvents(); // show it...
  }

  void GUIProgressLoggerImpl::setProgress(const SignedSize value, const int /* current_recursion_depth */) const
  {
    if (value < begin_ || value > end_)
    {
      std::cout << "ProgressLogger: Invalid progress value '" << value << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << std::endl;
    }
    else
    {
      if (dlg_)
      {
        dlg_->setValue((int)value);
        QApplication::processEvents(); // show it...
      }
      else
      {
        std::cout << "ProgressLogger warning: 'setProgress' called before 'startProgress'!" << std::endl;
      }
    }
  }
  SignedSize GUIProgressLoggerImpl::nextProgress() const
  {
    return ++current_;
  }

  void GUIProgressLoggerImpl::endProgress(const int /* current_recursion_depth */, UInt64 /* bytes_processed */) const
  {
    if (dlg_)
    {
      dlg_->setValue((int)end_);
    }
    else
    {
      std::cout << "ProgressLogger warning: 'endProgress' called before 'startProgress'!" << std::endl;
    }
  }

  GUIProgressLoggerImpl::~GUIProgressLoggerImpl()
  {
    delete dlg_;
  }
}
