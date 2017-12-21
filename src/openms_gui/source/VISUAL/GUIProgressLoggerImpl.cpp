// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QProgressDialog>
#include <iostream>

namespace OpenMS
{
  /// create new object (needed by Factory)
  ProgressLogger::ProgressLoggerImpl* GUIProgressLoggerImpl::create()
  {
    return new GUIProgressLoggerImpl();
  }

  /// name of the model (needed by Factory)
  const String GUIProgressLoggerImpl::getProductName()
  {
    return "GUI";
  }

  GUIProgressLoggerImpl::GUIProgressLoggerImpl() :
    dlg_(nullptr),
    begin_(0),
    end_(0)
  {
  }

  void GUIProgressLoggerImpl::startProgress(const SignedSize begin, const SignedSize end, const String& label, const int /* current_recursion_depth */) const
  {
    begin_ = begin;
    end_ = end;
    if (!dlg_)
    {
      dlg_ = new QProgressDialog(label.c_str(), QString(), int(begin), int(end));
    }
    dlg_->setWindowTitle(label.c_str());
    dlg_->setWindowModality(Qt::WindowModal);
    dlg_->show();
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
      }
      else
      {
        std::cout << "ProgressLogger warning: 'setValue' called before 'startProgress'!" << std::endl;
      }
    }
  }

  void GUIProgressLoggerImpl::endProgress(const int /* current_recursion_depth */) const
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
