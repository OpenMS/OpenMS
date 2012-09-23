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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QString>
#include <QtGui/QProgressDialog>

#include <iostream>

using namespace std;

namespace OpenMS
{
  int ProgressLogger::recursion_depth_ = 0;

  ProgressLogger::ProgressLogger() :
    type_(NONE),
    begin_(0),
    end_(0),
    value_(0),
    dlg_(0),
    stop_watch_(),
    last_invoke_()
  {
  }

  ProgressLogger::~ProgressLogger()
  {
    delete(dlg_);
  }

  void ProgressLogger::setLogType(LogType type) const
  {
    type_ = type;
  }

  ProgressLogger::LogType ProgressLogger::getLogType() const
  {
    return type_;
  }

  void ProgressLogger::startProgress(SignedSize begin, SignedSize end, const String & label) const
  {
    OPENMS_PRECONDITION(begin <= end, "ProgressLogger::init : invalid range!");
    last_invoke_ = time(NULL);

    switch (type_)
    {
    case CMD:
      begin_ = begin;
      end_ = end;
      if (recursion_depth_)
        cout << '\n';
      cout << string(2 * recursion_depth_, ' ') << "Progress of '" << label << "':" << endl;
      stop_watch_.reset();
      stop_watch_.start();
      break;

    case GUI:
      begin_ = begin;
      end_ = end;
      if (!dlg_)
        dlg_ = new QProgressDialog(label.c_str(), QString(), int(begin), int(end));
      dlg_->setWindowTitle(label.c_str());
      dlg_->setWindowModality(Qt::WindowModal);
      dlg_->show();
      break;

    case NONE:
      break;
    }
    ++recursion_depth_;
    return;
  }

  void ProgressLogger::setProgress(SignedSize value) const
  {
    // update only if at least 1 second has passed
    if (last_invoke_ == time(NULL))
      return;

    last_invoke_ = time(NULL);

    switch (type_)
    {
    case CMD:
      if (begin_ == end_)
      {
        cout << '.' << flush;
      }
      else if (value < begin_ || value > end_)
      {
        cout << "ProgressLogger: Invalid progress value '" << value
             << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << endl;
      }
      else
      {
        cout << '\r' << string(2 * recursion_depth_, ' ') << QString::number(Real(value - begin_) / Real(end_ - begin_) * 100.0, 'f', 2).toStdString()  << " %               ";
        cout << flush;
      }
      break;

    case GUI:
      if (value < begin_ || value > end_)
      {
        cout << "ProgressLogger: Invalid progress value '" << value << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << endl;
      }
      else
      {
        if (dlg_)
        {
          dlg_->setValue((int)value);
        }
        else
        {
          cout << "ProgressLogger warning: 'setValue' called before 'startProgress'!" << endl;
        }
      }
      break;

    case NONE:
      break;
    }
  }

  void ProgressLogger::endProgress() const
  {
    if (recursion_depth_)
      --recursion_depth_;
    switch (type_)
    {
    case CMD:
      stop_watch_.stop();
      if (begin_ == end_)
      {
        if (recursion_depth_)
          cout << '\n';
        cout << endl << string(2 * recursion_depth_, ' ') << "-- done [took " << String::number(stop_watch_.getCPUTime(), 3) << " s(CPU), " << String::number(stop_watch_.getClockTime(), 3) << " s(Wall)] -- " << endl;
      }
      else
      {
        cout << '\r' << string(2 * recursion_depth_, ' ') << "-- done [took " << String::number(stop_watch_.getCPUTime(), 3) << " s(CPU), " << String::number(stop_watch_.getClockTime(), 3) << " s(Wall)] -- " << endl;
      }
      break;

    case GUI:
      if (dlg_)
      {
        dlg_->setValue((int)end_);
      }
      else
      {
        cout << "ProgressLogger warning: 'endProgress' called before 'startProgress'!" << endl;
      }
      break;

    case NONE:
      break;
    }
  }

} //namespace OpenMS
