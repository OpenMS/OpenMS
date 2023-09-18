// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <QtCore/QString>

#include <iostream>

using namespace std;

namespace OpenMS
{
  class CMDProgressLoggerImpl :
    public ProgressLogger::ProgressLoggerImpl
  {
public:
    CMDProgressLoggerImpl() :
      stop_watch_()
      
    {
    }

    /// create new object (needed by Factory)
    static ProgressLogger::ProgressLoggerImpl* create()
    {
      return new CMDProgressLoggerImpl();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "CMD";
    }

    void startProgress(const SignedSize begin, const SignedSize end, const String& label, const int current_recursion_depth) const override
    {
      begin_ = begin;
      current_ = begin_;
      end_ = end;
      cout << string(2 * current_recursion_depth, ' ') << "Progress of '" << label << "':" << endl;
      stop_watch_.reset();
      stop_watch_.start();
    }

    void setProgress(const SignedSize value, const int current_recursion_depth) const override
    {
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
        cout << '\r' << string(2 * current_recursion_depth, ' ') << QString::number(float(value - begin_) / float(end_ - begin_) * 100.0, 'f', 2).toStdString()  << " %               ";
        cout << flush;
      }
    }
    SignedSize nextProgress() const override
    {
      #pragma omp atomic
      ++current_;
      return current_;
    }

    void endProgress(const int current_recursion_depth, UInt64 bytes_processed) const override
    {
      stop_watch_.stop();
      String IO_stats;
      if (bytes_processed)
      {
        IO_stats = " @ " + bytesToHumanReadable(bytes_processed / stop_watch_.getClockTime()) + "/s";
      }
      cout << '\r' << string(2 * current_recursion_depth, ' ') << "-- done [took " << StopWatch::toString(stop_watch_.getCPUTime()) << " (CPU), " << StopWatch::toString(stop_watch_.getClockTime()) << " (Wall)" << IO_stats << "] -- " << endl;
    }

private:
    mutable StopWatch stop_watch_;
    mutable SignedSize begin_{0};
    mutable SignedSize end_{0};
    mutable SignedSize current_{0};
  };

  class NoProgressLoggerImpl :
    public ProgressLogger::ProgressLoggerImpl
  {
public:
    /// create new object (needed by Factory)
    static ProgressLogger::ProgressLoggerImpl* create()
    {
      return new NoProgressLoggerImpl();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "NONE";
    }

    void startProgress(const SignedSize /* begin */, const SignedSize /* end */, const String& /* label */, const int /* current_recursion_depth */) const override
    {
    }

    void setProgress(const SignedSize /* value */, const int /* current_recursion_depth */) const override
    {
    }

    SignedSize nextProgress() const override
    {
      return 0;
    }
    
    void endProgress(const int /* current_recursion_depth */, UInt64 /*bytes_processed*/) const override
    {
    }

  };



  void ProgressLogger::ProgressLoggerImpl::registerChildren()
  {
    Factory<ProgressLogger::ProgressLoggerImpl>::registerProduct(CMDProgressLoggerImpl::getProductName(), &CMDProgressLoggerImpl::create);
    // this will only be registered by GUI base app in OpenMS_GUI
    // Factory<ProgressLogger::ProgressLoggerImpl>::registerProduct(GUIProgressLoggerImpl::getProductName(), &GUIProgressLoggerImpl::create);
    Factory<ProgressLogger::ProgressLoggerImpl>::registerProduct(NoProgressLoggerImpl::getProductName(), &NoProgressLoggerImpl::create);
  }

  int ProgressLogger::recursion_depth_ = 0;

  String ProgressLogger::logTypeToFactoryName_(ProgressLogger::LogType type)
  {
    switch (type)
    {
      case NONE:
      {
        return "NONE";
      }
      case CMD:
      {
        return "CMD";
      }
      case GUI:
      {
        return "GUI";
      }
    }

// should never happen but gcc emits a warning/error
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code-return"
    return "";
#pragma clang diagnostic pop
  }

  ProgressLogger::ProgressLogger() :
    type_(NONE),
    last_invoke_()
  {
    current_logger_ = Factory<ProgressLogger::ProgressLoggerImpl>::create(logTypeToFactoryName_(type_));
  }

  ProgressLogger::ProgressLogger(const ProgressLogger& other) :
    type_(other.type_),
    last_invoke_(other.last_invoke_)
  {
    // recreate our logger
    current_logger_ = Factory<ProgressLogger::ProgressLoggerImpl>::create(logTypeToFactoryName_(type_));
  }

  ProgressLogger& ProgressLogger::operator=(const ProgressLogger& other)
  {
    if (&other == this)
    {
      return *this;
    }

    this->last_invoke_ = other.last_invoke_;
    this->type_ = other.type_;

    // we clean our old logger
    delete current_logger_;

    // .. and get a new one
    current_logger_ = Factory<ProgressLogger::ProgressLoggerImpl>::create(logTypeToFactoryName_(type_));

    return *this;
  }

  ProgressLogger::~ProgressLogger()
  {
    delete current_logger_;
  }

  void ProgressLogger::setLogType(LogType type) const
  {
    type_ = type;
    // remove the old logger
    delete current_logger_;

    current_logger_ = Factory<ProgressLogger::ProgressLoggerImpl>::create(logTypeToFactoryName_(type_));
  }

  ProgressLogger::LogType ProgressLogger::getLogType() const
  {
    return type_;
  }

  void ProgressLogger::startProgress(SignedSize begin, SignedSize end, const String& label) const
  {
    OPENMS_PRECONDITION(begin <= end, "ProgressLogger::init : invalid range!");
    last_invoke_ = time(nullptr);
    current_logger_->startProgress(begin, end, label, recursion_depth_);
    ++recursion_depth_;
  }

  void ProgressLogger::setProgress(SignedSize value) const
  {
    // update only if at least 1 second has passed
    if (last_invoke_ == time(nullptr))
    {
      return;
    }
    last_invoke_ = time(nullptr);
    current_logger_->setProgress(value, recursion_depth_);
  }
  void ProgressLogger::nextProgress() const
  {
    auto p = current_logger_->nextProgress();
    // update only if at least 1 second has passed
    if (last_invoke_ == time(nullptr))
    {
      return;
    }

    last_invoke_ = time(nullptr);
    current_logger_->setProgress(p, recursion_depth_);
  }

  void ProgressLogger::endProgress(UInt64 bytes_processed) const
  {
    if (recursion_depth_)
    {
      --recursion_depth_;
    }
    current_logger_->endProgress(recursion_depth_, bytes_processed);
  }


} //namespace OpenMS
