// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>

class QProgressDialog;

namespace OpenMS
{
  /**
    @brief Implements a GUI version of the ProgressLoggerImpl.
  */
  class OPENMS_GUI_DLLAPI GUIProgressLoggerImpl :
    public ProgressLogger::ProgressLoggerImpl
  {
public:
    /// default c'tor.
    GUIProgressLoggerImpl();

    /**
      @brief Implement ProgressLoggerImpl::startProgress().
    */
    void startProgress(const SignedSize begin, const SignedSize end, const String& label, const int /* current_recursion_depth */) const override;

    /**
      @brief Implement ProgressLoggerImpl::setProgress().
    */
    void setProgress(const SignedSize value, const int /* current_recursion_depth */) const override;
    
    /**
      @brief Implement ProgressLoggerImpl::nextProgress().
    */
    SignedSize nextProgress() const override;
    
    /**
      @brief Implement ProgressLoggerImpl::endProgress().
    */
    void endProgress(const int current_recursion_depth, UInt64 bytes_processed = 0) const override;

    /// d'tor
    ~GUIProgressLoggerImpl() override;

private:
    mutable QProgressDialog* dlg_;
    mutable SignedSize begin_;
    mutable SignedSize end_;
    mutable SignedSize current_;
  };
}

