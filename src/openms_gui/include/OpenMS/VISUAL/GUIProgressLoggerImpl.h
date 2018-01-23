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

#ifndef OPENMS_VISUAL_GUIPROGRESSLOGGERIMPL_H
#define OPENMS_VISUAL_GUIPROGRESSLOGGERIMPL_H

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
    /// create new object (needed by Factory)
    static ProgressLogger::ProgressLoggerImpl* create();

    /// name of the model (needed by Factory)
    static const String getProductName();

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
      @brief Implement ProgressLoggerImpl::endProgress().
    */
    void endProgress(const int /* current_recursion_depth */) const override;

    /// d'tor
    ~GUIProgressLoggerImpl() override;

private:
    mutable QProgressDialog* dlg_;
    mutable SignedSize begin_;
    mutable SignedSize end_;
  };
}

#endif
