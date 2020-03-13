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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class PythonModuleRequirement;
}

namespace OpenMS
{
  namespace Internal
  {
    ///Preferences dialog for Spectrum3DWidget
    class OPENMS_GUI_DLLAPI PythonModuleRequirement : public QWidget
    {
      Q_OBJECT
      
    public:
      explicit PythonModuleRequirement(QWidget* parent = nullptr);
      ~PythonModuleRequirement();

      void setRequiredModules(const StringList& m);
      void setFreeText(const String& text);


    signals:
      /// emitted whenever the requirement check was executed...
      void valueChanged(StringList& valid_modules, StringList& missing_modules);


    private slots:
      void validate_(const String& python_exe);

    private:
      StringList required_modules_; ///< list of modules which are needed (order might be important -- know your Python...)
      String info_text_; ///< additional text to display for the user

      Ui::PythonModuleRequirement* ui_;
    };

  }
} // ns OpenMS

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a PythonModuleRequirement member
using PythonModuleRequirement = OpenMS::Internal::PythonModuleRequirement;
