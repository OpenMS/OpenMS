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
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

class QString; // declare this OUTSIDE of namespace OpenMS!
class QStringList;

namespace OpenMS
{

  /// Macro for Qt's connect() overload resolution (in case signals/slots are overloaded and we need to tell connect what overload to pick
  /// without repeating ourselves.
  /// This can be solved in Qt 5.7 by using qOverload<>
  /// @note: provide the brackets for 'args' yourself, since there might be multiple arguments, separated by comma
  /// Example: QObject::connect(spinBox, CONNECTCAST(QSpinBox, valueChanged, (double)), slider, &QSlider::setValue);
  #define CONNECTCAST(class,func,args) static_cast<void(class::*)args>(&class::func)

  /**
    @brief Class which holds static GUI-related helper functions.

    Since all methods are static, the c'tor is private.
    
    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI GUIHelpers
  {
  public:

    /// Open a folder in file explorer
    /// Will show a message box on failure
    static void openFolder(const QString& folder);

    /// Open TOPPView (e.g. from within TOPPAS)
    static void startTOPPView(const QStringList& args);

    /// Open a certain URL (in a browser)
    /// Will show a message box on failure
    static void openURL(const QString& target);

  private:
    /// private C'tor
    GUIHelpers();
  };

}
