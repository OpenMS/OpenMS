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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtGui/QDialog>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum1DGoToDialog.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

  /**
      @brief simple goto/set visible area dialog for exact placement of the viewing window

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI Spectrum1DGoToDialog :
    public QDialog,
    public Ui::Spectrum1DGoToDialogTemplate
  {
    Q_OBJECT

public:
    ///Constructor
    Spectrum1DGoToDialog(QWidget * parent = nullptr);
    ///Destructor
    ~Spectrum1DGoToDialog() override;

    ///Sets the m/z range displayed initially
    void setRange(float min, float max);

    ///Sets the m/z range displayed initially
    void setMinMaxOfRange(float min, float max);

    /// Fixes the currently stored range (i.e. ensure correct order of min-max; enforce minimum of 1 Da window IFF min==max
    void fixRange();

    ///Returns the lower m/z bound
    float getMin() const;
    ///Returns the upper m/z bound
    float getMax() const;
  };

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H
