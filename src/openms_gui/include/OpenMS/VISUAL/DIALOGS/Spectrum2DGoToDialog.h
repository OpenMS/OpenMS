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

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum2DGoToDialog.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

  class String;

  /**
      @brief GoTo dialog used to zoom to a m/z and retention time range or to a feature.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI Spectrum2DGoToDialog :
    public QDialog,
    public Ui::Spectrum2DGoToDialogTemplate
  {
    Q_OBJECT

public:
    ///Constructor
    Spectrum2DGoToDialog(QWidget * parent = nullptr);
    ///Destructor
    ~Spectrum2DGoToDialog() override;

    /// Returns if a range should be display (true) or if a feature should be displayed (false)
    bool showRange() const;

    /// Fixes the currently stored range (i.e. ensure correct order of min-max; enforce minimum of 1 Da (or 1 sec for RT) window IFF min==max
    void fixRange();

    ///@name Methods for ranges
    //@{
    ///Sets the data range to display initially
    void setRange(float min_rt, float max_rt, float min_mz, float max_mz);
    ///Sets the data range of the complete experiment for better navigation with the dialog
    void setMinMaxOfRange(float min_rt, float max_rt, float min_mz, float max_mz);
    ///Returns the lower RT bound
    float getMinRT() const;
    ///Returns the upper RT bound
    float getMaxRT() const;
    ///Returns the lower m/z bound
    float getMinMZ() const;
    ///Returns the upper m/z bound
    float getMaxMZ() const;
    //@}

    ///@name Methods for feature numbers
    //@{
    ///Returns the selected feature numbers. If a number is returned, the feature rather than the range should be displayed.
    String getFeatureNumber() const;
    ///Disables the feature number field
    void enableFeatureNumber(bool);
    //@}

  };

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H
