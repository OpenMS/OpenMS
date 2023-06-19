// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/PlotCanvas.h> // for AreaXYType

#include <QtWidgets/QDialog>

namespace Ui
{
  class Plot2DGoToDialogTemplate;
}

namespace OpenMS
{

  class String;

  /**
      @brief GoTo dialog used to zoom to a m/z and retention time range or to a feature.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI Plot2DGoToDialog :
    public QDialog
  {
    Q_OBJECT

public:
    using AreaXYType = PlotCanvas::AreaXYType;

    ///Constructor
    /// @param parent Parent widget
    /// @param x_name Name of the x_axis dimension
    /// @param y_name Name of the y_axis dimension
    Plot2DGoToDialog(QWidget* parent, std::string_view x_name, std::string_view y_name);
    ///Destructor
    ~Plot2DGoToDialog() override;

    /// Returns if a feature UID was set an a feature should be displayed (false), otherwise, show a range (true)
    bool showRange() const;

    bool checked();

    ///@name Methods for ranges
    //@{
    ///Sets the data range to display initially
    void setRange(const AreaXYType& range);
    ///Sets the data range of the complete experiment for better navigation with the dialog
    void setMinMaxOfRange(const AreaXYType& max_range);

    /// Query the range set by the user.
    /// If any dimension is <1, it is extended to at least 1 to ensure proper displaying.
    AreaXYType getRange();
    //@}

    ///@name Methods for feature numbers
    //@{
    ///Returns the selected feature numbers. If a number is returned, the feature rather than the range should be displayed.
    String getFeatureNumber() const;
    ///Disables the feature number field
    void enableFeatureNumber(bool);
    //@}

  private:
    Ui::Plot2DGoToDialogTemplate* ui_;

  };

}
