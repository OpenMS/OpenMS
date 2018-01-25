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

#ifndef OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPViewOpenDialog.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

class QAbstractButton;

namespace OpenMS
{
  class Param;
  class String;
  /**
      @brief Dataset opening options for TOPPView

      @ingroup TOPPView_elements
  */
  class OPENMS_GUI_DLLAPI TOPPViewOpenDialog :
    public QDialog,
    public Ui::TOPPViewOpenDialogTemplate
  {
    Q_OBJECT

public:
    /// Constructor
    TOPPViewOpenDialog(const String & data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent = nullptr);
    /// Destructor
    ~TOPPViewOpenDialog() override;

    /// Returns true, if 2D mode is to be used for maps
    bool viewMapAs2D() const;
    /// Returns true, if 1D mode is to be used for maps
    bool viewMapAs1D() const;
    /// Returns of the low intensity peaks should be hidden
    bool isCutoffEnabled() const;
    /// Returns true, if the data should be opened in a new window
    bool openAsNewWindow() const;
    ///Returns the index of the selected merge layer. If the option is not selected -1 is returned.
    Int getMergeLayer() const;

    /// Disables view dimension section and sets the selected option
    void disableDimension(bool as_2d);
    /// Disables cutoff section and sets the selected option
    void disableCutoff(bool cutoff_on);
    /// Disables opening location section and sets the selected option
    void disableLocation(bool window);
    /**
        @brief Sets the possible merge layers (index and name) and activates the the option

        It is deactivated by default and can be deactivated manually by passing an empty list.
    */
    void setMergeLayers(const Map<Size, String> & layers);

protected slots:
    ///slot that disables 2D/3D options, when as layer is selected
    void updateViewMode_(QAbstractButton * button);

protected:
    ///Stores if this option is disabled, to avoid activating it in updateViewMode_()
    bool map_as_2d_disabled_;
  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H
