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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASIOMAPPINGDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASIOMAPPINGDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <QtCore/QVector>


namespace OpenMS
{
  class TOPPASEdge;

  /**
      @brief Dialog which allows to configure the input/output parameter mapping of an edge.

      This dialog allows to select an output parameter of the source vertex and an input
      parameter of the target vertex. Only valid selections are allowed, i.e. the type
      (file or list of files) and at least one valid file type of either vertex must match.

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASIOMappingDialog :
    public QDialog,
    public Ui::TOPPASIOMappingDialogTemplate
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASIOMappingDialog(TOPPASEdge * parent);

public slots:

    /// Called instead of exec() after edge is constructed (in order to avoid showing the dialog if not necessary)
    int firstExec();

protected:

    /// Fills the table
    void fillComboBoxes_();

    /// The edge we are configuring
    TOPPASEdge * edge_;

    /// Vector storing the mapping of the target input combobox indices to param indices of edges
    QVector<int> target_input_param_indices_;

protected slots:

    /// Called when OK is pressed; checks if the selected parameters are valid
    void checkValidity_();

  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
