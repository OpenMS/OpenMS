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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_DATAFILTERDIALOG_H
#define OPENMS_VISUAL_DIALOGS_DATAFILTERDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_DataFilterDialog.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>


namespace OpenMS
{
  /**
      @brief Dialog for creating and changing a DataFilter

  */
  class OPENMS_GUI_DLLAPI DataFilterDialog :
    public QDialog,
    public Ui::DataFilterDialogTemplate
  {
    Q_OBJECT

public:
    /// constructor
    DataFilterDialog(DataFilters::DataFilter & filter, QWidget * parent);

protected slots:
    /// Checks if the settings are valid and writes them to filter_ if so
    void check_();
    /// Is called when field_ changes and enables/disables the meta data functionality as needed
    void field_changed_(const QString &);
    /// Is called when op_ changes and disables the value field, if operation is "exists", else enables it
    void op_changed_(const QString &);

protected:
    /// Reference to the filter that is modified
    DataFilters::DataFilter & filter_;

private:
    ///Not implemented
    DataFilterDialog();
  };

}
#endif // OPENMS_VISUAL_DIALOGS_OPENDIALOG_H
