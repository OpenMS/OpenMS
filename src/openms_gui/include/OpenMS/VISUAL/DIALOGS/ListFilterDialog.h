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

#include <QDialog>

namespace Ui
{
  class ListFilterDialog;
}

namespace OpenMS
{
  /**
      @brief Dialog for creating and changing a DataFilter

  */
  class OPENMS_GUI_DLLAPI ListFilterDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// constructor
    ListFilterDialog() = delete;

    /**
      @brief C'tor with items to show and select from

      @param parent Parent widget
      @param items A set of strings to show and select from. Can be filtered in the dialog
      @param items_prechosen A set of strings which are already chosen (on the right side) when first showing this dialog. This must be a subset of @p items

      @throws Exception::InvalidValue if any of @p items_prechosen is not contained in @p items

    **/     
    ListFilterDialog(QWidget* parent, const QStringList& items = QStringList(), const QStringList& items_prechosen = QStringList());

    /// destructor
    virtual ~ListFilterDialog();

    /// when pressing 'X' button in corner of the Window
    void closeEvent(QCloseEvent* event) override;

    /// A set of strings to show and select from. Can be filtered in the dialog
    /// @throws Exception::InvalidValue if any of @p items_prechosen is not contained in @p items
    void setItems(const QStringList& items);

    /// A set of strings which are already chosen (on the right side). Overwrites the currently chosen set.
    /// @throws Exception::InvalidValue if any of @p items_prechosen is not contained in @p items
    void setPrechosenItems(const QStringList& items_prechosen);

    /// get all items which where selected by the user
    QStringList getChosenItems() const;

protected slots:
    /// button '>>' clicked
    void BtnLRClicked_();
    /// button '> ALL >' clicked
    void BtnLRAllClicked_();
    /// button '<<' clicked
    void BtnRLClicked_();
    /// button '< ALL <' clicked
    void BtnRLAllClicked_();

private:
    Ui::ListFilterDialog* ui_;
  };

}
