// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
