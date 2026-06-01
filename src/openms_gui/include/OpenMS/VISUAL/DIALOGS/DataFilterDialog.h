// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/PROCESSING/MISC/DataFilters.h>

#include <QDialog>

namespace Ui
{
  class DataFilterDialogTemplate;
}

namespace OpenMS
{
  /**
      @brief Dialog for creating and changing a DataFilter

  */
  class OPENMS_GUI_DLLAPI DataFilterDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// constructor
    DataFilterDialog(DataFilters::DataFilter & filter, QWidget * parent);

    /// destructor
    virtual ~DataFilterDialog();

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
    /// Not implemented
    DataFilterDialog();

    Ui::DataFilterDialogTemplate* ui_;
  };

}
