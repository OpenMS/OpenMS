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

#include <OpenMS/KERNEL/Feature.h>

#include <QDialog>

namespace Ui
{
  class FeatureEditDialogTemplate;
}

namespace OpenMS
{
  /**
      @brief Dialog for editing a feature

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI FeatureEditDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    FeatureEditDialog(QWidget * parent);
    /// Destructor
    ~FeatureEditDialog() override;

    /// Sets the feature
    void setFeature(const Feature & feature);
    /// Returns the feature
    const Feature & getFeature() const;

protected:

    /// The feature to edit
    mutable Feature feature_;

private:
    ///Not implemented
    FeatureEditDialog();

    Ui::FeatureEditDialogTemplate* ui_;

  };

}
