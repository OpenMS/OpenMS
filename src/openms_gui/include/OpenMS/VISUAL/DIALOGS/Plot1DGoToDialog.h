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

#include <OpenMS/CONCEPT/Types.h>

#include <QtWidgets/QDialog>

namespace Ui
{
  class Plot1DGoToDialogTemplate;
}

namespace OpenMS
{

  /**
      @brief simple goto/set visible area dialog for exact placement of the viewing window

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI Plot1DGoToDialog :
    public QDialog
  {
    Q_OBJECT

public:
    ///Constructor
    Plot1DGoToDialog(QWidget * parent = nullptr);
    ///Destructor
    ~Plot1DGoToDialog() override;

    ///Sets the m/z range displayed initially
    void setRange(float min, float max);

    ///Sets the m/z range displayed initially
    void setMinMaxOfRange(float min, float max);


    bool checked();

    /// Fixes the currently stored range (i.e. ensure correct order of min-max; enforce minimum of 1 Da window IFF min==max
    void fixRange();

    ///Returns the lower m/z bound
    float getMin() const;
    ///Returns the upper m/z bound
    float getMax() const;

  private:
    Ui::Plot1DGoToDialogTemplate* ui_;

  };

}
