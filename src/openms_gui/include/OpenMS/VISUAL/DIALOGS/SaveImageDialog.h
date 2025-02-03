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

#include <QtWidgets/QDialog>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QCheckBox>

namespace OpenMS
{
  /**
      @brief Dialog for saving an image.

      @image html SaveImageDialog.png

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI SaveImageDialog :
    public QDialog
  {
    Q_OBJECT

public:
    ///Constructor
    SaveImageDialog(QWidget * parent = nullptr);
    ///set size and size ratio
    void setSize(int x, int y);
    ///accessors for the width
    int getXSize();
    ///accessors for the height
    int getYSize();
    ///accessors for the format
    QString getFormat();

public slots:
    ///changes width keeping proportions
    void xSizeChanged(const QString & s);
    ///changes height keeping proportions
    void ySizeChanged(const QString & s);
    ///set size ratio when proportions checkbox is activated
    void proportionsActivated(bool state);
    ///checks if the values for width and height are ok before accepting the dialog
    void checkSize();

private:
    //format
    QComboBox * format_;
    //size
    QLineEdit * size_x_;
    QLineEdit * size_y_;
    QCheckBox * size_proportions_;
    //ratio size_x_/size_y_
    float size_ratio_;

    //set the size ratio (width/height)
    void setSizeRatio_(float r);
  };
}
