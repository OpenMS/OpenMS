// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/MultiGradient.h>

//QT
#include <QtWidgets>

class QPaintEvent;
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

//std lib
#include <vector>

namespace OpenMS
{

  /**
      @brief A widget witch allows constructing gradients of multiple colors.

      @image html MultiGradientSelector.png

      The above example image shows a MultiGradientSelector.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI MultiGradientSelector :
    public QWidget
  {
    Q_OBJECT
public:
    ///Constructor
    MultiGradientSelector(QWidget * parent = nullptr);
    ///Destructor
    ~MultiGradientSelector() override;

    ///returns a const reference to the gradient
    const MultiGradient & gradient() const;
    ///returns a mutable reference to the gradient
    MultiGradient & gradient();

    /// sets the interpolation mode
    void setInterpolationMode(MultiGradient::InterpolationMode mode);
    /// returns the interpolation mode
    MultiGradient::InterpolationMode getInterpolationMode() const;

public slots:
    /// sets what interpolation mode is used
    void stairsInterpolation(bool state);

protected:
    ///@name re-implemented Qt events
    //@{
    void paintEvent(QPaintEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    void mouseMoveEvent(QMouseEvent * e) override;
    void mouseReleaseEvent(QMouseEvent * e) override;
    void mouseDoubleClickEvent(QMouseEvent * e) override;
    void keyPressEvent(QKeyEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    //@}

    // the actual gradient
    MultiGradient gradient_;

    // border margin
    Int margin_;
    // height of the gradient area
    Int gradient_area_width_;
    // height of the lever area
    Int lever_area_height_;

    //position (0-100) in the vector of the selected lever
    Int selected_;
    //color of the selected lever
    QColor selected_color_;

    //stores if the left button is pressed
    bool left_button_pressed_;
  };

}
