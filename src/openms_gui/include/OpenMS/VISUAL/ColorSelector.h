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

//QT
#include <QtWidgets/QWidget>

class QPaintEvent;
class QMouseEvent;

namespace OpenMS
{

  /**
      @brief A widget for selecting a color.

      It represents a color (displayed as background color) and allows changing the color.

      \image html ColorSelector.png

      The above example image shows four ColorSelector instances on the right side.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI ColorSelector :
    public QWidget
  {
    Q_OBJECT

public:
    /// Constructor
    ColorSelector(QWidget * parent = nullptr);

    /// Destructor
    ~ColorSelector() override;

    /// Returns the selected color
    const QColor & getColor();

    /// Sets the selected color
    void setColor(const QColor &);

    /// Qt size hint
    QSize sizeHint() const override;
protected:
    ///@name Reimplemented Qt events
    //@{
    void paintEvent(QPaintEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    //@}
    QColor color_;
  };

}
