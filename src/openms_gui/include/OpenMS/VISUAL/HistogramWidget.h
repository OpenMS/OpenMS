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

// QT
#include <QtWidgets>
#include <QPixmap>
class QPaintEvent;
class QResizeEvent;
class QMouseEvent;

//OpenMS
#include <OpenMS/MATH/STATISTICS/Histogram.h>

namespace OpenMS
{
  class AxisWidget;
  class String;

  /**
      @brief Widget which can visualize a histogram.

      @image html HistogramWidget.png

      It can also be used to define a left and right boundary inside the values.
  It supports normal and log scaling via the context menu.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI HistogramWidget :
    public QWidget
  {
    Q_OBJECT

public:
    /// Constructor
    HistogramWidget(const Math::Histogram<> & distribution, QWidget * parent = nullptr);

    /// Destructor
    ~HistogramWidget() override;

    /// Returns the value f the lower splitter
    double getLeftSplitter() const;

    /// Returns the value of the upper splitter
    double getRightSplitter() const;

    /// Set axis legends
    void setLegend(const String & legend);

public slots:
    /// Shows the splitters if @p on is true. Hides them otherwise.
    void showSplitters(bool on);

    /// Sets the value of the right splitter
    void setRightSplitter(double pos);

    /// Sets the value of the left splitter
    void setLeftSplitter(double pos);

    /// Enables/disables log mode
    void setLogMode(bool log_mode);

protected:
    /// The histogram to display
    Math::Histogram<> dist_;

    /// Flag that indicates if splitters are shown
    bool show_splitters_;

    /// Value of the right splitter
    double left_splitter_;

    /// Value of the right splitter
    double right_splitter_;

    /// The splitter that is currently dragged (0=none, 1=left, 2=right)
    UInt moving_splitter_;

    /// X axis
    AxisWidget * bottom_axis_;

    /// Margin around plot
    UInt margin_;

    /// Internal buffer for the double buffering
    QPixmap buffer_;

    /// Flag that indicates the current mode
    bool log_mode_;

    /// Repaints the contents to the buffer and calls update()
    void invalidate_();

    ///@name reimplemented Qt events
    //@{
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void resizeEvent(QResizeEvent *) override;
    //@}

protected slots:

    /// Context menu event
    void showContextMenu(const QPoint & pos);
  };
} // namespace OpenMS


