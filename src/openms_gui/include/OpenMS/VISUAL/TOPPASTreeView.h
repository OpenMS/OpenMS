// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//QT
#include <QtWidgets/QTreeWidget>
#include <QMouseEvent>
#include <QtCore/QPoint>

namespace OpenMS
{
  class String;

  /**
      @brief Tree view implementation for the list of TOPP tools

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI TOPPASTreeView :
    public QTreeWidget
  {
    Q_OBJECT

public:
    /// Constructor
    TOPPASTreeView(QWidget * parent = nullptr);
    /// Destructor
    ~TOPPASTreeView() override;

    /// Filter tree elements by name (case insensitive and partial=substring matches are valid)
    /// An empty filter shows all elements.
    /// If an element in a subtree is matched, all parents up to the root are also shown.
    void filter(const QString& must_match);

    /// expand all subtrees, i.e. make them visible
    void expandAll();

    /// collapse all subtrees; only show the uppermost level
    void collapseAll();

protected:
    ///@name Reimplemented Qt events
    //@{
    void mousePressEvent(QMouseEvent * e) override;
    void mouseMoveEvent(QMouseEvent * e) override;
    void keyPressEvent(QKeyEvent * e) override;
    void leaveEvent(QEvent * e) override;
    void enterEvent(QEnterEvent* e) override;
    //@}

    /// The drag start position
    QPoint drag_start_pos_;
  };

}
