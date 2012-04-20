// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASTREEVIEW_H
#define OPENMS_VISUAL_TOPPASTREEVIEW_H

#include <OpenMS/config.h>

//QT
#include <QtGui/QTreeWidget>
#include <QtGui/QMouseEvent>
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
    TOPPASTreeView(QWidget * parent = 0);
    /// Destructor
    ~TOPPASTreeView();

protected:
    ///@name Reimplemented Qt events
    //@{
    void mousePressEvent(QMouseEvent * e);
    void mouseMoveEvent(QMouseEvent * e);
    void keyPressEvent(QKeyEvent * e);
    void leaveEvent(QEvent * e);
    void enterEvent(QEvent * e);
    //@}

    /// The drag start position
    QPoint drag_start_pos_;
  };

}
#endif // OPENMS_VISUAL_TOPPASTREEVIEW_H
