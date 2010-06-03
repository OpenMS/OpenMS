// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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


#ifndef OPENMS_VISUAL_TOPPASWIDGET_H
#define OPENMS_VISUAL_TOPPASWIDGET_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QGraphicsView>

namespace OpenMS
{
	class TOPPASScene;
	class Param;

  /**
  	@brief Widget visualizing and allowing to edit TOPP pipelines.
  	
  	This class is a subclass of QGraphicsView and visualizes a TOPPASScene.
  	Several TOPPASWidgets can be opened in TOPPAS at the same time,
  	managed by a QWorkspace.
		
		@ingroup TOPPAS_elements
  */
  class OPENMS_DLLAPI TOPPASWidget
  	: public QGraphicsView
  {
      Q_OBJECT

    public:
    
      /// Default constructor
      TOPPASWidget(const Param& preferences, QWidget* parent = 0, const String& tmp_path = "");

      /// Destructor
      virtual ~TOPPASWidget();
      
      /// Widget id used as identifier
			Int window_id;
			
			/// Returns the scene
			TOPPASScene* getScene();
			/// Zooms in or out, depending on @p zoom_in
			void zoom(bool zoom_in);
		
		signals:
		
			/// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
			void sendStatusMessage(std::string message, OpenMS::UInt time);
			/// Emitted when the cursor position changes (for displaying e.g. in status bar)
			void sendCursorStatus(double x=0.0, double y=0.0);
			/// Message about the destruction of this widget
		  void aboutToBeDestroyed(int w_id);
		  /// Emitted when a drop event occurs
		  void toolDroppedOnWidget(double x = 0.0, double y = 0.0);
		  /// Emitted when a drop event occurs
      void pipelineDroppedOnWidget(const String& filename, bool new_window);
		
		protected:
		
			/// The scene visualized by this widget
			TOPPASScene* scene_;
			
			///@name reimplemented QT events
			//@{
			void wheelEvent(QWheelEvent* event);
			void keyPressEvent(QKeyEvent* e);
			void keyReleaseEvent(QKeyEvent* e);
			void leaveEvent(QEvent* e);
			void enterEvent(QEvent* e);
			void dragEnterEvent(QDragEnterEvent* event);
			void dragMoveEvent(QDragMoveEvent* event);
			void dropEvent(QDropEvent* event);
			void resizeEvent(QResizeEvent* event);
			void closeEvent(QCloseEvent* e);
			//@}
  };
}

#endif
