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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H
#define OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

//QT
#include <QtGui/QMainWindow>
#include <QtGui/QWorkspace>
#include <QtGui/QButtonGroup>
#include <QtCore/QProcess>
#include <QtGui/QSplashScreen>

class QToolBar;
class QListWidget;
class QTextEdit;
class QWorkspace;
class QLabel;
class QWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QWebView;
class QNetworkAccessManager;
class QNetworkReply;

namespace OpenMS
{
	class TOPPASWidget;
	class TOPPASScene;
	class TOPPASTabBar;
	class TOPPASResources;

  /**
  	@brief Main window of the IDEvaluation tool
  
    @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI IDEvaluationBase 
  	: public QMainWindow, 
  		public DefaultParamHandler
  {
      Q_OBJECT

    public:
    	
      ///Constructor
      IDEvaluationBase(QWidget* parent=0);
      ///Destructor
      virtual ~IDEvaluationBase();

      QSize sizeHint() const;

      void setVisibleArea(double low, double high);

      static StringList getSupportedImageFormats();

    public slots:
      
      void resetZoom();

      void setIntensityMode(int index);

			/// compute q-values from ids and store as vector of points for plotting
      bool getPoints(std::vector<PeptideIdentification>& peptides /* cannot be const, to avoid copy */, const std::vector< DoubleReal >& q_value_thresholds, MSSpectrum<>& points);

      /// opens the file in a new window
      void addSearchFile(const String& file_name);
    	/// shows the dialog for opening files
      void openFileDialog();
      void saveCurrentPipelineAs();
      /// saves the pipeline (determined by Qt's sender mechanism)
      void savePipeline();
      /// exports the current pipeline as image
      void exportAsImage(const QString& file_name, const QString& format = "");
    	/// changes the current path according to the currently active window/layer
      //void updateCurrentPath();
      /// Shows the 'About' dialog
      void showAboutDialog();
      /**
      	@brief Shows a status message in the status bar.
      	
      	If @p time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UInt time);
      /// shows x,y coordinates in the status bar
      //void showCursorStatus(double x, double y);
      /// closes the active window
      //void closeFile();
      /// updates the toolbar
      //void updateToolBar();

      void loadFiles(const StringList& list);

      void showURL();

    protected slots:
	  
		  /// enable/disable menu entries depending on the current state
    	//void updateMenu();
    	/// Shows the widget as window in the workspace (the special_id is only used for the first untitled widget (to be able to auto-close it later)
    	//void showAsWindow_(TOPPASWidget* sw, const String& caption, const int special_id = -1);
			



    protected:

			/// Log output window
      QTextEdit* log_;
      /// Workflow Description window
      QTextEdit* desc_;

      /// Main workspace
      QWorkspace* ws_;

      Spectrum1DWidget* spec_1d_;

      /// Label for messages in the status bar
      QLabel* message_label_;

      ///returns the window with id @p id
      TOPPASWidget* window_(int id) const;


      /// The current path (used for loading and storing).
      /// Depending on the preferences this is static or changes with the current window/layer.
      String current_path_;
      ///@name reimplemented Qt events
      //@{
      void closeEvent(QCloseEvent* event);
			void keyPressEvent(QKeyEvent* e);
			//@}
			
			///Log message states
			enum LogState
			{ 
				LS_NOTICE,   ///< Notice
				LS_WARNING,  ///< Warning
				LS_ERROR     ///< Fatal error
			};
			/// Shows a log message in the log_ window 
      void showLogMessage_(LogState state, const String& heading, const String& body);

      std::vector< DoubleReal > q_value_thresholds_;
      /** @name Toolbar
      */
      //@{
      QToolBar* tool_bar_;
      //common intensity modes

      QButtonGroup* intensity_button_group_;
      //1D specific stuff
      //@}

  }
  ; //class

} //namespace

#endif // OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H
