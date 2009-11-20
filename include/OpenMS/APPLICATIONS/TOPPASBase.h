// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOPPASBASE_H
#define OPENMS_APPLICATIONS_TOPPASBASE_H

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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

namespace OpenMS
{
	class TOPPASWidget;
	class TOPPASTabBar;
	
  /**
  	@brief Main window of the TOPPAS tool
  
    @ingroup TOPPAS_elements
  */
  class OPENMS_DLLAPI TOPPASBase 
  	: public QMainWindow, 
  		public DefaultParamHandler
  {
      Q_OBJECT

    public:
    	
      ///Constructor
      TOPPASBase(QWidget* parent=0);
      ///Destructor
      virtual ~TOPPASBase();
			
			/**
      	@brief Loads the preferences from the filename given.
      	
      	If the filename is empty, the application name + ".ini" is used as filename
      */
      void loadPreferences(String filename="");
      /// stores the preferences (used when this window is closed)
      void savePreferences();
			/// loads the files and updates the splashscreen
			void loadFiles(const StringList& list, QSplashScreen* splash_screen);
			/// opens the file in a new window
			void openFile(const String& file_name);

    public slots:
    	
    	/// shows the dialog for opening files
      void openFileDialog();
			/// shows the dialog for opening example files
			void openExampleDialog();
      /// shows the dialog for creating a new file
      void newFileDialog();
      /// shows the dialog for saving the current file
      void saveFileDialog();
			/// shows the dialog for saving the current file as new file
      void saveAsFileDialog(TOPPASWidget* tw = 0);
      /// shows the preferences dialog
      void preferencesDialog();
    	/// changes the current path according to the currently active window/layer
      void updateCurrentPath();
    	/// brings the tab corresponding to the active window in front
      void updateTabBar(QWidget* w);
      /// refreshes the definitions of the TOPP tools
      void refreshDefinitions();
      /// Shows the 'About' dialog
      void showAboutDialog();
      /// shows the URL stored in the data of the sender QAction
      void showURL();
      /**
      	@brief Shows a status message in the status bar.
      	
      	If @p time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UInt time);
      /// shows x,y coordinates in the status bar
      void showCursorStatus(double x, double y);
      /// closes the active window
      void closeFile();
      /// updates the toolbar
      void updateToolBar();
    	/// Runs the pipeline of the current window
			void runPipeline();
			/// Terminates the current pipeline
			void abortPipeline();
			/// Called when a tool is started
			void toolStarted();
			/// Called when a tool is finished
			void toolFinished();
			/// Called when a tool crashes
			void toolCrashed();
			/// Called when a tool execution fails
			void toolFailed();
			/// Called when a pipeline execution ended successfully in an output vertex
			void outputVertexFinished(const String& file);
			/// Called when a TOPP tool produces (error) output.
			void updateTOPPOutputLog(const QString& out);
			/// Called by the scene if the pipeline execution finishes successfully
			void showSuccessLogMessage();
			
    protected slots:
		
			/** @name Tabbar slots
      */
      //@{
    	/// Closes the window corresponding to the data of the tab with identifier @p id
      void closeByTab(int id);
      /// Raises the window corresponding to the data of the tab with identifier @p id
      void focusByTab(int id);
		  //@}
		  
		  /// enable/disable menu entries depending on the current state
    	void updateMenu();
    	/// Shows the widget as window in the workspace
    	void showAsWindow_(TOPPASWidget* sw, const String& caption);
			/// Inserts a new TOPP tool in the current window at (x,y)
			void insertNewVertex_(double x, double y, QTreeWidgetItem* item = 0);
			/// Inserts the @p item in the middle of the current window
			void insertNewVertexInCenter_(QTreeWidgetItem* item);
			
    protected:

			/// Log output window
      QTextEdit* log_;

      /** @name Toolbar
      */
      //@{
      QToolBar* tool_bar_;
      //@}

      /// Main workspace
      QWorkspace* ws_;

      ///Tab bar. The address of the corresponding window to a tab is stored as an int in tabData()
      TOPPASTabBar* tab_bar_;
      
      /// Tree view of all available TOPP tools
      QTreeWidget* tools_tree_view_;
      /// List of ready analysis pipelines
      QListWidget* blocks_list_;
      
      /** @name Status bar
      */
      //@{
      /// Label for messages in the status bar
      QLabel* message_label_;
			//@}
      
      ///returns the window with id @p id
      TOPPASWidget* window_(int id) const;
      
      
      /// The current path (used for loading and storing).
      /// Depending on the preferences this is static or changes with the current window/layer.
      String current_path_;
			
			/// The path for temporary files
			String tmp_path_;
      
      ///returns a pointer to the active TOPPASWidget (0 if none is active)
      TOPPASWidget* activeWindow_() const;
      
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
  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPASBASE_H
