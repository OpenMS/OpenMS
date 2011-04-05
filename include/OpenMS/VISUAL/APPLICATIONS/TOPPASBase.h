// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/TOPPASTreeView.h>

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
	class TOPPASScene;
	class TOPPASTabBar;
	class TOPPASResources;
	
  /**
  	@brief Main window of the TOPPAS tool
  
    @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASBase 
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

    public slots:    	
			/// opens the file in a new window
      void addTOPPASFile(const String& file_name, bool in_new_window = true);
    	/// shows the dialog for opening files
      void openFileDialog();
			/// shows the dialog for opening example files
			void openExampleDialog();
      /// shows the dialog for creating a new file
      void newPipeline();
      /// shows the dialog for including another workflow in the currently opened one
      void includePipeline();
			/// shows the dialog for saving the current file and updates the current tab caption
      void saveCurrentPipelineAs();
      /// saves the pipeline (determined by qt's sender mechanism)
      void savePipeline();
      /// shows a file dialog for selecting the resource file to load
      void loadPipelineResourceFile();
      /// shows a file dialog for selecting the resource file to write to
      void savePipelineResourceFile();
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
			/// Called when a file was successfully written to an output vertex
			void outputVertexFinished(const String& file);
			/// Called when a TOPP tool produces (error) output.
			void updateTOPPOutputLog(const QString& out);
			/// Called by the scene if the pipeline execution finishes successfully
      void showPipelineFinishedLogMessage();
			/// Saves @p scene to the clipboard
			void saveToClipboard(TOPPASScene* scene);
			/// Sends the clipboard content to the sender of the connected signal
			void sendClipboardContent();
      /// Refreshes the parameters of the TOPP tools of the current workflow and stores an updated workflow including the current parameters
      void refreshParameters();
      /// Open files in a new TOPPView instance
      void openFilesInTOPPView(QStringList all_files);
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
			
			/// Offset counter for new inserted nodes (to avoid invisible stacking)
			static int node_offset_;
      
			/// z-value counter for new inserted nodes (new nodes should be on top)
			static qreal z_value_;

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

			/// The clipboard
      TOPPASScene* clipboard_scene_;

    public:
      /// @name common functions unsed in TOPPAS and TOPPView
      //@{
      /// Creates and fills a tree widget with all available tools
      static TOPPASTreeView* createTOPPToolsTreeWidget(QWidget* parent_widget = 0);

      /// Saves the workflow in the provided TOPPASWidget to a user defined location.
      /// Returns the full file name or "" if no valid one is selected.
      static QString savePipelineAs(TOPPASWidget* w, QString current_path);

      /// Loads and sets the resources of the TOPPASWidget.
      static QString loadPipelineResourceFile(TOPPASWidget* w, QString current_path);

      /// Saves the resources of the TOPPASWidget.
      static QString savePipelineResourceFile(TOPPASWidget* w, QString current_path);

      /// Refreshes the TOPP tools parameters of the pipeline
      static QString refreshPipelineParameters(TOPPASWidget* tw, QString current_path);      
      //@}
  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPASBASE_H
