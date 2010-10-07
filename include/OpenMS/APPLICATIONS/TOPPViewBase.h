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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOPPVIEWBASE_H
#define OPENMS_APPLICATIONS_TOPPVIEWBASE_H

//OpenMS
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/SYSTEM/FileWatcher.h>

//STL
#include <map>

//QT
#include <QtGui/QMainWindow>
#include <QtGui/QButtonGroup>
#include <QtGui/QActionGroup>
#include <QtCore/QStringList>
#include <QtCore/QProcess>

class QAction;
class QComboBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QListWidgetItem;
class QTreeWidget;
class QTreeWidgetItem;
class QDockWidget;
class QToolButton;
class QCloseEvent;
class QTextEdit;
class QCheckBox;
class QSplashScreen;
class QWorkspace;
class QToolButton;

namespace OpenMS
{
  class MultiGradientSelector;
  class Spectrum1DWidget;
  class Spectrum2DWidget;
  class Spectrum3DWidget;
  class ToolsDialog;
  class DBConnection;
  class EnhancedTabBar;
  class FileWatcher;

  /**
    @brief Main window of TOPPView tool

    @improvement Use DataRepository singleton to share data between TOPPView and the canvas classes (Hiwi)

		@improvement For painting single mass traces with no width we currently paint each line twice (once going down, and then coming back up).
		This could be more efficient...

    @improvement Keep spectrum browser widgets of all layers in memory in order to avoid rebuilding the entire tree view every time the active layer changes (Hiwi, Johannes)

  	@todo Add TOPPView live-tutorial (Stephan, Marc)

    @ingroup TOPPView_elements
  */
  class OPENMS_DLLAPI TOPPViewBase
  	: public QMainWindow,
  		public DefaultParamHandler
  {
      Q_OBJECT
      
    friend class TestTOPPView;
      
    public:
    	///@name Type definitions
    	//@{
    	//Feature map type
    	typedef LayerData::FeatureMapType FeatureMapType;
      //Feature map managed type
      typedef LayerData::FeatureMapSharedPtrType FeatureMapSharedPtrType;

    	//Consensus feature map type
    	typedef LayerData::ConsensusMapType ConsensusMapType;
      //Consensus  map managed type
      typedef LayerData::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

    	//Peak map type
    	typedef LayerData::ExperimentType ExperimentType;
      //Main managed data type (experiment)
      typedef LayerData::ExperimentSharedPtrType ExperimentSharedPtrType;
    	///Peak spectrum type
    	typedef ExperimentType::SpectrumType SpectrumType;
    	//@}

      ///Constructor
      TOPPViewBase(QWidget* parent=0);
      ///Destructor
      ~TOPPViewBase();

      /**
      	@brief Opens and displays data from a file

      	Loads the data and adds it to the application by calling addData_()

      	@param filename The file to open
      	@param show_options If the options dialog should be shown (otherwise the defaults are used)
      	@param caption Sets the layer name and window caption of the data. If unset the file name is used.
      	@param add_to_recent If the file should be added to the recent files after opening
      	@param window_id in which window the file is opened if opened as a new layer (0 or default equals current window).
      	@param spectrum_id determines the spectrum to show in 1D view.
      */
      void addDataFile(const String& filename, bool show_options, bool add_to_recent, String caption="", UInt window_id=0, Size spectrum_id=0);
      /**
      	@brief Opens and displays a data from a database

      	Loads the data and adds it to the application by calling addData_()

      	@param db_id The id in the database
      	@param show_options If the options dialog should be shown (otherwise the defaults are used)
      	@param caption Sets the layer name and window caption of the data. If unset the file name is used.
      	@param window_id in which window the file is opened if opened as a new layer (0 or default equals current window).
      */
      void addDataDB(UInt db_id, bool show_options, String caption="", UInt window_id=0);

      /// opens all the files that are inside the handed over string list
      void loadFiles(const StringList& list, QSplashScreen* splash_screen);

      /**
      	@brief Loads the preferences from the filename given.

      	If the filename is empty, the application name + ".ini" is used as filename
      */
      void loadPreferences(String filename="");
      /// stores the preferences (used when this window is closed)
      void savePreferences();

			/// Returns the active Layer data (0 if no layer is active)
			const LayerData* getCurrentLayer() const;

    public slots:
      /// changes the current path according to the currently active window/layer
      void updateCurrentPath();
      /// shows the URL stored in the data of the sender QAction
      void showURL();
      /// shows the file dialog for opening files
      void openFileDialog();
      /// shows the file dialog for opening example files
      void openExampleDialog();
      /// shows the DB dialog for opening files
      void openDatabaseDialog();
      /// shows the goto dialog
      void gotoDialog();
      /// shows the preferences dialog
      void preferencesDialog();
      /// Shows statistics (count,min,max,avg) about Intensity, Quality, Charge and meta data
      void layerStatistics();
      /// lets the user edit the meta data of a layer
      void editMetadata();
      /// chooses searched spectrum
      void chooseSpectrumByUser(const QString& text);
      /// closes the active window
      void closeFile();
      /// updates the toolbar
      void updateToolBar();
      /*
      /// updates entries of loaded datas
      void updateDataBar();
      */
      /// adapts the layer bar to the active window
      void updateLayerBar();
      /// adapts the spectrum bar to the active window
      void updateSpectrumBar();
      /// adapts the filter bar to the active window
      void updateFilterBar();
      /// brings the tab corresponding to the active window in front
      void updateTabBar(QWidget* w);
      /// tile the open windows vertically
      void tileVertical();
      /// tile the open windows horizontally
      void tileHorizontal();
      /**
      	@brief Shows a status message in the status bar.

      	If @p time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UInt time);
      /// shows m/z and rt in the status bar
      void showCursorStatus(double mz, double rt);
      /// shows m/z and rt in the status bar (inverting RT and m/z)
      void showCursorStatusInvert(double mz, double rt);
      /// Apply TOPP tool
      void showTOPPDialog();
      /// Annotates current layer with ID data
      void annotateWithID();
      /// Shows the theoretical spectrum generation dialog
      void showSpectrumGenerationDialog();
      /// Shows the spectrum alignment dialog
      void showSpectrumAlignmentDialog();
      /// Shows the spectrum with index @p index of the active layer in 1D
      void showSpectrumAs1D(int index);
      /// Shows the current peak data of the active layer in 2D
      void showCurrentPeaksAs2D();
      /// Shows the current peak data of the active layer in 3D
      void showCurrentPeaksAs3D();
      /// Shows the 'About' dialog
      void showAboutDialog();
      /// Saves the whole current layer data
      void saveLayerAll();
      /// Saves the visible layer data
      void saveLayerVisible();
			/// Toggles the grid lines
      void toggleGridLines();
    	/// Toggles the axis legends
      void toggleAxisLegends();
			/// Shows current layer preferences
      void showPreferences();
			/// dialog for inspecting database meta data
			void metadataDatabaseDialog();
			/// dialog for inspecting file meta data
			void metadataFileDialog();
			/// Shows the selected spectrum
			void spectrumSelectionChange(QTreeWidgetItem* current, QTreeWidgetItem* previous);
			/// Opens a new 1D window and shows the spectrum, if not already in 1D
			void spectrumDoubleClicked(QTreeWidgetItem* current, int /*col*/);

    protected slots:
      /** @name Layer manager and filter manager slots
      */
      //@{
    	/// slot for layer manager selection change
    	void layerSelectionChange(int);
    	/// Enables/disables the data filters for the current layer
    	void layerFilterVisibilityChange(bool);
    	/// slot for layer manager context menu
    	void layerContextMenu(const QPoint& pos);
    	/// slot for spectrum manager context menu
    	void spectrumContextMenu(const QPoint& pos);
    	/// slot for spectrum browser column header context menu
    	void spectrumBrowserHeaderContextMenu(const QPoint& pos);
    	/// slot for log window context menu
    	void logContextMenu(const QPoint& pos);
    	/// slot for layer manager visibility change (check box)
    	void layerVisibilityChange(QListWidgetItem* item);
    	/// slot for filter manager context menu
    	void filterContextMenu(const QPoint& pos);
    	/// slot for editing a filter
    	void filterEdit(QListWidgetItem* item);
    	/// slot for editing the preferences of the current layer
    	void layerEdit(QListWidgetItem* /*item*/);
    	/// slot for the finished signal of the TOPP tools execution
    	void finishTOPPToolExecution(int exitCode, QProcess::ExitStatus exitStatus);
    	/// aborts the execution of a TOPP tool
    	void abortTOPPTool();
    	/// retuns the last invoked TOPP tool with the same parameters
			void rerunTOPPTool();
    	/// enabled/disabled menu entries depending on the current state
    	void updateMenu();
    	/// shows the spectrum browser and updates it
    	void showSpectrumBrowser();
      //@}

      /** @name Tabbar slots
      */
      //@{
    	/// Closes the window corresponding to the data of the tab with identifier @p id
      void closeByTab(int id);
      /// Raises the window corresponding to the data of the tab with identifier @p id
      void focusByTab(int id);
      /// Opens a file from the recent files menu
      void openRecentFile();
			/// Slot for drag-and-drop of layer manager to tabbar
			void copyLayer(const QMimeData* data, QWidget* source, int id=-1);
      //@}

      /** @name Toolbar slots
      */
      //@{
      void setDrawMode1D(int);
      void setIntensityMode(int);
      void changeLayerFlag(bool);
      void changeLabel(QAction*);
			void changeUnassigned(QAction*);
      void resetZoom();
      void toggleProjections();
      //@}

			/// Appends process output to log window
			void updateProcessLog();

      /// Shows the tutorial browser
      void showTutorial();      

      /// Called if a data file has been externally changed
      void fileChanged_(const String&);
    protected:
  		/**
  			@brief Adds a peak or feature map to the viewer

  			@param feature_map The feature data (empty if not feature data)
  			@param consensus_map The consensus feature data (empty if not consensus feature data)
				@param peptides The peptide identifications (empty if not ID data)
  			@param peak_map The peak data (empty if not peak data)
				@param data_type Type of the data
  			@param show_as_1d Force dataset to be opened in 1D mode (even if it contains several spectra)
  			@param show_options If the options dialog should be shown (otherwise the defaults are used)
  			@param filename source file name (if the data came from a file)
      	@param caption Sets the layer name and window caption of the data. If unset the file name is used. If set, the file is not monitored foro changes.
      	@param window_id in which window the file is opened if opened as a new layer (0 or default equals current
      	@param spectrum_id determines the spectrum to show in 1D view.
      */
      void addData_(FeatureMapSharedPtrType feature_map, ConsensusMapSharedPtrType consensus_map, std::vector<PeptideIdentification>& peptides, ExperimentSharedPtrType peak_map, LayerData::DataType data_type, bool show_as_1d, bool show_options, const String& filename="", const String& caption="", UInt window_id=0, Size spectrum_id=0);

      /// unique list of files referenced by all layers
       std::set<String> getFilenamesOfOpenFiles();

    	/// Tries to open a db connection (queries the user for the DB password)
    	void connectToDB_(DBConnection& db);
    	/**
    		@brief Shows a dialog where the user can select files
    	*/
    	QStringList getFileList_(const String& path_overwrite = "");

      ///Returns the parameters for a SpectrumCanvas of dimension @p dim
      Param getSpectrumParameters_(UInt dim);

      void showAsWindow_(SpectrumWidget* sw, const String& caption);
      ///returns the window with id @p id
      SpectrumWidget* window_(int id) const;
      ///returns a pointer to the active SpectrumWidget (0 if none is active)
      SpectrumWidget*  activeWindow_() const;
      ///returns a pointer to the active SpectrumCanvas (0 if none is active)
      SpectrumCanvas*  activeCanvas_() const;
      ///returns a pointer to the active Spectrum1DWidget (0 the active window is no Spectrum1DWidget or there is no active window)
      Spectrum1DWidget* active1DWindow_() const;
      ///returns a pointer to the active Spectrum2DWidget (0 the active window is no Spectrum2DWidget or there is no active window)
      Spectrum2DWidget* active2DWindow_() const;
      ///returns a pointer to the active Spectrum3DWidget (0 the active window is no Spectrum2DWidget or there is no active window)
      Spectrum3DWidget* active3DWindow_() const;

      /// Layer management widget
      QListWidget* layer_manager_;

      /*
      /// data management widget (shows loaded files)
      QTreeWidget* data_manager_view_;
      */
      /// Watcher that tracks file changes (in order to update the data in the different views)
      FileWatcher* watcher_;

      /// Holds the messageboxes for each layer that are currently popped up (to avoid popping them up again, if file changes again before the messagebox is closed)
      bool watcher_msgbox_;

      ///@name Spectrum selection widgets
      //@{
      QTreeWidget* spectrum_selection_;
      QDockWidget* spectrum_bar_;
      QLineEdit* spectrum_search_box_;
      QComboBox* spectrum_combo_box_;
      //@}

      ///@name Data filter widgets
      //@{
      QListWidget* filters_;
      QCheckBox* filters_check_box_;
      //@}

      /// Log output window
      QTextEdit* log_;

      /** @name Toolbar
      */
      //@{
      QToolBar* tool_bar_;
      //common intensity modes

      QButtonGroup* intensity_group_;
      //1D specific stuff

      QToolBar* tool_bar_1d_;
      QButtonGroup* draw_group_1d_;

      //2D specific stuff
      QToolBar* tool_bar_2d_peak_;
      QToolBar* tool_bar_2d_feat_;
      QToolBar* tool_bar_2d_cons_;
			QToolBar* tool_bar_2d_ident_;
      QAction* dm_precursors_2d_;
      QAction* dm_hull_2d_;
      QAction* dm_hulls_2d_;
      QToolButton* dm_label_2d_;
			QActionGroup* group_label_2d_;
      QToolButton* dm_unassigned_2d_;
			QActionGroup* group_unassigned_2d_;
      QAction* dm_elements_2d_;
      QAction* projections_2d_;
			QAction* dm_ident_2d_;
      //@}

      /// Main workspace
      QWorkspace* ws_;

      ///Tab bar. The address of the corresponding window to a tab is stored as an int in tabData()
      EnhancedTabBar* tab_bar_;

      /** @name Status bar
      */
      //@{
      /// Label for messages in the status bar
      QLabel* message_label_;
      /// m/z label for messages in the status bar
      QLabel* mz_label_;
      /// RT label for messages in the status bar
      QLabel* rt_label_;
			//@}

      /// @name Recent files
      //@{
      ///adds a Filename to the recent files
      void addRecentFile_(const String& filename);
      ///update the recent files menu
      void updateRecentMenu_();
      /// list of the recently opened files
      QStringList recent_files_;
			/// list of the recently opened files actions (menu entries)
			std::vector<QAction*> recent_actions_;
			//@}


      /// @name TOPP tool execution
      //@{
			/// Runs the TOPP tool according to the information in topp_
			void runTOPPTool_();
			///Information needed for execution of TOPP tools
			struct
			{
				Param param;
				String tool;
				String in;
				String out;
				String file_name;
				String layer_name;
				UInt window_id;
				Size spectrum_id;
				QProcess* process;
				bool visible;
			} topp_;
			//@}

      /// check if all avaiable preferences get set by the .ini file. If there are some missing entries fill them with default values.
      void checkPreferences_();
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

      ///Additional context menu for 2D layers
      QMenu* add_2d_context_;

  		/// Apply TOPP tool. If @p visible is true, only the visible data is used, otherwise the whole layer is used.
      void showTOPPDialog_(bool visible);

      /// The current path (used for loading and storing).
      /// Depending on the preferences this is static or changes with the current window/layer.
      String current_path_;

      // static helper functions
      public:        
        /// Returns true if @p contains at least one MS1 spectrum
        static bool containsMS1Scans(const ExperimentType& exp);
        
        /// Estimates the noise by evaluating n_scans random scans of MS level 1. Assumes that 4/5 of intensities is noise.
        float estimateNoiseFromRandomMS1Scans(const ExperimentType& exp, UInt n_scans = 10);

        /// Counts the number of exact zero valued intensities in all MS1 spectra
        UInt countZeros(const ExperimentType& exp);
  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPVIEWBASE_H
