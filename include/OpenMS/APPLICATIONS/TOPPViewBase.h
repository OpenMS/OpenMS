// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOPPVIEWBASE_H
#define OPENMS_APPLICATIONS_TOPPVIEWBASE_H

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/FORMAT/FileHandler.h>

//STL
#include <map>

//QT
#include <QtGui/QMainWindow>
#include <QtGui/QWorkspace>
#include <QtGui/QButtonGroup>
#include <QtCore/QStringList>
#include <QtCore/QProcess>

class QAction;
class QComboBox;
class QLabel;
class QListWidget;
class QListWidgetItem;
class QDockWidget;
class QToolButton;
class QCloseEvent;
class QTextEdit;
class QCheckBox;

namespace OpenMS
{
  class MultiGradientSelector;
  class Spectrum1DWidget;
  class Spectrum2DWidget;
  class Spectrum3DWidget;
  class ToolsDialog;

  /**
  	@brief Main window of TOPPView tool
		
		@todo Fix mysterious crash when the second file is not readable from command line (Marc)
		@todo Add layer context menu to main menu; Make canvas context menu extensible (Marc)
		@todo Rerun TOPP tool - add option to apply it on the visible data only (Marc)
  	@todo Projections: fix painting outside of widget boundaries, repaint when the user does not zoom/translate for X seconds, add splitter to resize (Marc)
		@todo Speed up 2D view: paint only highest point per pixel (Marc)
  	
  	@ingroup TOPPView_elements
  */
  class TOPPViewBase 
  	: public QMainWindow, 
  		public DefaultParamHandler
  {
      Q_OBJECT

    public:
      ///Constructor
      TOPPViewBase(QWidget* parent=0);
      ///Destructor
      ~TOPPViewBase();
      
      /**
      	@brief Opens and displays a spectrum form a file
      	
      	@param filename The file to open
      	@param as_new_window If the data is displayed in the current window or in a new window
      	@param maps_as_2d If maps are displayed 2D or 3D
      	@param use_mower If a mower should be used to suppress noise in the data
      	@param force_type File type to force
      	@param caption Sets the layer name and window caption of the data. If unset the file name is used.
      	@param window_id in which window the file is opend of opened as a new layer (0 or default equals current window).
      */
      void addSpectrum(const String& filename, bool as_new_window=true, bool maps_as_2d=true, bool use_mower=false, FileHandler::Type force_type=FileHandler::UNKNOWN, String caption="", UInt window_id=0);
      /**
      	@brief Opens and displays a spectrum form the database
      	
      	@param db_id The id in the database
      	@param as_new_window If the data is displayed in the current window or in a new window
      	@param maps_as_2d If maps are displayed 2D or 3D
      	@param use_mower If a mower should be used to suppress noise in the data
      */
      void addDBSpectrum(UInt db_id, bool as_new_window=true, bool maps_as_2d=true, bool use_mower=false);

      /// opens all the files that are inside the handed over string list
      void loadFiles(const StringList& list);

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
      /// shows the URL stored in the data of the sender QAction
      void showURL();
      /// shows the dialog for opening files
      void openFileDialog();
      /// shows the dialog for opening files
      void openDatabaseDialog();
      /// shows the goto dialog
      void gotoDialog();
      /// shows the preferences dialog
      void preferencesDialog();
      /// Shows statistics (count,min,max,avg) about Intensity, Quality, Charge and meta data
      void layerStatistics();
      /// lets the user edit the meta data of a layer
      void editMetadata();
      /// closes the active window
      void closeFile();
      /// updates the toolbar
      void updateToolBar();
      /// adapts the layer bar to the active window
      void updateLayerBar();
      /// adapts the filter bar to the active window
      void updateFilterBar();
      /// brings the tab corresponding to the active window in front
      void updateTabBar(QWidget* w);
      /// tile the open windows vertically
      void tileVertical();
      /// tile the open windows horizontally
      void tileHorizontal();
      /// Links/unlinks two spectra (for zooming)
      void linkActiveTo(int);
      /**
      	@brief Shows a status message in the status bar.
      	
      	If @p time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UInt time);
      /// shows m/z, intensity and rt in the status bar
      void showCursorStatus(double mz, double intensity, double rt);
      /// TOPP tool dialog
      void showTOPPDialog();
      /// Annotates current layer with ID data
      void annotateWithID();
      /// Shows the current peak data of the active layer in 3D 
      void showCurrentPeaksAs3D();
			/// Shows the spectrum with index @p index of the sctive layer in 1D
			void showSpectrumAs1D(int index);
      /// Shows the 'About' dialog
      void showAboutDialog();

    protected slots:
      /** @name Layer manager slots
      */
      //@{  
    	/// slot for layer manager selection change
    	void layerSelectionChange(int);
    	/// Enables/disables the data filters for the current layer
    	void layerFilterVisibilityChange(bool);
    	/// slot for layer manager context menu
    	void layerContextMenu(const QPoint& pos);
    	/// slot for layer manager visibility change (check box)
    	void layerVisibilityChange(QListWidgetItem* item);
    	/// slot for filter manager context menu
    	void filterContextMenu(const QPoint& pos);
    	/// slot for editing a filter
    	void filterEdit(QListWidgetItem* item);
    	/// slot for the finished signal of the TOPP tools execution
    	void finishTOPPToolExecution(int exitCode, QProcess::ExitStatus exitStatus);
    	/// aborts the execution of a TOPP tool
    	void abortTOPPTool();
    	/// enabled/disabled menu entries depending on the current state
    	void updateMenu();
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
      //@}
      
      /** @name Toolbar slots
      */
      //@{
      void setDrawMode1D(int);
      void setIntensityMode(int);
      void changeLayerFlag(bool);
      void resetZoom();
      void showProjections();
      //@}
		
			/// Appends process output to log window
			void updateProcessLog();
		
    protected:      
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
      ///Estimates the noise by evaluating 10 random scans of MS level 1
      float estimateNoise_(const SpectrumCanvas::ExperimentType& exp);

      /// Layer mangment widget
      QListWidget* layer_manager_;

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
      QComboBox* link_box_;
      //2D specific stuff
      QToolBar* tool_bar_2d_;
      QAction* dm_precursors_2d_;
      QAction* dm_hull_2d_;
      QAction* dm_hulls_2d_;
      QAction* dm_numbers_2d_;
      QAction* projections_2d_;
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
      /// Intensity label for messages in the status bar
      QLabel* int_label_;
      /// RT label for messages in the status bar
      QLabel* rt_label_;
			//@}
			
      /**
        @brief Map that stores linked pairs of 1D windows.
      	
      	Each link is stored twice (both directions).
      */
      std::map<int,int> link_map_;

      /** @name Recent files
      */
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
			
			///@name Members for execution of TOPP tools
			//@{
			QProcess* process_;
			ToolsDialog* tools_dialog_;
			String topp_filename_;
			String topp_layer_name_;
			UInt topp_window_id_;
			//@}

      /// check if all avaiable preferences get set by the .ini file. If there are some missing entries fill them with default values.
      void checkPreferences_();
      ///reimplemented Qt close event
      void closeEvent(QCloseEvent* event); 
  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPVIEWBASE_H
