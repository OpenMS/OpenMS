// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

//OpenGL
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

//STL
#include <map>

//QT
#include <QtGui/QMainWindow>
#include <QtGui/QWorkspace>
#include <QtCore/QStringList>

class QAction;
class QComboBox;
class QLabel;
class QListWidget;
class QListWidgetItem;
class QDockWidget;
class QToolButton;
class QCloseEvent;
class QTextEdit;

namespace OpenMS
{
  class MultiGradientSelector;
  class Spectrum1DWidget;
  class Spectrum2DWidget;
  class Spectrum3DWidget;

  /**
  	@brief Main window of TOPPView tool
  	
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
      	@param maximize If the window the new file was added to should be displayed maximized
      	@param use_mower If a mower should be used to suppress noise in the data
      	@param force_type File type to force
      */
      void addSpectrum(const String& filename, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER, FileHandler::Type force_type=FileHandler::UNKNOWN);
      /**
      	@brief Opens and displays a spectrum form the database
      	
      	@param db_id The id in the database
      	@param as_new_window If the data is displayed in the current window or in a new window
      	@param maps_as_2d If maps are displayed 2D or 3D
      	@param maximize If the window the new file was added to should be displayed maximized
      	@param use_mower If a mower should be used to suppress noise in the data
      */
      void addDBSpectrum(UInt db_id, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);

      /// maximizes the size of the active window
      void maximizeActiveSpectrum();
      /// opens all the files that are inside the handed over iterator range
      template <class StringListIterator>
      void loadFiles(const StringListIterator& begin, const StringListIterator& end)
      {
        //use mower?
        OpenDialog::Mower mow = OpenDialog::NO_MOWER;
        if ( (String)param_.getValue("preferences:intensity_cutoff")=="noise_estimator")
        {
          mow = OpenDialog::NOISE_ESTIMATOR;
        }
				
				bool last_was_plus = false;
        for (StringListIterator it=begin; it!=end; ++it)
        {
          if (*it=="+")
          {
          	last_was_plus = true;
          	continue;
        	}
        	else if (!last_was_plus)
        	{
        		addSpectrum(*it,true,(String)param_.getValue("preferences:default_map_view")=="2d",true,mow);
        	}
        	else
        	{
        		last_was_plus = false;
        		addSpectrum(*it,false,(String)param_.getValue("preferences:default_map_view")=="2d",true,mow);
        	}
        }
        maximizeActiveSpectrum();
      }

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
      /// shows the dialog for opening spectra from file or the database
      void openSpectrumDialog();
      /// shows the goto dialog
      void gotoDialog();
      /// shows the preferences dialog
      void preferencesDialog();
      /// lets the user edit the preferences of a layer
      void layerPreferencesDialog();
      /// Lets the user change the intensity distribution of a layer
      void layerIntensityDistribution();
			/// Changes the axis visibility of the current window
			void changeAxisVisibility();
      /// lets the user edit the meta data of a layer
      void editMetadata();
      /// saves the contents of the active window
      void saveImage();
      /// saves the content of active window to an image
      void print();
      /// closes the active window
      void closeFile();
      /// saves the current view of the current layer
      void saveLayer();
      /// renames a layer
      void renameLayer();
      /// updates the toolbar
      void updateToolbar();
      /// adapts the layer bar to the active window
      void updateLayerbar();
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
      	
      	If time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
      	Otherwise the message is displayed for @p time ms.
      */
      void showStatusMessage(std::string msg, OpenMS::UInt time);
      /// shows m/z, intensity and rt in the status bar
      void showCursorStatus(double mz, double intensity, double rt);
      /// TOPP tool dialog
      void showTOPPDialog();
      /// Annotates current layer with ID data
      void annotateWithID();
      /// INI File edit dialog
      void editParamDialog();
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
    	/// slot for layer manager context menu
    	void layerContextMenu(const QPoint& pos);
    	/// signal for layer manager visibility change (check box)
    	void layerVisibilityChange(QListWidgetItem* item);
      //@}
      
      /** @name Tabbar slots
      */
      //@{    	
    	/// Closes the window corresponding to the data of the tab with index @p index
      void closeByTab(int index);
      /// Raises the window corresponding to the data of the tab with index @p index
      void focusByTab(int index);
      /// Removes the tab with data @p id
      void removeTab(int id);
      /// Opens a file from the recent files menu
      void openRecentFile();
      //@}
      
      /** @name Toolbar slots
      */
      //@{
      void setActionMode(int);
      void setDrawMode1D(int);
      void setIntensityMode(int);
      void changeGridLines();
      void changeLayerFlag(bool);
      void resetZoom();
      void showProjections();
      //@}
		
		/// Appends process output to log window
		void updateProcessLog();
		
    protected:
      /// Adds a tab for the window in the tabbar
      void addTab_(SpectrumWidget*, const String&);
      /**
      	@brief Opens a SpectrumWidget as a new window
      	
      	Connect the slots/signals for status messages and mode changes (paint or mouse mode).
      	Adds a tab and the window caption..
      */
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

      /// Log output window
      QTextEdit* log_;

      /** @name Toolbar
      */
      //@{
      QToolBar* tool_bar_;
      //common actions
      QButtonGroup* action_group_;
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
      /// If this QAction is checked, the recent files are opened as a new layer
      QAction* recent_as_new_layer_;

			/// list of the recently opened files actions (menu entries)
			std::vector<QAction*> recent_actions_;
			//@}

      /// check if all avaiable preferences get set by the .ini file. If there are some missing entries fill them with default values.
      void checkPreferences_();
      
      void closeEvent(QCloseEvent* event); 
  }
  ; //class

} //namespace

#endif // OPENMS_APPLICATIONS_TOPPVIEWBASE_H
