// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: SpectrumMDIWindow.h,v 1.53 2006/06/09 22:32:48 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUMMDIWINDOW_H
#define OPENMS_VISUAL_SPECTRUMMDIWINDOW_H

//OpenGL
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>

//STL
#include <map>

//QT
#include <qmainwindow.h>
#include <qworkspace.h>
class QAction;
class QComboBox;
class QToolButton;
class QMenuBar;
class QLabel;
class QRadioButton;
class QActionGroup;
class QPopupMenu;

namespace OpenMS
{
	class MultiGradientSelector;
	class LayerManager;
	class Spectrum1DWindow;
	class Spectrum2DWindow;
	class Spectrum3DWindow;
	class SpectrumWindow;
	
	/**
		@brief MDI window for several SpectrumWindow instances
		
		@todo Add PeakPicker dialog (Eva)
		@todo Fix Feature Convex hull display (Marc)
		@todo Fix multi-dataset Display color problem (Marc)
		@todo How to display MS-MS/MS Scans (Marc)
		
		@ingroup spectrum_widgets
	*/	
	class SpectrumMDIWindow : public QMainWindow, public PreferencesManager 
	{
		Q_OBJECT
		
		public:
			static SpectrumMDIWindow* instance();
			
			/// Accessor for the workspace
			inline QWorkspace* getWorkspace() 
			{ 
				return ws_; 
			}
			
			/**
				@brief Opens and displays a spectrum form a file
				
				@param filename The file to open
				@param as_new_window If the data is displayed in the current window or in a new window
				@param maps_as_2d If maps are displayed 2D or 3D
				@param maximize If the window the new file was added to should be displayed maximized
				@param use_mower If a mower should be used to suppress noise in the data
			*/
			void addSpectrum(const String& filename, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);
			/**
				@brief Opens and displays a spectrum form the database
				
				@param db_id The id in the database
				@param as_new_window If the data is displayed in the current window or in a new window
				@param maps_as_2d If maps are displayed 2D or 3D
				@param maximize If the window the new file was added to should be displayed maximized
				@param use_mower If a mower should be used to suppress noise in the data
			*/
			void addDBSpectrum(UnsignedInt db_id, bool as_new_window=true, bool maps_as_2d=true, bool maximize=false, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);

			void addTab(SpectrumWindow*, const String&);
			/// maximizes the size of the active window
			void maximizeActiveSpectrum();
			/// opens all the files that are inside the handed over iterator range
			template <class StringListIterator>
			void loadFiles(const StringListIterator& begin, const StringListIterator& end)
			{
				for (StringListIterator it=begin; it!=end; ++it)
				{
					addSpectrum(*it,true,bool(getPrefAsInt("Preferences:DefaultMapView2D")),true);
				}
				maximizeActiveSpectrum();
				tab_bar_->setCurrentTab(PointerSizeInt(&(*ws_->activeWindow())));
			}
			/// connect the slots/signals for status messages and mode changes (paint or mouse mode)
			void connectWindowSignals(SpectrumWindow* sw);	
			/// returns selected peaks of the active spectrum framed by \c data_set_.begin() and the last peak BEFORE \c data_set_.end();
			std::vector<MSExperiment<>::SpectrumType::Iterator> getActiveSpectrumSelectedPeaks();
		
			
			///returns a pointer to the active SpectrumWindow (0 if none is sctive)
			SpectrumWindow*  activeWindow() const;
			///returns a pointer to the active Spectrum1DWindow (0 the active window is no Spectrum1DWindow or there is no active window)
			Spectrum1DWindow* active1DWindow() const;
			///returns a pointer to the active Spectrum2DWindow (0 the active window is no Spectrum2DWindow or there is no active window)
			Spectrum2DWindow* active2DWindow() const;
			///returns a pointer to the active Spectrum3DWindow (0 the active window is no Spectrum2DWindow or there is no active window)
			Spectrum3DWindow* active3DWindow() const;
			
			/**
				@brief Loads the preferences from the filename given.
				
				If the filename is empty, the application name + ".ini" is used as filename
			*/
			void loadPreferences(std::string filename="");
			/// stores the preferences (used when this window is closed)
			void savePreferences();

			/// PreferencesManager
			virtual PreferencesDialogPage* createPreferences(QWidget* parent);

		public slots:
			/// shows the dialog for opening spectra from file or the database
			void openSpectrumDialog();
			/// shows the goto dialog
		  void gotoDialog();
		  /// shows the preferences dialog
			void preferencesDialog();
			/// saves the contents of the active window
			void saveImage();
			/// saves the content of active window to an image
			void print(); 
			/// closes the active window
			void closeFile();
			/// updates the toolbar, when the active window changes
			void updateToolbar(QWidget* w);
			/// updates the toolbar, when the active layer of a 1D Window changes
			void update1DToolbar(QWidget* w);
			/// updates the toolbar, when the active layer of a 2D Window changes
			void update2DToolbar(QWidget* w);
		  /// updates the toolbar, when the active layer of a 2D Window changes
			void update3DToolbar(QWidget* w);
			/// adapts the layer bar to the active window
			void updateLayerbar();
			/// brings the tab corresponding to the active window in front
			void updateTabBar(QWidget* w);
			/// tile the open windows vertically
			void tileVertical();
			/// tile the open windows horizontally
			void tileHorizontal();
			//links or unlinks two spectra (for zooming)
			void linkActiveTo(const QString&);
			/**
				@brief Shows a status message in the status bar.
				
				If time is 0 the status message is displayed until showStatusMessage is called with an empty message or a new message.
				Otherwise the message is displayed for @p time ms.
			*/
			void showStatusMessage(std::string msg, OpenMS::UnsignedInt time);
			/// shows m/z, intensity and rt in the status bar 
			void showCursorStatus(double mz, double intensity, double rt);
		
			/// Shows a list all selected peaks
			void showPeaklistActiveSpectrum();
			//picks peaks in the active spectrum
			void pickActiveSpectrum();
			//finds features in the active spectrum
			void findFeaturesActiveSpectrum();
		
		private slots:
			void closeFileByTab(OpenMS::SignedInt);
			void focusSpectrumByAddress(int);
			void removeWidgetFromBar(QObject*);
			void setActionMode(QAction*);
			void setDrawMode(QAction*);
			void setSnapToMax(bool); 	
			void showGridLines(bool); 	
			void switchAxis(bool);
			void setMirroredXAxis(bool);
			void setMirroredYAxis(bool);
			void showPoints(bool);
			void showColors(bool);
			void showContours(bool);
		  void setIntensityScaledDots(bool on);
			void setActionMode2D(QAction*);
		  void setActionMode3D(QAction*);
			void resetZoom();
			void openRecentFile(int i);
		  void setBackView3D(bool on);
		  void setTopView3D(bool on);
		  void setResetZoomView3D(bool on);
		
		protected slots:
			///use this event to do the cleanup
		  virtual void closeEvent(QCloseEvent * e);
			/// Call whenever a window is closed
		  virtual void windowClosed();
		  
		protected:
			///singleton instance
			static SpectrumMDIWindow* instance_;
			///not accessable as this class is a singleton
			SpectrumMDIWindow(QWidget* parent=0, const char* name="SpectrumMDIWindow", WFlags f=0);
			///not accessable as this class is a singleton
			~SpectrumMDIWindow();     
		
			//toolbar
			void createToolBar_();			
			QToolBar* tool_bar_;
			//layers
			QToolBar* layer_bar_;
			LayerManager* layer_manager_;
			
			// 1DWidget
			QActionGroup* action_modes_;
			QActionGroup* draw_modes_;
			QAction* set_zoom_action_;
			QAction* set_translate_action_;
			QAction* set_pick_action_;
			QAction* set_peak_mode_;
			QAction* set_connected_lines_mode_;
			QToolButton* snap_button_;
			QToolButton* grid_button_;
			QToolButton* print_button_;
			QToolButton* switch_axis_button_;
			QToolButton* mirror_xaxis_button_;
			QToolButton* mirror_yaxis_button_;
			QToolButton* reset_zoom_button_;
			//Combobox for linking spectra
			QComboBox* link_box_;

			// 2DWidget
			QToolBar* tool_bar_2d_;
			QActionGroup* action_modes_2d_;
			QAction* set_zoom_action_2d_;
			QAction* set_translate_action_2d_;
			QAction* set_pick_action_2d_;
			QAction* set_measure_action_2d_;
			QToolButton* show_points_button_2d_;
			QToolButton* show_colors_button_2d_;
		  QToolButton* show_contours_button_2d_;
		  QToolButton* intensity_scaled_dots_button_2d_;
	  	QToolButton* reset_zoom_button_2d_;
		
		
		  //3DWidget
		  QToolBar* tool_bar_3d_;
			QActionGroup* action_modes_3d_;
			QAction* set_zoom_action_3d_;
		  QAction* set_pick_action_3d_;
			QToolButton* show_back_view_3d_;
			QToolButton* show_top_view_3d_;
		  QToolButton*  show_reset_view_3d_;
		
		  /// Main workspace
			QWorkspace* ws_;	
			
			///Tab bar
			EnhancedTabBar* tab_bar_;
			///map (maps int(&(*QWidget)) to *QWidget) used for toolbar and tabbar
			std::map<PointerSizeInt,SpectrumWindow*> id_map_;
			
			/// Label for messages in the status bar
			QLabel* message_label_;
			/// m/z label for messages in the status bar
			QLabel* mz_label_;
			/// Intensity label for messages in the status bar
			QLabel* int_label_;
			/// RT label for messages in the status bar
			QLabel* rt_label_;
			
			/**
			  @brief Map that stores linked pairs of 1D windows (uses int value of addresses to identify the widgets).
				
				Each link is stored twice (both directions).
			*/
			std::map<int,int> link_map_;

			///unlinks active spectrum (for zooming)
			void unlinkActive_();

			///adds a Filename to the recent files
			void addRecentFile_(const String& filename);
		
			///update the recent files menu
			void updateRecentMenu_();
		
		  /// check if all avaiable preferences get set by the .ini file. If there are some missing entries fill them with default values.
			void checkPreferences_(); 
			/// list of the recently opened files
			std::vector<String> recent_files_;
		
			/// Pointer to "Tools" menu: so that derived classes can add stuff into it
			QPopupMenu* tools_menu_;
			/// pointer to the recent files menu
			QPopupMenu* recent_menu_;
			
			///add 1D spectrum
			SpectrumWindow* addSpectrum1D_(const String& filename, const String& caption, bool as_new_window=true, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);
			///add 2D spectrum
			SpectrumWindow* addSpectrum2D_(const String& filename, const String& caption, bool as_new_window=true, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);
			///add 3D spectrum
			SpectrumWindow* addSpectrum3D_(const String& filename, const String& caption, bool as_new_window=true, OpenDialog::Mower use_mower=OpenDialog::NO_MOWER);
			
			/// Adds the result of featurefinding to the canvas
			void setFeatureMap_(Spectrum2DCanvas* canvas, Spectrum2DCanvas::ExperimentType& exp, String caption);

	}; //class
	
} //namespace
#endif
