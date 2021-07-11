// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/VISUAL/EnhancedTabBar.h>
#include <OpenMS/VISUAL/EnhancedWorkspace.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/FilterList.h>
#include <OpenMS/VISUAL/RecentFilesMenu.h>
#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/TOPPViewMenu.h>
#include <OpenMS/VISUAL/TVToolDiscovery.h>

//STL
#include <map>

//QT
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QActionGroup>
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
class QCheckBox;
class QSplashScreen;
class QToolButton;

namespace OpenMS
{
  class DataSelectionTabs;
  class FileWatcher;
  class LogWindow;
  class LayerListView;
  class MultiGradientSelector;
  class Plot1DWidget;
  class Plot2DWidget;
  class Plot3DWidget;
  class ToolsDialog;

  /**
    @brief Main window of TOPPView tool

    Uses the default QMainWindow layout (see Qt documentation) with a central
    widget in the middle (consistent of a EnhancedTabBar and an
    EnhancedWorkspace) and multiple docked widgets around it (to the right and
    below) and multiple tool bars. On top and bottom are a menu bar and a
    status bar.

    The main layout is using 
    - Central Widget: 
      - EnhancedTabBar: tab_bar_
      - EnhancedWorkspace: ws_
    - Docked to the right:
      - layer_dock_widget_
      - views_dockwidget_
      - filter_dock_widget_
    - Docked to the bottom:
      - log_bar (only connected through slots)

    The views_dockwidget_ internally holds a tab widget views_tabwidget_ which
    holds the two different views on the data (spectra and identification view)
    which are implemented using idview_behaviour_ and
    spectraview_behavior_.

    @improvement Use DataRepository singleton to share data between TOPPView and the canvas classes (Hiwi)

    @improvement For painting single mass traces with no width we currently paint each line twice (once going down, and then coming back up).
    This could be more efficient...

    @improvement Keep spectrum browser widgets of all layers in memory in order to avoid rebuilding the entire tree view every time the active layer changes (Hiwi, Johannes)

    @ingroup TOPPView_elements
  */
  class OPENMS_GUI_DLLAPI TOPPViewBase :
    public QMainWindow,
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
    //Main on-disc managed data type (experiment)
    typedef LayerData::ODExperimentSharedPtrType ODExperimentSharedPtrType;
    ///Peak spectrum type
    typedef ExperimentType::SpectrumType SpectrumType;
    //@}

    /// Used for deciding whether new tool/util params should be generated or reused from TOPPView's ini file
    enum class TOOL_SCAN
    {
      /**
         TVToolDiscovery does not generate params for each tool/util unless they are absolutely needed and could not be
         extracted from TOPPView's ini file. This may be useful for testing.
      */
      SKIP_SCAN,
      /// Only generate params for each tool/util if TOPPView's last ini file has an older version. (Default behaviour)
      SCAN_IF_NEWER_VERSION,
      /// Forces TVToolDiscovery to generate params and using them instead of the params in TOPPView's ini file
      FORCE_SCAN
    };

    ///Constructor
    explicit TOPPViewBase(TOOL_SCAN scan_mode = TOOL_SCAN::SCAN_IF_NEWER_VERSION, QWidget* parent = nullptr);
    ///Destructor
    ~TOPPViewBase() override;

    enum class LOAD_RESULT
    {
      OK,
      FILE_NOT_FOUND,       ///< file did not exist
      FILETYPE_UNKNOWN,     ///< file exists, but type could no be determined                                                
      FILETYPE_UNSUPPORTED, ///< filetype is known, but the format not supported as layer data
      LOAD_ERROR            ///< an error occurred while loading the file
    };

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
    LOAD_RESULT addDataFile(const String& filename, bool show_options, bool add_to_recent, String caption = "", UInt window_id = 0, Size spectrum_id = 0);

    /**
      @brief Adds a peak or feature map to the viewer

      @param feature_map The feature data (empty if not feature data)
      @param consensus_map The consensus feature data (empty if not consensus feature data)
      @param peptides The peptide identifications (empty if not ID data)
      @param peak_map The peak data (empty if not peak data)
      @param on_disc_peak_map The peak data managed on disc (empty if not peak data)
      @param data_type Type of the data
      @param show_as_1d Force dataset to be opened in 1D mode (even if it contains several spectra)
      @param show_options If the options dialog should be shown (otherwise the defaults are used)
      @param filename source file name (if the data came from a file)
      @param caption Sets the layer name and window caption of the data. If unset the file name is used. If set, the file is not monitored for changes.
      @param window_id in which window the file is opened if opened as a new layer (0 or default equals current
      @param spectrum_id determines the spectrum to show in 1D view.
    */
    void addData(FeatureMapSharedPtrType feature_map,
                 ConsensusMapSharedPtrType consensus_map,
                 std::vector<PeptideIdentification>& peptides,
                 ExperimentSharedPtrType peak_map,
                 ODExperimentSharedPtrType on_disc_peak_map,
                 LayerData::DataType data_type,
                 bool show_as_1d,
                 bool show_options,
                 bool as_new_window = true,
                 const String& filename = "",
                 const String& caption = "",
                 UInt window_id = 0,
                 Size spectrum_id = 0);

    /// Opens all the files in the string list
    void loadFiles(const StringList& list, QSplashScreen* splash_screen);

    /**
      @brief Loads the preferences from the filename given.

      If the filename is empty, the application name + ".ini" is used as filename
    */
    void loadPreferences(String filename = "");

    /// Stores the preferences (used when this window is closed)
    void savePreferences();

    /// Returns the parameters for a PlotCanvas of dimension @p dim
    Param getSpectrumParameters(UInt dim);

    /// Returns the active Layer data (0 if no layer is active)
    const LayerData* getCurrentLayer() const;

    /// Returns the active Layer data (0 if no layer is active)
    LayerData* getCurrentLayer();

    //@name Accessors for the main gui components.
    //@brief The top level enhanced workspace and the EnhancedTabWidgets resing in the EnhancedTabBar.
    //@{
    /// returns a pointer to the EnhancedWorkspace containing PlotWidgets
    EnhancedWorkspace* getWorkspace();

    /// returns a pointer to the active PlotWidget (0 if none is active)
    PlotWidget* getActivePlotWidget() const;

    /// returns a pointer to the active Plot1DWidget (0 the active window is no Plot1DWidget or there is no active window)
    Plot1DWidget* getActive1DWidget() const;

    /// returns a pointer to the active Plot2DWidget (0 the active window is no Plot2DWidget or there is no active window)
    Plot2DWidget* getActive2DWidget() const;

    /// returns a pointer to the active Plot3DWidget (0 the active window is no Plot2DWidget or there is no active window)
    Plot3DWidget* getActive3DWidget() const;
    //@}

    /// returns a pointer to the active PlotCanvas (0 if none is active)
    PlotCanvas* getActiveCanvas() const;

    /// Opens the provided spectrum widget in a new window
    void showPlotWidgetInWindow(PlotWidget* sw, const String& caption);

public slots:
    /// changes the current path according to the currently active window/layer
    void updateCurrentPath();
    /// shows the file dialog for opening files (a starting directory, e.g. for the example files can be provided; otherwise, uses the current_path_)
    void openFilesByDialog(const String& initial_directory = "");
    /// shows the DB dialog for opening files
    void showGoToDialog();
    /// shows the preferences dialog
    void preferencesDialog();
    /// Shows statistics (count,min,max,avg) about Intensity, Quality, Charge and meta data
    void layerStatistics();
    /// lets the user edit the meta data of a layer
    void editMetadata();
    /// gets called if a layer got activated
    void layerActivated();
    /// gets called when a layer changes in zoom
    void layerZoomChanged();
    /// link the zoom of individual windows
    void linkZoom();
    /// gets called if a layer got deactivated
    void layerDeactivated();
    /// closes the active window
    void closeTab();

    /// returns the last invoked TOPP tool with the same parameters
    void rerunTOPPTool();

    /// calls update*Bar and updateMenu_() to make sure the interface matches the current data
    void updateBarsAndMenus();
    /// updates the toolbar
    void updateToolBar();
    /// adapts the layer bar to the active window
    void updateLayerBar();
    /// adapts view bar to the active window
    void updateViewBar();
    /// activates/deactivates menu entries
    void updateMenu();
    /// adapts the filter bar to the active window
    void updateFilterBar();
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
    /// Annotates current layer with ID data from AccurateMassSearch
    void annotateWithAMS();
    /// Annotates current layer with ID data
    void annotateWithID();
    /// Annotates current chromatogram layer with ID data
    void annotateWithOSW();
    /// Shows the theoretical spectrum generation dialog
    void showSpectrumGenerationDialog();
    /// Shows the spectrum alignment dialog
    void showSpectrumAlignmentDialog();
    /// Shows the current peak data of the active layer in 2D
    void showCurrentPeaksAs2D();
    /// Shows the current peak data of the active layer in 3D
    void showCurrentPeaksAs3D();
    /// Shows the current peak data of the active layer as ion mobility
    void showCurrentPeaksAsIonMobility();
    /// Shows the current peak data of the active layer as DIA data
    void showCurrentPeaksAsDIA();
    /// Saves the whole current layer data
    void saveLayerAll();
    /// Saves the visible layer data
    void saveLayerVisible();
    /// Toggles the grid lines
    void toggleGridLines();
    /// Toggles the axis legends
    void toggleAxisLegends();
    /// Toggles drawing of interesting MZs
    void toggleInterestingMZs();
    /// Shows current layer preferences
    void showPreferences();
    /// dialog for inspecting database meta data
    void metadataFileDialog();

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

    /// list of the recently opened files
    /// called when RecentFileMenu items is clicked
    void openFile(const String& filename);

    /// Enables/disables the data filters for the current layer
    void layerFilterVisibilityChange(bool);

    /// shows a spectrum's metadata with index @p spectrum_index from the currently active canvas
    void showSpectrumMetaData(int spectrum_index);

protected slots:
    /// slot for the finished signal of the TOPP tools execution
    void finishTOPPToolExecution(int exitCode, QProcess::ExitStatus exitStatus);
    /// aborts the execution of a TOPP tool
    void abortTOPPTool();
    /// shows the spectrum browser and updates it
    void showSpectrumBrowser();

    /** @name Tabbar slots
    */
    //@{
    /// Closes the window corresponding to the data of the tab with identifier @p id
    void closeByTab(int id);
    /// Raises the window corresponding to the data of the tab with identifier @p id
    void showWindow(int id);
    /// Slot for drag-and-drop of layer manager to tabbar
    void copyLayer(const QMimeData* data, QWidget* source, int id = -1);
    //@}

    /// Appends process output to log window
    void updateProcessLog();

    /// Called if a data file has been externally changed
    void fileChanged_(const String&);
protected:
    /// Initializes the default parameters on TOPPView construction.
    void initializeDefaultParameters_();

    /**
        @brief Shows a dialog where the user can select files
    */
    QStringList chooseFilesDialog_(const String& path_overwrite = "");

    ///@name dock widgets
    //@{
    QDockWidget* layer_dock_widget_;
    QDockWidget* views_dockwidget_;
    QDockWidget* filter_dock_widget_;
    //@}

    /// Layer management widget
    LayerListView* layers_view_;

    DataSelectionTabs* selection_view_;

    ///@name Filter widget
    //@{
    FilterList* filter_list_;
    //@}

    /// Watcher that tracks file changes (in order to update the data in the different views)
    FileWatcher* watcher_ = nullptr;

    /// Holds the messageboxes for each layer that are currently popped up (to avoid popping them up again, if file changes again before the messagebox is closed)
    bool watcher_msgbox_ = false;

    /// Stores whether the individual windows should zoom together (be linked) or not
    bool zoom_together_ = false;

    /// Log output window
    LogWindow* log_;

    /// Determines TVToolDiscovery scans for tool/utils and generates new params.
    TOOL_SCAN scan_mode_;
    /// Scans for tools/utils and generates a param for each.
    TVToolDiscovery tool_scanner_;

    /** @name Toolbar
    */
    //@{
    QToolBar* tool_bar_;

    // common intensity modes
    QButtonGroup* intensity_button_group_;

    // 1D specific stuff
    QToolBar* tool_bar_1d_;
    QButtonGroup* draw_group_1d_;

    // 2D specific stuff
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
    EnhancedWorkspace ws_;  // not a pointer, but an actual object, so it gets destroyed before the DefaultParamhandler (on which it depends)
    /// Tab bar. The address of the corresponding window to a tab is stored as an int in tabData()
    EnhancedTabBar tab_bar_;
    /// manages recent list of filenames and the menu that goes with it
    RecentFilesMenu recent_files_;  // needs to be declared before 'menu_', because its needed there
    /// manages the menu items (active/inactive) and recent files etc
    TOPPViewMenu menu_;

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
    /// adds a Filename to the recent files
    void addRecentFile_(const String& filename);

    //@}


    /// @name TOPP tool execution
    //@{
    /// Runs the TOPP tool according to the information in topp_
    void runTOPPTool_();
    /// Information needed for execution of TOPP tools
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
      QProcess* process = nullptr;
      QTime timer;
      bool visible_area_only;
    } topp_;
    //@}

    /// check if all available preferences get set by the .ini file. If there are some missing entries fill them with default values.
    void checkPreferences_();
    ///@name reimplemented Qt events
    //@{
    void closeEvent(QCloseEvent* event) override;
    //@}

    ///Additional context menu for 2D layers
    QMenu* add_2d_context_;

    /// Apply TOPP tool. If @p visible is true, only the visible data is used, otherwise the whole layer is used.
    void showTOPPDialog_(bool visible);

    /// The current path (used for loading and storing).
    /// Depending on the preferences this is static or changes with the current window/layer.
    String current_path_;

    /// Adds tool/util params to param_ object by querying them from TVToolDiscovery
    void addToolParamsToIni_();

private:
    /// Suffix appended to caption of tabs when layer is shown in 3D
    static const String CAPTION_3D_SUFFIX_;
  }; //class

} //namespace

