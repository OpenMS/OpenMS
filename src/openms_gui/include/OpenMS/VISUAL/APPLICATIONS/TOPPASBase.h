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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/VISUAL/EnhancedWorkspace.h>
#include <OpenMS/VISUAL/TOPPASTreeView.h>

//QT
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMdiArea>
#include <QtNetwork/QNetworkReply>
#include <QtWidgets/QSplashScreen>

class QToolBar;
class QListWidget;
class QTextEdit;
class QMdiArea;
class QLabel;
class QWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QWebView;
class QNetworkAccessManager;


namespace OpenMS
{
  class EnhancedWorkSpace;
  class EnhancedTabBar;
  class TOPPASWidget;
  class TOPPASScene;
  class LogWindow;
  class TOPPASResources;

  /**
    @brief Main window of the TOPPAS tool

    @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASBase :
    public QMainWindow,
    public DefaultParamHandler
  {
    Q_OBJECT

public:

    ///Constructor
    TOPPASBase(QWidget* parent = nullptr);
    ///Destructor
    ~TOPPASBase() override;

    /**
@brief Loads the preferences from the filename given.

If the filename is empty, the application name + ".ini" is used as filename
*/
    void loadPreferences(String filename = "");
    /// stores the preferences (used when this window is closed)
    void savePreferences();
    /// loads the files and updates the splash screen
    void loadFiles(const StringList& list, QSplashScreen* splash_screen);

public slots:
    /// opens the file in a new window
    void addTOPPASFile(const String& file_name, bool in_new_window = true);
    /// shows the dialog for opening files
    void openFilesByDialog();
    /// shows the dialog for opening example files
    void openExampleDialog();
    /// creates a new tab
    void newPipeline();
    /// shows the dialog for including another workflow in the currently opened one
    void includePipeline();
    /// shows the dialog for saving the current file and updates the current tab caption
    void saveCurrentPipelineAs();
    /// saves the pipeline (determined by Qt's sender mechanism)
    void savePipeline();
    /// exports the current pipeline as image
    void exportAsImage();
    /// shows a file dialog for selecting the resource file to load
    void loadPipelineResourceFile();
    /// shows a file dialog for selecting the resource file to write to
    void savePipelineResourceFile();
    /// opens the OpenMS Homepage to download example workflows
    void openOnlinePipelineRepository();
    /// shows the preferences dialog
    void preferencesDialog();
    /// changes the current path according to the currently active window/layer
    void updateCurrentPath();
    /// brings the tab corresponding to the active window in front
    void updateTabBar(QMdiSubWindow* w);
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
    /// Opens a toppas file
    void openToppasFile(QString filename);
protected slots:

    /** @name Tab bar slots
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
    void insertNewVertex_(double x, double y, QTreeWidgetItem* item = nullptr);
    /// Inserts the @p item in the middle of the current window
    void insertNewVertexInCenter_(QTreeWidgetItem* item);

    /// triggered when user clicks a link - if it ends in .TOPPAS we're done
    void downloadTOPPASfromHomepage_(const QUrl& url);
    /// triggered when download of .toppas file is finished, so we can store & open it
    void toppasFileDownloaded_(QNetworkReply* r);
    /// debug...
    void TOPPASreadyRead();

    /// user edited the workflow description
    void descriptionUpdated_();

protected:

    /// Log output window
    LogWindow* log_;
    /// Workflow Description window
    QTextEdit* desc_;

    /** @name Toolbar
    */
    //@{
    QToolBar* tool_bar_;
    //@}

    /// Main workspace
    EnhancedWorkspace* ws_;

    /// OpenMS homepage workflow browser
    QWebView* webview_;
    /// download .toppas files from homepage
    QNetworkAccessManager* network_manager_;
    /// the content of the network request
    QNetworkReply* network_reply_;

    ///Tab bar. The address of the corresponding window to a tab is stored as an int in tabData()
    EnhancedTabBar* tab_bar_;

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
    TOPPASWidget* activeSubWindow_() const;

    ///@name reimplemented Qt events
    //@{
    void closeEvent(QCloseEvent* event) override;
    void keyPressEvent(QKeyEvent* e) override;
    //@}

    /// The clipboard
    TOPPASScene* clipboard_scene_;


public:
    /// @name common functions used in TOPPAS and TOPPView
    //@{
    /// Creates and fills a tree widget with all available tools
    static TOPPASTreeView* createTOPPToolsTreeWidget(QWidget* parent_widget = nullptr);

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
  }; //class

} //namespace

