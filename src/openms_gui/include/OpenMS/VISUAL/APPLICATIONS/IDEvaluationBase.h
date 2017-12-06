// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H
#define OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

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
  // OpenMS forward declarations
  class TOPPASWidget;
  class TOPPASScene;
  class TOPPASTabBar;
  class TOPPASResources;

  class Spectrum1DWidget;

  /**
    @brief Main window of the IDEvaluation tool

    @htmlinclude OpenMS_IDEvaluationBase.parameters

    @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI IDEvaluationBase :
    public QMainWindow,
    public DefaultParamHandler
  {
    Q_OBJECT

public:

    ///Constructor
    IDEvaluationBase(QWidget * parent = nullptr);
    ///Destructor
    ~IDEvaluationBase() override;

    QSize sizeHint() const override;

    void setVisibleArea(double low, double high);

    const PeakMap& getPoints() const;

    static StringList getSupportedImageFormats();

public slots:

    void resetZoom();

    void setIntensityMode(int index);

    /// compute q-values from ids and store as vector of points for plotting
    /// returns false on error, the return vector 'points' will also be empty in this case
    bool getPoints(std::vector<PeptideIdentification> & peptides /* cannot be const, to avoid copy */, const std::vector<double> & q_value_thresholds, MSSpectrum & points);

    /// calls 'getPoints()' after loading the idXML file and returns the result
    bool loadCurve(const String& file_name, MSSpectrum& points);
    /// opens the file in a new window
    /// @return false on error (no idXML file or missing information preventing FDR computation)
    bool addSearchFile(const String & file_name);
    /// shows the dialog for opening files
    void openFileDialog();
    /// saves the plot - querying for a filename first
    void saveImageAs();
    /// exports the current pipeline as image
    /// returns true on success, otherwise false + @p error_message
    bool exportAsImage(const QString & file_name, String & error_message, const QString & format = "");
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

    /// load Target/Decoy annotated files, return FALSE if any of these files did
    /// not contain target/decoy information or any other error which prevents FDR calculation
    bool loadFiles(const StringList & list);

    void showURL();

protected slots:

    /// enable/disable menu entries depending on the current state
    //void updateMenu();
    /// Shows the widget as window in the workspace (the special_id is only used for the first untitled widget (to be able to auto-close it later)
    //void showAsWindow_(TOPPASWidget* sw, const String& caption, const int special_id = -1);




protected:

    /// Log output window
    QTextEdit * log_;
    /// Workflow Description window
    QTextEdit * desc_;

    /// Main workspace
    QWorkspace * ws_;

    Spectrum1DWidget * spec_1d_;

    /// Label for messages in the status bar
    QLabel * message_label_;

    ///returns the window with id @p id
    TOPPASWidget * window_(int id) const;


    /// The current path (used for loading and storing).
    /// Depending on the preferences this is static or changes with the current window/layer.
    String current_path_;
    ///@name reimplemented Qt events
    //@{
    void closeEvent(QCloseEvent * event) override;
    void keyPressEvent(QKeyEvent * e) override;
    //@}

    ///Log message states
    enum LogState
    {
      LS_NOTICE,             ///< Notice
      LS_WARNING,            ///< Warning
      LS_ERROR               ///< Fatal error
    };
    /// Shows a log message in the log_ window
    void showLogMessage_(LogState state, const String & heading, const String & body);

    std::vector<double> q_value_thresholds_;

    // holds the computed curves for easy export to outside
    PeakMap data_;

    /** @name Toolbar
    */
    //@{
    QToolBar * tool_bar_;
    //common intensity modes

    QButtonGroup * intensity_button_group_;
    //1D specific stuff
    //@}

  }; //class

} //namespace

#endif // OPENMS_VISUAL_APPLICATIONS_IDEVALUATIONBASE_H
