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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_APPLICATIONS_INIFILEEDITORWINDOW_H
#define OPENMS_VISUAL_APPLICATIONS_INIFILEEDITORWINDOW_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtGui/QMainWindow>

class QToolBar;
class QAction;
class QString;
class QFileDialog;

namespace OpenMS
{
  /**
      @brief shows the ParamEditor widget in a QMainWindow with a toolbar
  */
  class OPENMS_GUI_DLLAPI INIFileEditorWindow :
    public QMainWindow
  {
    Q_OBJECT

public:
    /// menu is created here
    INIFileEditorWindow(QWidget * parent = nullptr);
    /// when user closes window a message box asks the user if he wants to save
    void closeEvent(QCloseEvent * event) override;

public slots:
    ///loads the xml-file into a Param object and loads Param into ParamEditor
    bool openFile(const String & filename = "");
    /// saves the users changes in a xml-file if the Param object is valid
    bool saveFile();
    /// like saveFile but with a file dialog to choose a filename
    bool saveFileAs();
    /// if the user changes data in ParamEditor the title shows a '*'
    void updateWindowTitle(bool);

private:
    /// ParamEditor object for visualization
    ParamEditor * editor_;
    /// Param object for storing data
    Param param_;
    /// filename of xml-file to store the Param object
    QString filename_;
    /// path used as next default location of the load/store dialogs
    String current_path_;
  };
}

#endif
