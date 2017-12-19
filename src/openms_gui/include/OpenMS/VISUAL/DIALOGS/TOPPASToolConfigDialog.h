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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_TOPPASTOOLCONFIGDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASTOOLCONFIGDIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

class QComboBox;
class QPushButton;
class QRadioButton;
class QString;

#include <QtGui/QDialog>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/LayerData.h>

namespace OpenMS
{
  class ParamEditor;

  /**
      @brief TOPP tool configuration dialog

      In the dialog, the user can set the parameters for the tool

      This information can then be used to execute the tool.

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASToolConfigDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /**
        @brief Constructor

        @param parent Qt parent widget
        @param param The param we are editing
        @param default_dir The default directory for loading and storing INI files
        @param tool_name The name of the tool
        @param tool_type The type of the tool
@param hidden_entries List of entries that are used already in edges etc and should not be shown
    */
    TOPPASToolConfigDialog(QWidget * parent, Param & param, const String& default_dir, const String& tool_name, const String& tool_type, const String& tool_desc, const QVector<String>& hidden_entries);
    ///Destructor
    ~TOPPASToolConfigDialog() override;

private:
    /// ParamEditor for reading ini-files
    ParamEditor * editor_;
    /// The param we are editing
    Param * param_;
    /// Param for loading the ini-file
    Param arg_param_;
    /// default-dir of ini-file to open
    String default_dir_;
    /// name of ini-file
    QString filename_;
    /// The name of the tool
    String tool_name_;
    /// The type of the tool
    String tool_type_;
    /// The parameters already explained by in edges
    QVector<String> hidden_entries_;

protected slots:
    /// Slot for OK button
    void ok_();
    /// loads an ini-file into the editor_
    void loadINI_();
    /// stores an ini-file from the editor_
    void storeINI_();
  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
