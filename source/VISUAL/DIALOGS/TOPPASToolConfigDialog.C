// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QStringList>
#include <QtGui/QPushButton>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QMessageBox>
#include <QtGui/QRadioButton>
#include <QtGui/QFileDialog>
#include <QtGui/QCheckBox>

using namespace std;

namespace OpenMS
{
  TOPPASToolConfigDialog::TOPPASToolConfigDialog(QWidget * parent, Param & param, String default_dir, String tool_name, String tool_type, QVector<String> hidden_entries) :
    QDialog(parent),
    param_(&param),
    default_dir_(default_dir),
    tool_name_(tool_name),
    tool_type_(tool_type),
    hidden_entries_(hidden_entries)
  {
    QGridLayout * main_grid = new QGridLayout(this);

    //Add advanced mode check box
    editor_ = new ParamEditor(this);
    main_grid->addWidget(editor_, 0, 0, 1, 1);

    QHBoxLayout * hbox = new QHBoxLayout;
    QPushButton * load_button = new QPushButton(tr("&Load"));
    connect(load_button, SIGNAL(clicked()), this, SLOT(loadINI_()));
    hbox->addWidget(load_button);
    QPushButton * store_button = new QPushButton(tr("&Store"));
    connect(store_button, SIGNAL(clicked()), this, SLOT(storeINI_()));
    hbox->addWidget(store_button);
    hbox->addStretch();

    // cancel button
    QPushButton * cancel_button = new QPushButton(tr("&Cancel"));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    hbox->addWidget(cancel_button);

    // ok button
    QPushButton * ok_button_ = new QPushButton(tr("&Ok"));
    connect(ok_button_, SIGNAL(clicked()), this, SLOT(ok_()));
    hbox->addWidget(ok_button_);

    main_grid->addLayout(hbox, 1, 0, 1, 1);

    setLayout(main_grid);

    editor_->load(*param_);

    String str;

    editor_->setFocus(Qt::MouseFocusReason);

    setWindowTitle(tool_name.toQString() + " " + tr("configuration"));
  }

  TOPPASToolConfigDialog::~TOPPASToolConfigDialog()
  {

  }

  void TOPPASToolConfigDialog::ok_()
  {
    if (editor_->isModified())
    {
      editor_->store();
      accept();
    }
    else
    {
      reject();
    }
  }

  void TOPPASToolConfigDialog::loadINI_()
  {
    QString string;
    filename_ = QFileDialog::getOpenFileName(this, tr("Open ini file"), default_dir_.c_str(), tr("ini files (*.ini);; all files (*.*)"));
    //not file selected
    if (filename_.isEmpty())
    {
      return;
    }
    if (!arg_param_.empty())
    {
      arg_param_.clear();
      param_->clear();
      editor_->clear();
    }
    try
    {
      arg_param_.load(filename_.toStdString());
    }
    catch (Exception::BaseException & e)
    {
      QMessageBox::critical(this, "Error", (String("Error loading INI file: ") + e.getMessage()).c_str());
      arg_param_.clear();
      return;
    }
    //Extract the required parameters
    *param_ = arg_param_.copy(tool_name_ + ":1:", true);
    //param_->remove("log");
    //param_->remove("no_progress");
    //param_->remove("debug");

    //remove parameters already explained by edges and the "type" parameter
    foreach(const String &name, hidden_entries_)
    {
      param_->remove(name);
    }

    //load data into editor
    editor_->load(*param_);
  }

  void TOPPASToolConfigDialog::storeINI_()
  {
    //nothing to save
    if (param_->empty())
      return;

    filename_ = QFileDialog::getSaveFileName(this, tr("Save ini file"), default_dir_.c_str(), tr("ini files (*.ini)"));
    //not file selected
    if (filename_.isEmpty())
    {
      return;
    }
    if (!filename_.endsWith(".ini"))
      filename_.append(".ini");
    editor_->store();
    arg_param_.insert(tool_name_ + ":1:", *param_);
    try
    {
      QString tmp_ini_file = File::getTempDirectory().toQString() + QDir::separator() + "TOPPAS_" + tool_name_.toQString() + "_";
      if (tool_type_ != "")
      {
        tmp_ini_file += tool_type_.toQString() + "_";
      }
      tmp_ini_file += File::getUniqueName().toQString() + "_tmp.ini";
      //store current parameters
      arg_param_.store(tmp_ini_file.toStdString());
      //restore other parameters that might be missing
      String call = String("\"") + File::getExecutablePath() + tool_name_ + "\"" + " -write_ini " + String(filename_) + " -ini " + String(tmp_ini_file);
      if (tool_type_ != "")
      {
        call += " -type " + tool_type_;
      }

      if (system(call.c_str()) != 0)
      {
        QMessageBox::critical(0, "Error", (String("Could not execute '") + call + "'!\n\nMake sure the TOPP tools are present in '" + File::getExecutablePath() + "', that you have permission to write to the temporary file path, and that there is space left in the temporary file path.").c_str());
        return;
      }
    }
    catch (Exception::BaseException & e)
    {
      QMessageBox::critical(this, "Error", (String("Error storing INI file: ") + e.getMessage()).c_str());
      return;
    }
  }

}
