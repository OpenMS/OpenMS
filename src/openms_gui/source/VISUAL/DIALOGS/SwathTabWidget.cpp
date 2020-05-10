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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/SwathTabWidget.h>
#include <ui_SwathTabWidget.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/DIALOGS/PythonModuleRequirement.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>

#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <QMessageBox>
#include <QProcess>
#include <QSignalBlocker>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    String getOSWExe()
    {
      return File::getExecutablePath() + "OpenSwathWorkflow";
    }

    QString getDefaultOutDir()
    {
      auto dir = QDir::homePath().append("/SwathWizardOut");
      if (!QDir().exists(dir)) QDir().mkpath(dir);
      return dir;
    }

    SwathTabWidget::SwathTabWidget(QWidget* parent) :
        QTabWidget(parent),
        ui(new Ui::SwathTabWidget)
    {
        ui->setupUi(this);
        
        auto py_selector = (PythonSelector*)ui->py_selector;

        auto py_pyprophet = (PythonModuleRequirement*)ui->py_pyprophet;
        py_pyprophet->setRequiredModules({"pyprophet", "stats"});
        py_pyprophet->setFreeText("In order to run PyProphet after OpenSWATH, the above modules need to be installed\n" \
                                  "Once they are available, the 'pyProphet' tab will become active and configurable.");
        py_pyprophet->setTitle("External: PyProphet tool");
        connect(py_selector, &PythonSelector::valueChanged, py_pyprophet, &PythonModuleRequirement::validate);
        
        // call once to update py_pyprophet canvas 
        // alternative: load latest data from .ini and set py_selector (will update py_pyprophet via above signal/slot)
        py_pyprophet->validate(py_selector->getLastPython().toQString());

        ui->input_tr->setFileFormatFilter("Transition sqLite file (*.pqp)");
        ui->input_iRT->setFileFormatFilter("Transition sqLite file (*.pqp)");
        ui->out_dir->setDirectory(getDefaultOutDir());

        // create a default config from OpenSwathWorkflow
        String tmp_file = File::getTemporaryFile();
        QProcess qp;
        qp.start(getOSWExe().toQString(), QStringList() << "-write_ini" << tmp_file.toQString());
        qp.waitForFinished();
        
        ParamXMLFile().load(tmp_file, swath_param_);
        swath_param_ = swath_param_.copy("OpenSwathWorkflow:1:", true);
        // parameters to show:
        StringList extract = {"mz_extraction_window", "rt_extraction_window", "threads"};
        
        for (const auto& name : extract) swath_param_wizard_.setValue(name, ""); // create a dummy param, just so we can use ::copySubset
        swath_param_wizard_ = swath_param_.copySubset(swath_param_wizard_);
                
        ui->list_editor->load(swath_param_wizard_);

        // keep the group of input widgets in sync with respect to their current-working-dir, when browsing for new files
        connect(ui->input_mzMLs, &InputFileList::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
        connect(ui->input_iRT, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
        connect(ui->input_tr, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
        connect(ui->input_swath_windows, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
    }

    SwathTabWidget::~SwathTabWidget()
    {
        delete ui;
    }

    String infileToOSW(const String& infile)
    {
      return FileHandler::swapExtension(File::basename(infile), FileTypes::OSW);
    }

    void SwathTabWidget::on_run_swath_clicked()
    {
      if (!checkInputReady_()) return;

      updateSwathParamFromWidgets_();
      Param tmp_param;
      tmp_param.insert("OpenSwathWorkflow:1:", swath_param_);
      String tmp_ini = File::getTemporaryFile();
      ParamXMLFile().store(tmp_ini, tmp_param);
      QProcess qp;
      StringList in_mzMLs = ui->input_mzMLs->getFilenames();
      ui->tab_run->setEnabled(false); // grey out the Wizard until OSW returns...
      writeLog_(QString("Starting OpenSwathWorkflow with %1 mzML file(s)").arg(in_mzMLs.size()), true);
      QString out_dir(ui->out_dir->dirNameValid() ?
        ui->out_dir->getDirectory() :
        getDefaultOutDir());
      for (const auto& mzML : in_mzMLs)
      {
        qp.start(getOSWExe().toQString(), QStringList() << "-ini" << tmp_ini.toQString() << "-in" << mzML.toQString() << "-out_osw" << out_dir + "/" + infileToOSW(mzML).toQString());
        const bool success = qp.waitForFinished(-1);
        if (qp.error() == QProcess::FailedToStart)
        {
          QMessageBox::critical(this, "Error", QString("Process '").append(getOSWExe().toQString()).append("' failed to start. Does it exist? Is it executable?"));
          return;
        }
        const QString external_sout(qp.readAllStandardOutput());
        const QString external_serr(qp.readAllStandardError());
        if (!external_sout.isEmpty()) writeLog_("Standard output: " + external_sout, true);
        if (!external_serr.isEmpty()) writeLog_("Standard error: " + external_serr);
        writeLog_(("Exit code: " + String(qp.exitCode()).toQString()));
      
        bool any_failure = (success == false || qp.exitStatus() != 0 || qp.exitCode() != 0);
        if (any_failure)
        {
          QMessageBox::critical(this, "Error", QString("Process '").append(getOSWExe().toQString()).append("' did not finish successfully. Please check the log."));
        }
      } // mzML loop
      ui->tab_run->setEnabled(true);
    }

    void SwathTabWidget::on_edit_advanced_parameters_clicked()
    {
      // refresh 'swath_param_' from data within the Wizards controls
      updateSwathParamFromWidgets_();
      Param tmp_param = swath_param_;

      // remove all input and output parameters from the user interface we are about to show
      StringList to_remove;
      for (Param::ParamIterator it = tmp_param.begin(); it != tmp_param.end(); ++it)
      {
        if (it->tags.count("input file") || it->tags.count("output file"))
        {
          to_remove.push_back(it->name); // do not remove right away.. does not work
        }
      }
      for (const auto& p : to_remove) 
      {
        tmp_param.remove(p);
        if (tmp_param.exists(p + "_type")) tmp_param.remove(p + "_type"); // for good measure of related input/output parameters
      }
      // show the parameters to the user
      String executable = File::getExecutablePath() + "INIFileEditor";
      String tmp_file = File::getTemporaryFile();
      ParamXMLFile().store(tmp_file, tmp_param);
      QProcess qp;
      qp.start(executable.toQString(), QStringList() << tmp_file.toQString());
      ui->tab_run->setEnabled(false); // grey out the Wizard until INIFileEditor returns...
      qp.waitForFinished(-1);
      ui->tab_run->setEnabled(true);
      ParamXMLFile().load(tmp_file, tmp_param);
      swath_param_.update(tmp_param, false);

      // refresh controls
      updateWidgetsfromSwathParam_();
    }

    void SwathTabWidget::updateSwathParamFromWidgets_()
    {
      // refresh 'swath_param_wizard_' which is linked into ParamEditor
      ui->list_editor->store();
      // ... and merge into main param
      swath_param_.update(swath_param_wizard_, false);

      Param tmp;
      // grab the files
      tmp.setValue("tr", ui->input_tr->getFilename());
      tmp.setValue("tr_irt", ui->input_iRT->getFilename());
      // do not set 'in' because it allows for one file only, while we have more and need to iterate manually
      String swath_windows = ui->input_swath_windows->getFilename();
      if (!swath_windows.empty()) tmp.setValue("swath_windows_file", swath_windows);
      // do not set '-out_osw' because we might have multiple -in's and have to iterate manually

      // update; do NOT write directly to swath_param_, because 'setValue(name, value)' will loose the description and the tags, i.e. input-file etc. We need this information though!
      swath_param_.update(tmp, false, false, true, true, OpenMS_Log_warn);
    }

    void SwathTabWidget::updateWidgetsfromSwathParam_()
    {
      swath_param_wizard_.update(swath_param_, false, false, true, false, OpenMS_Log_warn);
      ui->list_editor->load(swath_param_wizard_);
    }

    void SwathTabWidget::writeLog_(const QString& text, bool new_section)
    {
      if (new_section)
      {
        ui->log_text->append(QString(10, '#').append(QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss")).append(QString(10, '#')).append("\n"));
      }
      ui->log_text->append(text);
    }
    
    bool SwathTabWidget::checkInputReady_()
    {
      if (ui->input_mzMLs->getFilenames().empty())
      {
        QMessageBox::critical(this, "Error", "Input mzML file(s) are missing! Please provide at least one!");
        return false;
      }
      if (ui->input_tr->getFilename().isEmpty())
      {
        QMessageBox::critical(this, "Error", "Input file 'Transition Library' is missing! Please provide one!");
        return false;
      }
      if (ui->input_iRT->getFilename().isEmpty())
      {
        QMessageBox::critical(this, "Error", "Input file 'iRT Library' is missing! Please provide one!");
        return false;
      }

      // swath_windows_file is optional... no need to check

      return true;
    }

    void SwathTabWidget::broadcastNewCWD_(const QString& new_cwd)
    {
      QSignalBlocker blocker(this); // RAII to avoid infinite loop (setCWD signals updatedCWD which is connected to slot broadcastNewCWD_)
      ui->input_mzMLs->setCWD(new_cwd);
      ui->input_iRT->setCWD(new_cwd);
      ui->input_tr->setCWD(new_cwd);
      ui->input_swath_windows->setCWD(new_cwd);
    }

  }   //namespace Internal
} //namspace OpenMS




