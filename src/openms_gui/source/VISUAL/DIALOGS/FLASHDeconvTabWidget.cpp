// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/DIALOGS/FLASHDeconvTabWidget.h>
#include <OpenMS/VISUAL/DIALOGS/WizardHelper.h>
#include <ui_FLASHDeconvTabWidget.h>

#include <QDesktopServices>
#include <QMessageBox>
#include <QProcess>
#include <QProgressDialog>
#include <QSignalBlocker>
#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    template class WizardGUILock<FLASHDeconvTabWidget>;

    String getFLASHDeconvExe()
    {
      return File::findSiblingTOPPExecutable("FLASHDeconv");
    }

    QString getFDDefaultOutDir()
    {
      auto dir = QDir::homePath().append("/FLASHDeconvOut");
      if (!QDir().exists(dir))
        QDir().mkpath(dir);
      return dir;
    }

    FLASHDeconvTabWidget::FLASHDeconvTabWidget(QWidget* parent) :
        QTabWidget(parent),
        ui(new Ui::FLASHDeconvTabWidget),
        ep_([&](const String& out) { writeLog_(out.toQString()); },
            [&](const String& out) { writeLog_(out.toQString()); })
    {
      ui->setupUi(this);

      writeLog_(QString("Welcome to the Wizard!"), Qt::darkGreen, true);

      // keep the group of input widgets in sync with respect to their current-working-dir, when browsing for new files
      connect(ui->input_mzMLs, &InputFileList::updatedCWD, this, &FLASHDeconvTabWidget::broadcastNewCWD_);

      // check the "checkbox_spec" true (output files for masses per spectrum)
      ui->checkbox_spec->setCheckState(Qt::Checked);

      // param setting
      setWidgetsfromFDDefaultParam_();

      ui->out_dir->setDirectory(getFDDefaultOutDir());
    }

    FLASHDeconvTabWidget::~FLASHDeconvTabWidget()
    {
      delete ui;
    }

    String infileToFDoutput(const String& infile)
    {
      return FileHandler::swapExtension(File::basename(infile), FileTypes::TSV);
    }

    StringList FLASHDeconvTabWidget::getMzMLInputFiles() const
    {
      return ui->input_mzMLs->getFilenames();
    }

    void FLASHDeconvTabWidget::on_run_fd_clicked()
    {
      if (!checkFDInputReady_())
        return;

      WizardGUILock lock(this); // forbid user interaction

      // get parameter
      updateFLASHDeconvParamFromWidgets_();
      updateOutputParamFromWidgets_();
      Param fd_param;
      fd_param.insert("FLASHDeconv:1:", flashdeconv_param_);
      String tmp_ini = File::getTemporaryFile();
      StringList in_mzMLs = getMzMLInputFiles();
      writeLog_(QString("Starting FLASHDeconv with %1 mzML file(s)").arg(in_mzMLs.size()), Qt::darkGreen, true);

      QProgressDialog progress("Running FLASHDeconv ", "Abort ...", 0, (int)in_mzMLs.size(), this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0); // show immediately
      progress.setValue(0);
      int step = 0;

      for (const auto& mzML : in_mzMLs)
      {
        updateOutputParamFromPerInputFile(mzML.toQString());
        Param tmp_param = Param(fd_param);
        tmp_param.insert("FLASHDeconv:1:", flashdeconv_param_outputs_);

        ParamXMLFile().store(tmp_ini, tmp_param);

        auto r = ep_.run(this,
                         getFLASHDeconvExe().toQString(),
                         QStringList() << "-ini" << tmp_ini.toQString()
                                       << "-in" << mzML.toQString()
                                       << "-out" << getCurrentOutDir_() + "/" + infileToFDoutput(mzML).toQString(),
                         "",
                         true);
        if (r != ExternalProcess::RETURNSTATE::SUCCESS)
          break;
        if (progress.wasCanceled())
          break;
        progress.setValue(++step);
      } // mzML loop

      progress.close();
    }

    void FLASHDeconvTabWidget::on_edit_advanced_parameters_clicked()
    {
      // refresh 'flashdeconv_param_' from data within the Wizards controls
      updateFLASHDeconvParamFromWidgets_();

      Param tmp_param = flashdeconv_param_;

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
      flashdeconv_param_.update(tmp_param, false);
    }

    void FLASHDeconvTabWidget::on_open_output_directory_clicked()
    {
      QDesktopServices::openUrl(QUrl::fromLocalFile(getCurrentOutDir_()));
    }

    void FLASHDeconvTabWidget::updateFLASHDeconvParamFromWidgets_()
    {
      ui->list_editor->store();
    }

    void FLASHDeconvTabWidget::updateOutputParamFromWidgets_()
    {
      // refresh output params with default values
      flashdeconv_output_tags_.clear();

      // get checkbox results from the Wizard control
      if (ui->checkbox_spec->isChecked())
      {
        flashdeconv_output_tags_.push_back("out_spec");
      }
      if (ui->checkbox_mzml->isChecked())
      {
        flashdeconv_output_tags_.push_back("out_mzml");
        flashdeconv_output_tags_.push_back("out_annotated_mzml");
      }
      if (ui->checkbox_promex->isChecked())
      {
        flashdeconv_output_tags_.push_back("out_promex");
      }
      if (ui->checkbox_topfd->isChecked())
      {
        flashdeconv_output_tags_.push_back("out_topFD");
        flashdeconv_output_tags_.push_back("out_topFD_feature");
      }

      // optional FLASHIda support part
      if (ui->checkbox_readlogfile->isChecked())
      {
        flashdeconv_output_tags_.push_back("in_log");
      }
    }

    void FLASHDeconvTabWidget::updateOutputParamFromPerInputFile(const QString& input_file_name)
    {
      const Size max_ms_level = flashdeconv_param_.getValue("max_MS_level");
      std::string filepath_without_ext = getCurrentOutDir_().toStdString() + "/" + FileHandler::stripExtension(File::basename(input_file_name));

      for (const auto& param : flashdeconv_param_outputs_)
      {
        const std::string tag = param.name;

        std::string org_desc = param.description;
        auto org_tags = flashdeconv_param_outputs_.getTags(tag);

        // if this output format is requested by the user
        bool is_requested = false;
        if (!flashdeconv_output_tags_.empty() && std::find(flashdeconv_output_tags_.begin(), flashdeconv_output_tags_.end(), tag) != flashdeconv_output_tags_.end())
        {
          is_requested = true;
        }

        if (tag == "out_mzml" || tag == "out_annotated_mzml" || tag == "out_promex" || tag == "in_log") //  params having string values //  params having string values
        {
          // if not requested, set default value
          if (!is_requested)
          {
            flashdeconv_param_outputs_.setValue(tag, "", org_desc, org_tags);
            continue;
          }

          // if requested, set file path accordingly
          String out_path = filepath_without_ext;
          if (tag == "out_mzml")
          {
            out_path += "_deconv.mzML";
          }
          else if (tag == "out_annotated_mzml")
          {
            out_path += "_annotated.mzML";
          }
          else if (tag == "out_promex")
          {
            out_path += ".ms1ft";
          }
          else // (tag == "in_log")
          {
            String dir_path_only = File::path(input_file_name);
            String file_name_only = FileHandler::stripExtension(File::basename(input_file_name));
            out_path = dir_path_only + '/' + "IDALog_" + file_name_only + ".log";
          }
          flashdeconv_param_outputs_.setValue(tag, out_path, org_desc, org_tags);
        }
        else // Params with values as stringList
        {
          // if not requested, set default value
          if (!is_requested)
          {
            std::vector<std::string> tmp;
            flashdeconv_param_outputs_.setValue(tag, tmp, org_desc, org_tags);
            continue;
          }

          // if requested, set file path accordingly
          std::string out_extension = "";
          if (tag == "out_spec")
          {
            out_extension = ".tsv";
          }
          if (tag == "out_topFD")
          {
            out_extension = ".msalign";
          }
          if (tag == "out_topFD_feature")
          {
            out_extension = ".feature";
          }
          std::vector<std::string> files_paths;
          for (Size i = 0; i < max_ms_level; ++i)
          {
            files_paths.push_back(filepath_without_ext + "_ms" + std::to_string(i + 1) + out_extension);
          }
          flashdeconv_param_outputs_.setValue(tag, files_paths, org_desc, org_tags);
        }
      }
    }

    void FLASHDeconvTabWidget::setWidgetsfromFDDefaultParam_()
    {
      // create a default INI of FLASHDeconv
      String tmp_file = File::getTemporaryFile();
      if (ep_.run(this, getFLASHDeconvExe().toQString(), QStringList() << "-write_ini" << tmp_file.toQString(), "", true) != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        exit(1);
      }
      ParamXMLFile().load(tmp_file, flashdeconv_param_);
      flashdeconv_param_ = flashdeconv_param_.copy("FLASHDeconv:1:", true);

      // parameters to show in default mode : flashdeconv_param_wizard_
      flashdeconv_param_.remove("log");
      flashdeconv_param_.remove("no_progress");
      flashdeconv_param_.remove("debug");
      flashdeconv_param_.remove("in");
      flashdeconv_param_.remove("out");

      // parameters for different output format & in_log
      StringList out_params = {"out_spec", "out_annotated_mzml", "out_mzml", "out_promex", "out_topFD", "out_topFD_feature", "in_log"};
      for (const auto& name : out_params)
        flashdeconv_param_outputs_.setValue(name, ""); // create a dummy param, just so we can use ::copySubset
      flashdeconv_param_outputs_ = flashdeconv_param_.copySubset(flashdeconv_param_outputs_);

      // remove output format params from global parameter set
      for (const auto& name : out_params)
        flashdeconv_param_.remove(name);

      ui->list_editor->load(flashdeconv_param_);
    }

    QString FLASHDeconvTabWidget::getCurrentOutDir_() const
    {
      QString out_dir(ui->out_dir->dirNameValid() ? ui->out_dir->getDirectory() : getFDDefaultOutDir());
      return out_dir;
    }

    void FLASHDeconvTabWidget::writeLog_(const QString& text, const QColor& color, bool new_section)
    {
      QColor tc = ui->log_text->textColor();
      if (new_section)
      {
        ui->log_text->setTextColor(Qt::darkBlue);
        ui->log_text->append(QString(10, '#').append(QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss")).append(QString(10, '#')).append("\n"));
        ui->log_text->setTextColor(tc);
      }

      ui->log_text->setTextColor(color);
      ui->log_text->append(text);
      ui->log_text->setTextColor(tc); // restore old color
    }

    void FLASHDeconvTabWidget::writeLog_(const String& text, const QColor& color, bool new_section)
    {
      writeLog_(text.toQString(), color, new_section);
    }

    bool FLASHDeconvTabWidget::checkFDInputReady_()
    {
      if (ui->input_mzMLs->getFilenames().empty())
      {
        QMessageBox::critical(this, "Error", "Input mzML file(s) are missing! Please provide at least one!");
        return false;
      }

      return true;
    }

    void FLASHDeconvTabWidget::broadcastNewCWD_(const QString& new_cwd)
    {
      // RAII to avoid infinite loop (setCWD signals updatedCWD which is connected to slot broadcastNewCWD_)
      QSignalBlocker blocker1(ui->input_mzMLs);
      ui->input_mzMLs->setCWD(new_cwd);
    }
  } // namespace Internal
} // namespace OpenMS
