// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/VISUAL/DIALOGS/FLASHQuantTabWidget.h>
#include <OpenMS/VISUAL/DIALOGS/WizardHelper.h>
#include <ui_FLASHQuantTabWidget.h>

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
    template class WizardGUILock<FLASHQuantTabWidget>;

    String getFLASHQuantExe()
    {
      return File::findSiblingTOPPExecutable("FLASHQuant");
    }

    QString getFQDefaultOutDir()
    {
      auto dir = QDir::homePath().append("/FLASHQuantOut");
      if (!QDir().exists(dir)) QDir().mkpath(dir);
      return dir;
    }

    String getTopDownConsensusFeatureGroupExe()
    {
      return File::findSiblingTOPPExecutable("TopDownConsensusFeatureGroup");
    }

    FLASHQuantTabWidget::FLASHQuantTabWidget(QWidget* parent) :
        QTabWidget(parent),
        ui(new Ui::FLASHQuantTabWidget),
        ep_([&](const String& out) {writeLog_(out.toQString());},
            [&](const String& out) {writeLog_(out.toQString());})
    {
      ui->setupUi(this);

      writeLog_(QString("Welcome to the Wizard!"), Qt::darkGreen, true);

      // keep the group of input widgets in sync with respect to their current-working-dir, when browsing for new files
      connect(ui->input_mzMLs, &InputFileList::updatedCWD, this, &FLASHQuantTabWidget::broadcastNewCWD_);

      // param setting
      setWidgetsfromFQDefaultParam_();

      ui->out_dir->setDirectory(getFQDefaultOutDir());
    }

    FLASHQuantTabWidget::~FLASHQuantTabWidget()
    {
      delete ui;
    }

    QString FLASHQuantTabWidget::infileToFQoutput(const String& infile, const String& extension) const
    {
      String file_name = FileHandler::stripExtension(File::basename(infile)) + ".fq." + extension;
      return getCurrentOutDir_() + "/" + file_name.toQString();
    }

    StringList FLASHQuantTabWidget::getMzMLInputFiles() const
    {
      return ui->input_mzMLs->getFilenames();
    }

    void FLASHQuantTabWidget::on_run_fq_clicked()
    {
      if (!checkFQInputReady_())
        return;

      WizardGUILock lock(this); // forbid user interaction

      // get parameter
      updateFLASHQuantParamFromWidgets_();
      updateOutputParamFromWidgets_();
      Param fq_param;
      fq_param.insert("FLASHQuant:1:", flashquant_param_);

      String tmp_ini = File::getTemporaryFile();

      StringList in_mzMLs = getMzMLInputFiles();
      writeLog_(QString("Starting FLASHQuant with %1 mzML file(s)").arg(in_mzMLs.size()), Qt::darkGreen, true);

      QProgressDialog progress("Running FLASHQuant ", "Abort ...", 0, (int)in_mzMLs.size(), this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0); // show immediately
      progress.setValue(0);
      int step = 0;

      for (const auto& mzML : in_mzMLs)
      {
        Param tmp_param = Param(fq_param);
        ParamXMLFile().store(tmp_ini, tmp_param);
        QStringList full_param_string = QStringList() << "-ini" << tmp_ini.toQString()
                                                      << "-in" << mzML.toQString()
                                                      << "-out" << infileToFQoutput(mzML, "tsv");
        if (featurexml_output_)
        {
          full_param_string << "-out_feat" << infileToFQoutput(mzML, "featureXML");
        }

        auto r = ep_.run(this,
                         getFLASHQuantExe().toQString(),
                         full_param_string,
                         "",
                         true);
        if (r != ExternalProcess::RETURNSTATE::SUCCESS) break;
        if (progress.wasCanceled()) break;
        progress.setValue(++step);
      } // mzML loop

      progress.close();

      /// consensus feature group calculation
      if (ui->checkbox_consensus->isChecked())
      {
        runTopDownConsensusFeatureGroup_();
      }
    }

    void FLASHQuantTabWidget::on_open_output_directory_clicked()
    {
      QDesktopServices::openUrl( QUrl::fromLocalFile(getCurrentOutDir_()) );
    }

    void FLASHQuantTabWidget::on_consensus_names_check_clicked()
    {
      // TODO: if checkbox is checked, edit parameter part comes down...

      /// check if Checkbox is clicked & the prefix is given
      if( !ui->checkbox_consensus->isChecked() )
      {
        QMessageBox qmsg_box;
        QString msg = "Consensus output was not requested. Please check the checkbox to request.";
        qmsg_box.critical(nullptr, "Error", msg);
        qmsg_box.show();
        return;
      }
      if( ui->consensus_name_edit->text().isEmpty() )
      {
        QMessageBox qmsg_box;
        QString msg = "Prefix of replicate is not given. Consensus feature group will be reported from all input LS-MS files.\n"
                      "Place your cursor over \"Prefix of replicate\" to check the details.";
        qmsg_box.warning(nullptr, "Error", msg);
        qmsg_box.show();
        return;
      }

      String given_prefix = ui->consensus_name_edit->text();

      /// find replicate groups
      std::vector<std::vector<String>> file_groups = groupReplicateFiles_(given_prefix);

      /// output message
      QMessageBox qmsg_box;
      String msg = "Given prefix: " + given_prefix;
      msg += "\nTheses files in each group will be considered as replicate:\n";
      for (Size i = 0; i < file_groups.size(); ++i)
      {
        msg+="Group" + to_string(i+1) + "\n";
        for (auto &f : file_groups[i])
        {
          msg+=f + "\n";
        }
      }
      qmsg_box.setText(msg.toQString());
      qmsg_box.exec();
    }

    void FLASHQuantTabWidget::updateFLASHQuantParamFromWidgets_()
    {
      ui->list_editor->store();
    }

    void FLASHQuantTabWidget::updateOutputParamFromWidgets_()
    {
      // refresh output params with default values
      featurexml_output_ = false;

      // get checkbox results from the Wizard control
      if( ui->checkbox_featxml->isChecked() )
      {
        featurexml_output_ = true;
      }
    }

    void FLASHQuantTabWidget::setWidgetsfromFQDefaultParam_()
    {
      // create a default INI of FLASHQuant
      String tmp_file = File::getTemporaryFile();
      if (ep_.run(this, getFLASHQuantExe().toQString(), QStringList() << "-write_ini" << tmp_file.toQString(), "", true) != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        exit(1);
      }
      ParamXMLFile().load(tmp_file, flashquant_param_);
      flashquant_param_ = flashquant_param_.copy("FLASHQuant:1:", true);

      // parameters to show in default mode : flashquant_param_wizard_
      flashquant_param_.remove("log");
      flashquant_param_.remove("no_progress");
      flashquant_param_.remove("debug");
      flashquant_param_.remove("in");
      flashquant_param_.remove("out");

      // parameters for feature_xml
      flashquant_param_.remove("out_feat");
      featurexml_output_ = false;

      ui->list_editor->load(flashquant_param_);
    }

    QString FLASHQuantTabWidget::getCurrentOutDir_() const
    {
      QString out_dir(ui->out_dir->dirNameValid() ?
        ui->out_dir->getDirectory() :
        getFQDefaultOutDir());
      return out_dir;
    }

    void FLASHQuantTabWidget::writeLog_(const QString& text, const QColor& color, bool new_section)
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

    void FLASHQuantTabWidget::writeLog_(const String& text, const QColor& color, bool new_section)
    {
      writeLog_(text.toQString(), color, new_section);
    }

    bool FLASHQuantTabWidget::checkFQInputReady_()
    {
      if (ui->input_mzMLs->getFilenames().empty())
      {
        QMessageBox::critical(this, "Error", "Input mzML file(s) are missing! Please provide at least one!");
        return false;
      }

      return true;
    }

    std::vector<std::vector<String>> FLASHQuantTabWidget::groupReplicateFiles_(String delimiter)
    {
      std::vector<String> candidates(ui->input_mzMLs->getFilenames());
      std::vector<std::vector<String>> file_groups;

      // if no delimiter is given, put all LC-MS inputs together in one group
      if (delimiter.empty())
      {
        file_groups.push_back(candidates);
        return file_groups;
      }

      while(candidates.size() > 0)
      {
        // start with the first element
        String reference_file = candidates[0];
        String prefix = reference_file.substr(0, reference_file.find(delimiter));

        candidates.erase(candidates.begin()); // remove the reference

        // find the replicates
        std::vector<String> collected;
        for (auto& file : candidates)
        {
          if (file.rfind(prefix, 0) == 0)
          {
            collected.push_back(file);
          }
        }

        // if no replicates are found, finish
        if (collected.size() == 0)
        {
          continue;
        }

        // remove collected files from candidates
        for (auto& f : collected)
        {
          candidates.erase(std::remove(candidates.begin(), candidates.end(), f), candidates.end());
        }

        collected.push_back(reference_file);
        std::sort(collected.begin(), collected.end());
        file_groups.push_back(collected);
      }
      return file_groups;
    }

    void FLASHQuantTabWidget::runTopDownConsensusFeatureGroup_()
    {
      writeLog_(QString("Starting TopDownConsensusFeatureGroup..."), Qt::darkGreen, true);

      /// find replicate groups
      String given_prefix = ui->consensus_name_edit->text();
      std::vector<std::vector<String>> file_groups = groupReplicateFiles_(given_prefix);

      /// small dialog in front of main window
      QProgressDialog progress("Running TopDownConsensusFeatureGroup ", "Abort ...", 0, (int)file_groups.size(), this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0); // show immediately
      progress.setValue(0);
      int step = 0;

      // run TopDownConsensusFeatureGroup
      for (auto& group : file_groups)
      {
        // get FDQ result file names
        QStringList fdq_results;
        for (auto &mzml : group)
        {
          fdq_results.push_back(infileToFQoutput(mzml, "tsv"));
        }

        // output file name
        String tmp_infile = infileToFQoutput(group[0], "tsv");
        String output_prefix = tmp_infile.substr(0, tmp_infile.find(given_prefix));
        String consensus_path = output_prefix + ".fq.consensus.tsv";

        // run
        QStringList params = QStringList() << "-in" << fdq_results << "-out" << consensus_path.toQString();
        auto r = ep_.run(this,
                         getTopDownConsensusFeatureGroupExe().toQString(),
                         params,
                         "",
                         true);
        if (r != ExternalProcess::RETURNSTATE::SUCCESS) break;
        if (progress.wasCanceled()) break;
        progress.setValue(++step);
      }
      progress.close();
    }

    void FLASHQuantTabWidget::broadcastNewCWD_(const QString& new_cwd)
    {
      // RAII to avoid infinite loop (setCWD signals updatedCWD which is connected to slot broadcastNewCWD_)
      QSignalBlocker blocker1(ui->input_mzMLs);
      ui->input_mzMLs->setCWD(new_cwd);
    }
  }   //namespace Internal
} //namspace OpenMS




