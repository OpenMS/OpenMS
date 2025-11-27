// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/SwathTabWidget.h>
#include <ui_SwathTabWidget.h>
#include <OpenMS/VISUAL/DIALOGS/WizardHelper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/DIALOGS/PythonModuleRequirement.h>

#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <QMessageBox>
#include <QProcess>
#include <QProgressDialog>
#include <QSignalBlocker>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    template class WizardGUILock<SwathTabWidget>;

    String getOSWExe()
    {
      return File::findSiblingTOPPExecutable("OpenSwathWorkflow");
    }

    QString getDefaultOutDir()
    {
      auto dir = QDir::homePath().append("/SwathWizardOut");
      if (!QDir().exists(dir)) QDir().mkpath(dir);
      return dir;
    }

    SwathTabWidget::SwathTabWidget(QWidget* parent) :
        QTabWidget(parent),
        ui(new Ui::SwathTabWidget),
        ep_([&](const String& out) {writeLog_(out.toQString());},
            [&](const String& out) {writeLog_(out.toQString());})
    {
      ui->setupUi(this);

      writeLog_(QString("Welcome to the Wizard!"), Qt::darkGreen, true);

      auto py_selector = (PythonSelector*)ui->py_selector;

      auto py_pyprophet = (PythonModuleRequirement*)ui->py_pyprophet;
      py_pyprophet->setRequiredModules( { "pyprophet", "msproteomicstoolslib" });
      py_pyprophet->setFreeText("In order to run PyProphet and TRIC after OpenSWATH, the above modules need to be installed\n" \
                                "Once they are available, the 'PyProphet and TRIC' tab will become active and configurable.\n" \
                                "To install the modules, visit the <a href=\"http://www.openswath.org\">openswath.org homepage</a> and follow the installation instructions.");
      py_pyprophet->setTitle("External: PyProphet and TRIC tools");
      connect(py_selector, &PythonSelector::valueChanged, py_pyprophet, &PythonModuleRequirement::validate);
        
      // call once to update py_pyprophet canvas 
      // alternative: load latest data from .ini and set py_selector (will update py_pyprophet via above signal/slot)
      py_pyprophet->validate(py_selector->getLastPython().toQString());

      ui->input_tr->setFileFormatFilter("Transition sqLite file (*.pqp)");
      ui->input_iRT->setFileFormatFilter("Transition sqLite file (*.pqp)");
        
      // create a default INI of OpenSwathWorkflow
      String tmp_file = File::getTemporaryFile();
      if (ep_.run(this, getOSWExe().toQString(), QStringList() << "-write_ini" << tmp_file.toQString(), "", true) != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        exit(1);
      }
        
      ParamXMLFile().load(tmp_file, swath_param_);
      swath_param_ = swath_param_.copy("OpenSwathWorkflow:1:", true);
      // parameters to show within the Wizard:
      StringList extract = {"mz_extraction_window", "rt_extraction_window", "threads"};
        
      for (const auto& name : extract) swath_param_wizard_.setValue(name, ""); // create a dummy param, just so we can use ::copySubset
      swath_param_wizard_ = swath_param_.copySubset(swath_param_wizard_);
                
      ui->list_editor->load(swath_param_wizard_);

      // keep the group of input widgets in sync with respect to their current-working-dir, when browsing for new files
      connect(ui->input_mzMLs, &InputFileList::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
      connect(ui->input_iRT, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
      connect(ui->input_tr, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
      connect(ui->input_swath_windows, &InputFile::updatedCWD, this, &SwathTabWidget::broadcastNewCWD_);
      
      // update information on assay libraries
      connect(ui->input_iRT, &InputFile::updatedFile, ui->input_iRT_stats, &SwathLibraryStats::updateFromFile);
      connect(ui->input_tr, &InputFile::updatedFile, ui->input_tr_stats, &SwathLibraryStats::updateFromFile);


      // if out_dir (user defined output directory for OSW results) changes ...
      connect(ui->out_dir, &OutputDirectory::directoryChanged, this, &SwathTabWidget::checkPyProphetInput_);
      // ... or the mzML input changes --> update the input list of tbl_py_osws (which depend in the basenames of mzMLs and the OSW output directory)
      connect(ui->input_mzMLs, &InputFileList::updatedCWD, this, &SwathTabWidget::checkPyProphetInput_);

      // prepare table for pyProphet input files
      auto& tbl = *(ui->tbl_py_osws);
      tbl.setHeaders(QStringList() << "filename");
      tbl.setSortingEnabled(false);
      connect(ui->tbl_py_osws, &TableView::resized, [&]() {
        ui->tbl_py_osws->setColumnWidth(0, ui->tbl_py_osws->width()); // full size column
      });

      ui->out_dir->setDirectory(getDefaultOutDir());
    }

    SwathTabWidget::~SwathTabWidget()
    {
      delete ui;
    }

    String infileToOSW(const String& infile)
    {
      return FileHandler::swapExtension(File::basename(infile), FileTypes::OSW);
    }

    String infileToChrom(const String& infile)
    {
      return FileHandler::swapExtension(File::basename(infile), FileTypes::SQMASS);
    }

    StringList SwathTabWidget::getMzMLInputFiles() const
    {
      return ui->input_mzMLs->getFilenames();
    }

    void SwathTabWidget::on_run_swath_clicked()
    {
      if (!checkOSWInputReady_()) return;

      WizardGUILock lock(this); // forbid user interaction

      updateSwathParamFromWidgets_();
      Param tmp_param;
      tmp_param.insert("OpenSwathWorkflow:1:", swath_param_);
      String tmp_ini = File::getTemporaryFile();
      ParamXMLFile().store(tmp_ini, tmp_param);
      StringList in_mzMLs = getMzMLInputFiles();
      writeLog_(QString("Starting OpenSwathWorkflow with %1 mzML file(s)").arg(in_mzMLs.size()), Qt::darkGreen, true);
      
      QProgressDialog progress("Running OpenSwath", "Abort ...", 0, (int)in_mzMLs.size(), this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0); // show immediately
      progress.setValue(0);
      int step = 0;

      writeLog_(QString("Running OpenSwathWorkflow (%1 files total): ").arg(in_mzMLs.size()), Qt::darkGreen, true);
      for (const auto& mzML : in_mzMLs)
      {
        auto r = ep_.run(this, 
                         getOSWExe().toQString(), 
                         QStringList() << "-ini" << tmp_ini.toQString() 
                                       << "-in" << mzML.toQString()
                                       << "-out_osw" << getCurrentOutDir_() + "/" + infileToOSW(mzML).toQString()
                                       << "-out_chrom" << getCurrentOutDir_() + "/" + infileToChrom(mzML).toQString(),
                         "",
                         true);
        if (r != ExternalProcess::RETURNSTATE::SUCCESS) break;
        if (progress.wasCanceled()) break;
        progress.setValue(++step);
      } // mzML loop
      
      progress.close();
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
      tmp.setValue("tr", ui->input_tr->getFilename().toStdString());
      tmp.setValue("tr_irt", ui->input_iRT->getFilename().toStdString());
      // do not set 'in' because it allows for one file only, while we have more and need to iterate manually
      String swath_windows = ui->input_swath_windows->getFilename();
      if (!swath_windows.empty()) tmp.setValue("swath_windows_file", swath_windows);
      // do not set '-out_osw' because we might have multiple -in's and have to iterate manually

      // call update(); do NOT write directly to swath_param_ using 'setValue(name, value)' because that will loose the description and the tags, i.e. input-file etc. We need this information though!
      swath_param_.update(tmp, false, false, true, true, OpenMS_Log_warn);
    }

    void SwathTabWidget::updateWidgetsfromSwathParam_()
    {
      swath_param_wizard_.update(swath_param_, false, false, true, false, OpenMS_Log_warn);
      ui->list_editor->load(swath_param_wizard_);
    }

    QString SwathTabWidget::getCurrentOutDir_() const
    {
      QString out_dir(ui->out_dir->dirNameValid() ?
        ui->out_dir->getDirectory() :
        getDefaultOutDir());
      return out_dir;
    }

    vector<pair<String, bool>> SwathTabWidget::getPyProphetInputFiles() const
    {
      vector<pair<String, bool>> files;
      String dir = getCurrentOutDir_();
      for (const auto& file : getMzMLInputFiles())
      {
        // predict output OSW filenames
        const String file_osw = dir + '/' + infileToOSW(file);
        // check if exists
        files.emplace_back(file_osw, File::exists(file_osw));
      }
      return files;
    }

    const char* msg_no_osws = "select mzML input files in 'LC-MS files' tab first and pick an output directory in 'Run OpenSwath' tab";

    void SwathTabWidget::checkPyProphetInput_()
    {
      // populate the file list widget for pyProphet input
      auto& tbl = *(ui->tbl_py_osws);
      tbl.setRowCount(0); // clear the table; clear() does not resize the table!
      auto files = getPyProphetInputFiles();
      if (files.empty())
      {
        tbl.appendRow();
        tbl.setAtBottomRow(msg_no_osws, 0, Qt::white, Qt::gray);
      }
      else 
      {
        for (const auto& file : files)
        {
          tbl.appendRow();
          tbl.setAtBottomRow(file.first.c_str(), 0, Qt::white, 
                             file.second ? Qt::black : Qt::red)
              ->setCheckState(Qt::Unchecked);
        }
      }
      // set label at bottom
      ui->lbl_pyOutDir->setText("Results can be found in '" + getCurrentOutDir_() + 
                                "'. If pyProphet ran, there will be PDF files with model statistics and TRIC will "
                                "generate TSV files (tric_aligned.tsv and tric_aligned_matrix.tsv) for downstream processing.\n To view results interactively, open them in TOPPView.");
    }

    void SwathTabWidget::writeLog_(const QString& text, const QColor& color, bool new_section)
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

    void SwathTabWidget::writeLog_(const String& text, const QColor& color, bool new_section)
    {
      writeLog_(text.toQString(), color, new_section);
    }
    
    bool SwathTabWidget::checkOSWInputReady_()
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
      // RAII to avoid infinite loop (setCWD signals updatedCWD which is connected to slot broadcastNewCWD_)
      QSignalBlocker blocker1(ui->input_mzMLs);
      QSignalBlocker blocker2(ui->input_iRT);
      QSignalBlocker blocker3(ui->input_tr);
      QSignalBlocker blocker4(ui->input_swath_windows);
      ui->input_mzMLs->setCWD(new_cwd);
      ui->input_iRT->setCWD(new_cwd);
      ui->input_tr->setCWD(new_cwd);
      ui->input_swath_windows->setCWD(new_cwd);
    }

    bool SwathTabWidget::findPythonScript_(const String& path_to_python_exe, String& script_name)
    {
      String path = File::path(path_to_python_exe);
      String script_backup = script_name;
      script_name = path + "/Scripts/" + script_backup; // Windows uses the Script subdirectory
      if (File::readable(script_name)) return true;
      writeLog_("Warning: Could not find " + script_backup + " at " + script_name + ".", Qt::red, true);
      script_name = path + "/" + script_backup;
      if (File::readable(script_name)) return true;
      writeLog_("Warning: Could not find " + script_backup + " at " + script_name + ".", Qt::red, true);
      return false;
    }

    QStringList SwathTabWidget::getPyProphetOutputFileNames() const
    {
      auto in = getPyProphetInputFiles();
      QStringList out;
      for (auto& f : in)
      {
        auto s = (FileHandler::stripExtension(f.first) + "_pyProphet_out.osw");
        out << s.toQString();
      }
      return out;
    }
    


    void SwathTabWidget::on_btn_runPyProphet_clicked()
    {
      if (!ui->py_pyprophet->isReady())
      {
        QMessageBox::warning(this, "Error", "Could not find all requirements for 'pyprophet & tric' (see 'Config' tab). Install modules via 'pip install <modulename>' and make sure it's available in $PATH");
        return;
      }
      WizardGUILock lock(this); // forbid user interaction

      auto inputs = getPyProphetInputFiles();
      if (inputs.empty())
      {
        QMessageBox::warning(this, "Error", "Provide at least one input file for pyProphet and TRIC in the 'LC-MS files' tab.");
        return;
      }
      QStringList osws = getPyProphetOutputFileNames();
      QStringList osws_orig;
      QStringList osws_reduced;
      QStringList tsvs;
      for (const auto& file : inputs)
      {
        if (file.second == false)
        {
          QMessageBox::warning(this, "Error", String("Required input file '" + file.first + "' not found. Please run OpenSwathWorkflow first to create it").toQString());
          return;
        }
        osws_orig << file.first.toQString();
        osws_reduced << (FileHandler::stripExtension(file.first) + "_pyProphet.oswr").toQString();
        tsvs << (FileHandler::swapExtension(file.first, FileTypes::TSV)).toQString();
      }
      // check presence of template
      QString library = ui->input_tr->getFilename();
      if (library.isEmpty())
      {
        QMessageBox::warning(this, "Error", String("The assay library is not specified. Please go to the 'database' tab and specify it.").toQString());
        return;
      }
#ifdef OPENMS_WINDOWSPLATFORM
      String pp = "pyprophet.exe"; // we need the full path for findPythonScript_
#else
      String pp = "pyprophet";
#endif
      if (!findPythonScript_(ui->py_selector->getLastPython(), pp)) // searches Script in Python installation
      {
        QMessageBox::warning(this, "Error", String("Could not find 'pyprophet' in the python installation '" + ui->py_selector->getLastPython() + "'. Please make sure it is installed. Visit http://openswath.org/en/latest/docs/tric.html for details.").toQString());
        return;
      }
      // list of calls to make: exe, args, [optional] list of args to append one-by-one in a loop
      std::vector<Command> calls;
      // merge all osws ...
      calls.emplace_back(pp, QStringList() << "merge" << "--template=" + library << "--out=model.osw" << osws, ArgLoop{});
      // to build/learn a common model --> creates merged_ms1ms2_report.pdf
      calls.emplace_back(pp, QStringList() << "score" << "--in=model.osw" << "--level=ms1ms2", ArgLoop{});
      // apply model in loop
      calls.emplace_back(pp, QStringList() << "score" << "--apply_weights=model.osw" << "--level=ms1ms2" << "--in" << "%1", ArgLoop{ Args{osws, 4} });
      // reduce (required to avoid https://github.com/PyProphet/pyprophet/issues/85)
      calls.emplace_back(pp, QStringList() << "reduce" << "--in" << "%1" << "--out" << "%1", ArgLoop{ Args{osws, 2}, Args{osws_reduced, 4} });
      // merge again for peptide and protein error rate control
      calls.emplace_back(pp, QStringList() << "merge" << "--template=model.osw" << "--out=model_global.osw" << osws_reduced, ArgLoop{});
      calls.emplace_back(pp, QStringList() << "peptide" << "--in=model_global.osw" << "--context=global", ArgLoop{});
      calls.emplace_back(pp, QStringList() << "protein" << "--in=model_global.osw" << "--context=global", ArgLoop{});
      // backpropagate in loop
      calls.emplace_back(pp, QStringList() << "backpropagate" << "--apply_scores=model_global.osw" << "--in" << "%1", ArgLoop{ Args{osws, 3} });
      // prepare for TRIC
      calls.emplace_back(pp, QStringList() << "export" << "--format=legacy_merged" << "--max_global_peptide_qvalue=0.01" << "--max_global_protein_qvalue=0.01" 
                                            << "--in=%1" << "--out=%1", ArgLoop{ Args{osws, 4}, Args{tsvs, 5} });
      String feature_alignment_py = "feature_alignment.py";
      if (!findPythonScript_(ui->py_selector->getLastPython(), feature_alignment_py)) // searches Script in Python installation
      {
        QMessageBox::warning(this, "Error", String("Could not find 'feature_alignment.py' from the msproteomicstool package in the python installation '" + ui->py_selector->getLastPython() + "'. Please make sure it is installed. Visit http://openswath.org/en/latest/docs/tric.html for details.").toQString());
        return;
      }
      calls.emplace_back(ui->py_selector->getLastPython(), QStringList() << feature_alignment_py.toQString() << "--in" << tsvs
                                            << "--out" << "tric_aligned.tsv" << "--out_matrix" << "tric_aligned_matrix.tsv" 
                                            << "--method" << "LocalMST" << "--realign_method" << "lowess" << "--max_rt_diff" << "90" 
                                            << "--fdr_cutoff" << QString::number(ui->tric_FDR_threshold->value()) 
                                            << "--alignment_score" << QString::number(ui->tric_RTmax->value()), ArgLoop{});
        
      QProgressDialog progress("Running pyprophet and TRIC", "Abort ...", 0, (int)calls.size(), this);
      progress.setWindowModality(Qt::ApplicationModal);
      progress.setMinimumDuration(0); // show immediately
      progress.setValue(0);

      // first - copy all original osw files, since augmenting them once with model information will lead to crashes when doing a second run on them
      for (int i = 0; i < osws_orig.size(); ++i)
      {
        QFile::remove(osws[i]); // copy() will not overwrite existing files :/
        QFile::copy(osws_orig[i], osws[i]);
      }

      int step = 0;
      for (const auto& call : calls)
      { 
        // this might just be one loop... depending on the call...
        for (size_t i_loop = 0; i_loop < call.getLoopCount(); ++i_loop)
        {
          auto returnstate = ep_.run(this, call.exe.toQString(), call.getArgs(i_loop), getCurrentOutDir_(), true);
          if (returnstate != ExternalProcess::RETURNSTATE::SUCCESS)
          {
            QMessageBox::warning(this, "Error", String("Running pyprophet/TRIC failed at step " + String(step) + "/" + String(calls.size()) + ". Please see log for details").toQString());
            return;
          }
          if (progress.wasCanceled())
          {
            return;
          }
        }

        progress.setValue(++step);
      }
      
      progress.close();
    }

    void SwathTabWidget::on_btn_pyresults_clicked()
    {
      GUIHelpers::openFolder(getCurrentOutDir_());
    }

    void SwathTabWidget::on_pushButton_clicked()
    {
      auto& tbl = *(ui->tbl_py_osws);
      int selected_rows = 0;
      QStringList missing_osw_files;
      QStringList args;
      auto raw_files = getMzMLInputFiles(); // mzML's for now; will be translated to sqMass below
      auto osw_files = getPyProphetOutputFileNames();
      if (tbl.rowCount() == 1 && tbl.item(0, 0)->data(Qt::DisplayRole).toString() == msg_no_osws)
      {
        QMessageBox::information(this, "Error", "No files are selected from the list above! Make sure to select mzML files in the 'LC-MS files' tab first.");
        return;
      }
      if (size_t(tbl.rowCount()) != raw_files.size() || tbl.rowCount() != osw_files.count())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Something went wrong in populating the input file window");
      }
      for (int i = 0; i < tbl.rowCount(); ++i)
      {
        if (tbl.item(i, 0)->checkState() == Qt::CheckState::Checked)
        {
          ++selected_rows;
          args << infileToChrom(raw_files[i]).toQString() << "!" << osw_files[i];
          if (!File::exists(osw_files[i])) missing_osw_files << File::basename(osw_files[i]).toQString();
        }
      }
      if (selected_rows == 0)
      {
        QMessageBox::information(this, "Error", "No files are selected from the list above! Select the files you want to open and try again.");
        return;
      }
      if ( ! missing_osw_files.isEmpty())
      {
        QMessageBox::information(this, "Error", "The following selected files to not yet have a pyProphet result file:\n" + missing_osw_files.join("\n") + "\nPlease run pyProphet first");
        return;
      }

      if (QMessageBox::question(this, "Confirm", (String("Confirm opening ") + selected_rows + " raw files in TOPPView").toQString(), 
                                QMessageBox::StandardButton::Ok, QMessageBox::StandardButton::Cancel) 
                             == QMessageBox::StandardButton::Ok)
      {
        // create on heap, to avoid QProcess being closed when application exits 
        // This is a small memleak, but QProcess::startDetached() is only available from Qt 5.10 on -- switch to that once its supported everywhere
        QProcess* qp = new QProcess;
        qp->setWorkingDirectory(getCurrentOutDir_());
        auto tv = File::findSiblingTOPPExecutable("TOPPView");
        qp->start(tv.toQString(), args);
        if (!qp->waitForStarted(2000))
        {
          QMessageBox::warning(this, "Error", String("Could not open TOPPView executable from '" + tv + "'").toQString());
          return;
        }
        // TOPPView is running now ... detached
      }
    }

  }   //namespace Internal
} //namspace OpenMS




