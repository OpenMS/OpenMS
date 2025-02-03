// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QMessageBox>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QRegularExpression>

#include <QSvgRenderer>
#include <map>

namespace OpenMS
{
  struct NameComponent
  {
    String prefix, suffix;
    int counter = -1;
    NameComponent() = default;

    NameComponent(const String& r_prefix, const String& r_suffix)
    : prefix(r_prefix),
      suffix(r_suffix)
    {}

    String toString() const
    {
      return (prefix + (counter != -1 ? String("_") + String(counter).fillLeft('0', 3) : String()) + "." + suffix);
    }
  };

  TOPPASToolVertex::TOPPASToolVertex()
    : TOPPASToolVertex("", "")
  {
  }

  TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type) :
    name_(name),
    type_(type)
  {
    brush_color_ = brush_color_.lighter(130); // make TOPP tools more white compared to all other nodes
    initParam_();
    connect(this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
    connect(this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
    connect(this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
    connect(this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
  }

  TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs) :
    TOPPASVertex(rhs),
    name_(rhs.name_),
    type_(rhs.type_),
    param_(rhs.param_),
    status_(rhs.status_),
    tool_ready_(rhs.tool_ready_)
    
  {
  }

  TOPPASToolVertex& TOPPASToolVertex::operator=(const TOPPASToolVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);

    param_ = rhs.param_;
    name_ = rhs.name_;
    type_ = rhs.type_;
    finished_ = rhs.finished_;
    status_ = rhs.status_;
    breakpoint_set_ = false;

    return *this;
  }

  std::unique_ptr<TOPPASVertex> TOPPASToolVertex::clone() const
  {
    return std::make_unique<TOPPASToolVertex>(*this);
  }

  bool TOPPASToolVertex::initParam_(const QString& old_ini_file)
  {
    // this is the only exception for writing directly to the tmpDir, instead of a subdir of tmpDir, as scene()->getTempDir() might not be available yet
    QString ini_file = File::getTemporaryFile().toQString();
    QString program = File::findSiblingTOPPExecutable(name_).toQString();
    QStringList arguments;
    arguments << "-write_ini" << ini_file;

    if (!type_.empty())
    {
      arguments << "-type";
      arguments << type_.toQString();
    }
    // allow for update using old parameters
    if (old_ini_file != "")
    {
      if (!File::exists(old_ini_file))
      {
        String msg = String("Could not open old INI file '") + old_ini_file + "'! File does not exist!";
        if (getScene_() && getScene_()->isGUIMode())
        {
          QMessageBox::critical(nullptr, "Error", msg.c_str());
        }
        else
        {
          OPENMS_LOG_ERROR << msg << std::endl;
        }
        tool_ready_ = false;
        return false;
      }
      arguments << "-ini" << old_ini_file;
    }

    // actually request the INI
    QProcess p;
    p.start(program, arguments);
    if (!p.waitForFinished(-1) || p.exitStatus() != 0 || p.exitCode() != 0)
    {
      String msg = String("Error! Call to '") + program + "' '" + String(arguments.join("' '")) +
          " returned with exit code (" + String(p.exitCode()) + "), exit status (" + String(p.exitStatus()) + ")." +
          "\noutput:\n" + String(QString(p.readAll())) +
          "\n";
      if (getScene_() && getScene_()->isGUIMode())
      {
        QMessageBox::critical(nullptr, "Error", msg.c_str());
      }
      else
      {
        OPENMS_LOG_ERROR << msg << std::endl;
      }
      tool_ready_ = false;
      return false;
    }
    if (!File::exists(ini_file))
    { // it would be weird to get here, since the TOPP tool ran successfully above, so INI file should exist, but nevertheless:
      String msg = String("Could not open '") + ini_file + "'! It does not exist!";
      if (getScene_() && getScene_()->isGUIMode())
      {
        QMessageBox::critical(nullptr, "Error", msg.c_str());
      }
      else
      {
        OPENMS_LOG_ERROR << msg << std::endl;
      }
      tool_ready_ = false;
      return false;
    }

    Param tmp_param;
    ParamXMLFile().load(String(ini_file).c_str(), tmp_param);
    // remember the parameters of this tool
    param_ = tmp_param.copy(name_ + ":1:", true); // get first instance (we never use more -- this is a legacy layer in paramXML)
    param_.setValue("no_progress", "true"); // by default, we do not want each tool to report loading/status statistics (would clutter the log window)
    // the user is free however, to re-enable it for individual nodes

    // write to disk to see if anything has changed
    writeParam_(param_, ini_file);
    bool changed = false;
    if (old_ini_file != "")
    {
      //check if INI file has changed (quick & dirty by file size)
      QFile q_ini(ini_file);
      QFile q_old_ini(old_ini_file);
      changed = q_ini.size() != q_old_ini.size();
    }
    setToolTip(String(param_.getSectionDescription(name_)).toQString());

    return changed;
  }

  void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
  {
    editParam();
  }

  void TOPPASToolVertex::editParam()
  {
    // use a copy for editing
    Param edit_param(param_);

    QVector<String> hidden_entries;
    // remove entries that are handled by edges already, user should not see them
    QVector<IOInfo> input_infos = getInputParameters();
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      int index = (*it)->getTargetInParam();
      if (index < 0)
      {
        continue;
      }

      const String& name = input_infos[index].param_name;
      if (edit_param.exists(name))
      {
        hidden_entries.push_back(name);
      }
    }

    QVector<IOInfo> output_infos = getOutputParameters();
    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
    {
      int index = (*it)->getSourceOutParam();
      if (index < 0)
      {
        continue;
      }

      const String& name = output_infos[index].param_name;
      if (edit_param.exists(name))
      {
        hidden_entries.push_back(name);
      }
    }

    // remove entries explained by edges
    foreach(const String &name, hidden_entries)
    {
      edit_param.remove(name);
    }

    // edit_param no longer contains tool description, take it from the node tooltip
    QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
    String default_dir;
    TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, name_, type_, toolTip(), hidden_entries);
    if (dialog.exec())
    {
      // take new values
      param_.update(edit_param);
      reset(true);
      emit parameterChanged(doesParamChangeInvalidate_());
    }

    getScene_()->updateEdgeColors();
  }

  TOPPASScene* TOPPASToolVertex::getScene_() const
  {
    return qobject_cast<TOPPASScene*>(scene());
  }

  bool TOPPASToolVertex::doesParamChangeInvalidate_()
  {
    return status_ == TOPPASToolVertex::TOOL_SCHEDULED || // all stati that will not tolerate a change in parameters
           status_ == TOPPASToolVertex::TOOL_RUNNING ||
           status_ == TOPPASToolVertex::TOOL_SUCCESS;
  }

  bool TOPPASToolVertex::invertRecylingMode()
  {
    allow_output_recycling_ = !allow_output_recycling_;
    emit parameterChanged(doesParamChangeInvalidate_()); // using 'true' is very conservative but safe. One could override this in child classes.
    return allow_output_recycling_;
  }

  QVector<TOPPASToolVertex::IOInfo> TOPPASToolVertex::getInputParameters() const
  {
    return getParameters_(true);
  }

  QVector<TOPPASToolVertex::IOInfo> TOPPASToolVertex::getOutputParameters() const
  {
    return getParameters_(false);
  }

  QVector<TOPPASToolVertex::IOInfo> TOPPASToolVertex::getParameters_(bool input_params) const
  {
    QVector<IOInfo> io_infos;
    auto add_params = [&](const String& search_tag) {
      for (Param::ParamIterator it = param_.begin(); it != param_.end(); ++it)
      {
        if (! it->tags.count(search_tag)) continue; // skip irrelevant parameters

        StringList valid_types(ListUtils::toStringList<std::string>(it->valid_strings));
        for (Size i = 0; i < valid_types.size(); ++i)
        {
          if (! valid_types[i].hasPrefix("*."))
          {
            std::cerr << "Invalid restriction \"" + valid_types[i] + "\"" + " for parameter \"" + it->name + "\"!" << std::endl;
            break;
          }
          valid_types[i] = valid_types[i].suffix(valid_types[i].size() - 2);
        }

        IOInfo io_info;
        io_info.param_name = it.getName();
        io_info.valid_types = valid_types;
        if (it->value.valueType() == ParamValue::STRING_LIST)
        { 
          io_info.type = IOInfo::IOT_LIST;
        }
        else if (it->value.valueType() == ParamValue::STRING_VALUE)
        {
          io_info.type = search_tag == TOPPBase::TAG_OUTPUT_DIR ?IOInfo::IOT_DIR : IOInfo::IOT_FILE;
        }
        else { std::cerr << "TOPPAS: Unexpected parameter value!" << std::endl; }
        io_infos.push_back(io_info);
      }
    };
    if (input_params)
    {
      add_params(TOPPBase::TAG_INPUT_FILE);
    }
    else
    {
      add_params(TOPPBase::TAG_OUTPUT_FILE);
      add_params(TOPPBase::TAG_OUTPUT_DIR);
    }

    // order in param can change --> sort
    std::sort(io_infos.begin(), io_infos.end());
    return io_infos;
  }

  void TOPPASToolVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget, false);

    QString draw_str = (type_.empty() ? name_ : name_ + " (" + type_ + ")").toQString();
    for (int i = 0; i < 10; ++i)
    {
      QString prev_str = draw_str;
      draw_str = toolnameWithWhitespacesForFancyWordWrapping_(painter, draw_str);
      if (draw_str == prev_str)
      {
        break;
      }
    }

    QRectF text_boundings = painter->boundingRect(QRectF(-65, -35, 130, 70), Qt::AlignCenter | Qt::TextWordWrap, draw_str);
    painter->drawText(text_boundings, Qt::AlignCenter | Qt::TextWordWrap, draw_str);

    if (status_ != TOOL_READY)
    {
      QString text = QString::number(round_counter_) + " / " + QString::number(round_total_);

      QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
      painter->drawText((int)(62.0 - text_boundings.width()), 48, text);
    }

    // progress light
    painter->setPen(Qt::black);
    QColor progress_color;
    switch (status_)
    {
    case TOOL_READY:
      progress_color = Qt::lightGray; break;

    case TOOL_SCHEDULED:
      progress_color = Qt::darkBlue; break;

    case TOOL_RUNNING:
      progress_color = Qt::yellow; break;

    case TOOL_SUCCESS:
      progress_color = Qt::green; break;

    case TOOL_CRASH:
      progress_color = Qt::red; break;

    default:
      progress_color = Qt::magenta; break; // signal weird status by color
    }
    painter->setBrush(progress_color);
    painter->drawEllipse(46, -52, 14, 14);

    // breakpoint set?
    if (breakpoint_set_)
    {
      QSvgRenderer* svg_renderer = new QSvgRenderer(QString(":/stop_sign.svg"), nullptr);
      painter->setOpacity(0.35);
      svg_renderer->render(painter, QRectF(-60, -60, 120, 120));
    }
  }

  QString TOPPASToolVertex::toolnameWithWhitespacesForFancyWordWrapping_(QPainter* painter, const QString& str)
  {
    qreal max_width = 130;
    QStringList parts = str.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
    QStringList new_parts;

    for(const QString& part : parts)
    {
      QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter | Qt::TextWordWrap, part);
      if (text_boundings.width() <= max_width)
      {
        //word not too long
        new_parts.append(part);
      }
      else
      {
        //word too long -> insert space at reasonable position -> Qt::TextWordWrap can break the line there
        int last_capital_index = 1;
        for (int i = 1; i <= part.size(); ++i)
        {
          QString tmp_str = part.left(i);
          //remember position of last capital letter
          if (tmp_str.at(i - 1).isUpper())
          {
            last_capital_index = i;
          }
          QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter | Qt::TextWordWrap, tmp_str);
          if (text_boundings.width() > max_width)
          {
            //break line at next capital letter before this position
            new_parts.append(part.left(last_capital_index - 1) + "-");
            new_parts.append(part.right(part.size() - last_capital_index + 1));
            break;
          }
        }
      }
    }

    return new_parts.join(" ");
  }

  QRectF TOPPASToolVertex::boundingRect() const
  {
    return QRectF(-71, -61, 142, 122);
  }

  String TOPPASToolVertex::getName() const
  {
    return name_;
  }

  const String& TOPPASToolVertex::getType() const
  {
    return type_;
  }

  void TOPPASToolVertex::run()
  {
    __DEBUG_BEGIN_METHOD__

    //check if everything ready (there might be more than one upstream node - ALL need to be ready)
    if (!isUpstreamFinished())
    {
      return;
    }

    if (finished_)
    {
      OPENMS_LOG_ERROR << "This should not happen. Calling an already finished node '" << this->name_ << "' (#" << this->getTopoNr() << ")!" << std::endl;
      throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    TOPPASScene* ts = getScene_();

    QString ini_file = ts->getTempDir()
                       + QDir::separator()
                       + getOutputDir().toQString()
                       + QDir::separator()
                       + name_.toQString();
    if (!type_.empty())
    {
      ini_file += "_" + type_.toQString();
    }
    // do not write the ini yet - we might need to alter it

    RoundPackages pkg;
    String error_msg("");
    bool success = buildRoundPackages(pkg, error_msg);
    if (!success)
    {
      OPENMS_LOG_ERROR << "Could not retrieve input files from upstream nodes...\n";
      emit toolFailed(-1, error_msg.toQString());
      return;
    }

    // all inputs are ready --> GO!
    if (!updateCurrentOutputFileNames(pkg, error_msg)) // based on input, we prepare output names
    {
      emit toolFailed(-1, error_msg.toQString());
      return;
    }

    createDirs();

    //emit toolStarted(); //disabled! Every signal emitted here does only mean the process is queued(!) not that its executed right away

    /// update round status
    round_total_ = (int) pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0; // once round_counter_ reaches round_total_, we are done

    QStringList shared_args;
    if (!type_.empty())
    {
      shared_args << "-type" << type_.toQString();
    }
    // get *all* input|output file parameters (regardless if edge exists)
    QVector<IOInfo> in_params = getInputParameters(), out_params = getOutputParameters();

    bool ini_round_dependent = false; // indicates if we need a new INI file for each round (usually GenericWrapper issue)

    // maximum number of filenames per TOPP parameter file-list to put on the commandline
    // If more filenames are needed, e.g. for MapAligner's -in/-out etc., they are put in the .INI file
    // to avoid exceeding the 8KB length limit of the Windows commandline
    static constexpr int MAX_FILES_CMDLINE {10};

    for (int round = 0; round < round_total_; ++round)
    {
      debugOut_(String("Enqueueing process nr ") + round + "/" + round_total_);
      QStringList args = shared_args;

      // we might need to modify input/output file parameters before storing to INI
      Param param_tmp = param_;

      /// INCOMING EDGES
      for (RoundPackageConstIt ite = pkg[round].begin();
           ite != pkg[round].end();
           ++ite)
      {
        TOPPASEdge incoming_edge = *(ite->second.edge);

        int param_index = incoming_edge.getTargetInParam();
        if (param_index < 0 || param_index >= in_params.size())
        {
          OPENMS_LOG_ERROR << "TOPPAS: Input parameter index out of bounds!" << std::endl;
          return;
        }

        String param_name = in_params[param_index].param_name;

        const QStringList& file_list = ite->second.filenames.get();

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        // OR if there are a lot of input files (which might exceed the 8k length limit of cmd.exe on Windows)
        if (param_name.hasPrefix("ETool:") || file_list.size() > MAX_FILES_CMDLINE)
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }

        if (!store_to_ini)
        {
          args << "-" + param_name.toQString() << file_list;
        }
        else
        {
          if (param_tmp.getValue(param_name).valueType() == ParamValue::STRING_LIST)
          {
            param_tmp.setValue(param_name, ListUtils::create<std::string>(StringListUtils::fromQStringList(file_list)));
          }
          else
          {
            if (file_list.size() > 1)
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            }
            param_tmp.setValue(param_name, String(file_list[0]));
          }
        }
      }

      // OUTGOING EDGES (output files and output folders)
      // ;output names are already prepared by 'updateCurrentOutputFileNames()'
      typedef RoundPackage::iterator EdgeIndexIt;
      for (EdgeIndexIt it_edge  = output_files_[round].begin();
           it_edge != output_files_[round].end();
           ++it_edge)
      {
        int param_index = it_edge->first;
        String param_name = out_params[param_index].param_name;

        bool store_to_ini = false;
        
        const QStringList& output_files = output_files_[round][param_index].filenames.get();
        
        // check for GenericWrapper input/output files and put them in INI file:
        // OR if there are a lot of input files (which might exceed the 8k length limit of cmd.exe on Windows)
        if (param_name.hasPrefix("ETool:") || output_files.size() > MAX_FILES_CMDLINE)
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }

        
        if (!store_to_ini)
        {
          args << "-" + param_name.toQString() << output_files;
        }
        else
        {
          if (param_tmp.getValue(param_name).valueType() == ParamValue::STRING_LIST)
          {
            param_tmp.setValue(param_name, ListUtils::create<std::string>(StringListUtils::fromQStringList(output_files)));
          }
          else
          {
            if (output_files.size() > 1)
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            }
            param_tmp.setValue(param_name, String(output_files[0]));
          }
        }
      }

      // each iteration might have different params (input/output items which are registered in subsections (GenericWrapper stuff))
      QString ini_file_iteration;
      if (ini_round_dependent)
      {
        ini_file_iteration = QDir::toNativeSeparators(ini_file + QString::number(round) + ".ini");
      }
      else
      {
        ini_file_iteration = QDir::toNativeSeparators(ini_file + ".ini");
      }
      writeParam_(param_tmp, ini_file_iteration);
      args << "-ini" << ini_file_iteration;

      // create process
      QProcess* p;
      if (!ts->isDryRun())
      {
        p = new QProcess();
      }
      else
      {
        p = new FakeProcess();
      }

      p->setProcessChannelMode(QProcess::MergedChannels);
      connect(p, SIGNAL(readyReadStandardOutput()), this, SLOT(forwardTOPPOutput()));
      connect(ts, SIGNAL(terminateCurrentPipeline()), p, SLOT(kill()));
      // let this node know that round is done
      connect(p, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(executionFinished(int, QProcess::ExitStatus)));

      // enqueue process
      String msg_enqueue = String("\nEnqueue: \"") + File::getExecutablePath() + name_ + "\" \"" + String(args.join("\" \"")) + "\"\n";
      if (round == 0)
      {
        // active if TOPPAS is run with --debug; will print to console
        OPENMS_LOG_DEBUG << msg_enqueue << std::endl;
        // show sys-call in logWindow of TOPPAS (or console for non-gui)
        if ((int) param_tmp.getValue("debug") > 0)
        {
          ts->logTOPPOutput(msg_enqueue.toQString());
        }
      }
      toolScheduledSlot();
      ts->enqueueProcess(TOPPASScene::TOPPProcess(p, File::findSiblingTOPPExecutable(name_).toQString(), args, this));
    }

    // run pending processes
    ts->runNextProcess();

    __DEBUG_END_METHOD__
  }

  void TOPPASToolVertex::emitToolStarted()
  {
    emit toolStarted();
  }

  void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
  {
    __DEBUG_BEGIN_METHOD__

    TOPPASScene* ts = getScene_();
    QProcess* p = qobject_cast<QProcess*>(QObject::sender());

    RAIICleanup clean([&]() {
      // clean up at end
      if (p)
      {
        delete p;
      }

      ts->processFinished();
    });

    //** ERROR handling
    if (es != QProcess::NormalExit)
    {
      emit toolCrashed();
    }
    else if (ec != 0)
    {
      emit toolFailed(ec);
    }
    else
    {
      //** no error ... proceed
      ++round_counter_;
      //std::cout << (String("Increased iteration_nr_ to ") + round_counter_ + " / " + round_total_ ) << " for " << this->name_ << std::endl;

      if (round_counter_ == round_total_) // all iterations performed --> proceed in pipeline
      {
        debugOut_("All iterations finished!");

        if (finished_)
        {
          OPENMS_LOG_ERROR << "SOMETHING is very fishy. The vertex is already set to finished, yet there was still a thread spawning..." << std::endl;
          throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        if (!ts->isDryRun())
        {
          renameOutput_(); // rename generated files by content
          emit toolFinished();
        }
        finished_ = true;

        if (!breakpoint_set_)
        {
          // call all children, proceed in pipeline
          for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
          {
            TOPPASVertex* tv = (*it)->getTargetVertex();
            debugOut_(String("Starting child ") + tv->getTopoNr());
            tv->run();
          }
          debugOut_("All children started!");
        }
      }
    }

   

    __DEBUG_END_METHOD__
  }

  bool TOPPASToolVertex::renameOutput_()
  {
    // get all output names
    QStringList files = this->getFileNames();

    std::map<String, NameComponent> name_old_to_new;
    std::map<String, int> name_new_count, name_new_idx; // count occurrence (for optional counter infix)

    // a first round to find which filenames are not unique (and require augmentation with a counter)

    for (const QString& file : files)
    {
      if (File::isDirectory(file)) continue; // skip output directories

      String new_prefix = FileHandler::stripExtension(file);
      String new_suffix = FileTypes::typeToName(FileHandler::getTypeByContent(file)); // this might replace bla.fasta with bla.FASTA ... which is the same file on Windows
      if (file.endsWith(new_suffix.toQString(), Qt::CaseInsensitive)) // --> use the native suffix (to avoid deleting the source file when renaming)
      {
        new_suffix = String(file).suffix(new_suffix.size());
      }
      NameComponent nc(new_prefix, new_suffix);
      name_old_to_new[file] = nc;
      ++name_new_count[nc.toString()];
    }
    // for all names which occur more than once, introduce a counter  
    for (const QString& file : files)
    {
      if (name_new_count[name_old_to_new[file].toString()] > 1) // candidate for counter
      {
        name_old_to_new[file].counter = ++name_new_idx[name_old_to_new[file].toString()]; // start at index 1
      }
    }


    for (Size i = 0; i < output_files_.size(); ++i)
    {
      for (RoundPackageIt it = output_files_[i].begin();
           it != output_files_[i].end();
           ++it)
      {
        for (int fi = 0; fi < it->second.filenames.size(); ++fi)
        {
          // skip output directories
          if (File::isDirectory(it->second.filenames[fi]))
          { 
            continue;
          }

          // rename file and update record
          String old_filename = QDir::toNativeSeparators(it->second.filenames[fi]);
          String new_filename = QDir::toNativeSeparators(name_old_to_new[it->second.filenames[fi]].toString().toQString());
          if (QFileInfo(old_filename.toQString()).canonicalFilePath() == QFileInfo(new_filename.toQString()).canonicalFilePath())
          { // source and target are identical -- no action required
            continue;
          }
          QFile file(old_filename.toQString());
          if (File::exists(new_filename))
          { // rename only works if the target file does not exist: delete it first
            bool success = File::remove(new_filename);
            if (!success)
            {
              OPENMS_LOG_ERROR << "Could not remove '" << new_filename << "'.\n";
              throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_filename);
            }
          }
          bool success = file.rename(new_filename.toQString());
          if (!success)
          {
            OPENMS_LOG_ERROR << "Could not rename '" << String(it->second.filenames[fi]) << "' to '" << new_filename << "'\n";
            throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_filename);
          }
          it->second.filenames.set(new_filename.toQString(), fi);
        }
      }
    }
    return true;
  }

  const Param& TOPPASToolVertex::getParam()
  {
    return param_;
  }

  void TOPPASToolVertex::setParam(const Param& param)
  {
    param_ = param;
  }

  TOPPASToolVertex::TOOLSTATUS TOPPASToolVertex::getStatus() const
  {
    return status_;
  }

  bool TOPPASToolVertex::updateCurrentOutputFileNames(const RoundPackages& pkg, String& error_msg)
  {
    const auto number_of_rounds = pkg.size();
    if (pkg.empty())
    {
      error_msg = "Less than one round received from upstream tools. Something is fishy!\n";
      OPENMS_LOG_ERROR << error_msg;
      return false;
    }

    QVector<IOInfo> out_params = getOutputParameters();
    // check if this tool outputs a list of files, or only single files
    bool has_only_singlefile_output = !IOInfo::isAnyList(out_params);

    // look for the input with the most files in round 0 (as this is the maximal number of output files we can produce)
    // we assume the number of files is equal in all rounds...
    int max_size_index(-1);
    int max_size(-1);

    // iterate over input edges
    for (RoundPackageConstIt it  = pkg[0].begin();
          it != pkg[0].end();
          ++it)
    {
      if (it->second.edge->getSourceVertex()->isRecyclingEnabled())
      { // skip recycling input nodes
        continue;
      }

      // we only need to find a good upstream node with a single file -- since we only output single files
      if (has_only_singlefile_output)
      { // .. take any non-recycled input edge, preferably from 'in' and/or single inputs
        if ((max_size < 1 || (it->second.edge->getTargetInParamName() == "in") || it->second.filenames.size() == 1))
        {
          max_size_index = it->first;
          max_size       = 1;
        }

      }
      else if ((it->second.filenames.size() > max_size) ||   // either just larger 
          // ... or it's from '-in' (which we prefer as naming source).. only for non-recycling -in though
          ((it->second.filenames.size () == max_size) && (it->second.edge->getTargetInParamName() == "in")))
      {
        max_size_index = it->first;
        max_size       = it->second.filenames.size();
      }
    }

    if (max_size_index == -1)
    {
      error_msg = "Did not find upstream nodes with un-recycled names. Something is fishy!\n";
      OPENMS_LOG_ERROR << error_msg;
      return false;
    }


    // now we construct output filenames for this node
    // use names from the selected upstream vertex (hoping that this is the maximal number of files we are going to produce)
    std::vector<QStringList> per_round_basenames;
    for (Size i = 0; i < number_of_rounds; ++i)
    {
      QStringList filenames = pkg[i].find(max_size_index)->second.filenames.get();
      //
      // remove suffix to avoid chaining .mzML.idxml.tsv
      // a new suffix is added later, depending on edge-type etc
      //
      // try to find the type (only by looking at the suffix); not doing it manually, since it could be .mzXML.gz
      for (QString& filename : filenames)
      {
        filename = FileHandler::stripExtension(filename).toQString();
      }
      per_round_basenames.push_back(filenames);
      //std::cerr << "  output filenames (round " << i  <<"): " << per_round_basenames.back().join(", ") << std::endl;
    }

    // maybe we find something more unique, e.g. last base directory if all filenames are equal
    smartFileNames_(per_round_basenames);

    // clear output file list
    output_files_.clear();
    output_files_.resize(number_of_rounds);

    const TOPPASScene* ts = getScene_();
    
    // output names for each outgoing edge
    for (int i = 0; i < out_params.size(); ++i)
    {
      // search for an out edge for this parameter (not required to exist)
      int param_index;
      TOPPASEdge* edge_out(nullptr);
      for (ConstEdgeIterator it_edge = outEdgesBegin(); it_edge != outEdgesEnd(); ++it_edge)
      {
        param_index = (*it_edge)->getSourceOutParam();
        if (i == param_index) // corresponding out edge found
        {
          edge_out = *it_edge;
          break;
        }
      }
      if (!edge_out)
      {
        continue;
      }

      // determine output file format if possible (for suffix)
      String file_suffix;
      if (out_params[i].type == IOInfo::IOT_DIR)
      {
        file_suffix = "_dir"; // we need something non-empty
      }
      // Single file or list of files
      else if (out_params[i].valid_types.size() == 1)
      { // only one format allowed
        auto t = FileTypes::nameToType(out_params[i].valid_types[0]);
        if (t != FileTypes::UNKNOWN) 
        { // only use canonical names for suffixes, i.e. exact upper/lowercase match, i.e. always use "consensusXML", not "ConsensusXML" or "cOnsensusXML"
          // (TOPPAS requires this, to avoid errors when copying temporary files)
          file_suffix = "." + FileTypes::typeToName(t);
        }
        else
        { // unknown type... just use it as it is
          file_suffix = "." + out_params[i].valid_types[0];
        }
      }
      else if (String p_out_format = out_params[i].param_name + "_type"; // expected parameter name which determines output format
               param_.exists(p_out_format))
      { // 'out_type' or alike is specified
        if (!param_.getValue(p_out_format).toString().empty())
        {
          file_suffix = "." + param_.getValue(p_out_format).toString();
        }
        else
        {
          OPENMS_LOG_WARN << "TOPPAS cannot determine output file format for param '" << out_params[i].param_name
                             << "' of Node " + this->name_ + "(" + String(this->getTopoNr()) + "). Format is ambiguous. Use parameter '" + p_out_format + "' to name intermediate output correctly!\n";
        }
      }
      if (file_suffix.empty())
      { 
        // Are we FileMerger? If so we can recover the out_type from the type of our input files
        if (name_ == "FileMerger")
        {
          // For this very specific case we know that all the upstream nodes have to have the same types
          file_suffix = "." + param_.getValue("in_type").toString();
        }
        // tag as unknown (TOPPAS will try to rename the output file once its written - see renameOutput_())
        else 
        {
          OPENMS_LOG_DEBUG << " unknown extension for : " << out_params[i].param_name << " in: " << name_ <<"\n";
          file_suffix = ".unknown";
        }
      }
      //std::cerr << "suffix is: " << file_suffix << "\n\n";

      // create common path of output files
      QString path = ts->getTempDir()
                     + QDir::separator()
                     + getOutputDir().toQString() // includes TopoNr
                     + QDir::separator()
                     + out_params[param_index].param_name.remove(':').toQString().left(50) // max 50 chars per subdir
                     + QDir::separator();

      VertexRoundPackage vrp;
      vrp.edge = edge_out;

      std::set<QString> filename_output_set; // verify that output files are unique (avoid overwriting)
      assert(per_round_basenames.size() == number_of_rounds);
      for (Size round = 0; round < number_of_rounds; ++round)
      {
        // store edge for this param for all rounds
        output_files_[round][param_index] = vrp; // index by index of source-out param

        // list --> single file (e.g. IDMerger or FileMerger)
        bool list_to_single = (per_round_basenames[round].size() > 1 && out_params[param_index].type == IOInfo::IOT_FILE);
        for (const QString &input_file : per_round_basenames[round])
        {
          QString fn = path + QFileInfo(input_file).fileName(); // out_path + filename
          OPENMS_LOG_DEBUG << "Single:" << fn.toStdString() << "\n";
          if (out_params[param_index].type == IOInfo::IOT_DIR)
          { // output is a directory
            fn = QDir::toNativeSeparators(path);
            if (number_of_rounds > 1)
            { // use a different output folder for each round if multiple rounds are present
              fn += QFileInfo(input_file).baseName(); 
            }
              
            output_files_[round][param_index].filenames.push_back(fn);
            OPENMS_LOG_DEBUG << "Dir:" << fn.toStdString() << "\n";
            break; // only one iteration required (there is only one output dir per output param, irrespective of #input files)
          }
          else if (list_to_single)
          {
            if (fn.contains(QRegularExpression(".*_to_.*_mrgd")))
            {
              fn = fn.left(fn.indexOf("_to_"));
              OPENMS_LOG_DEBUG << "  first merge in merge: " << fn.toStdString() << "\n";
            }
            QString fn_last = QFileInfo(per_round_basenames[round].last()).fileName();
            if (fn_last.contains(QRegularExpression(".*_to_.*_mrgd")))
            {
              int i_start = fn_last.indexOf("_to_") + 4;
              fn_last = fn_last.mid(i_start, fn_last.indexOf("_mrgd", i_start) - i_start);
              OPENMS_LOG_DEBUG << "  last merge in merge: " << fn_last.toStdString() << "\n";
            }
            fn += "_to_" + fn_last + "_mrgd";
            OPENMS_LOG_DEBUG << "  List: ..." << "_to_" + fn_last.toStdString() + "_mrgd" << "\n";
          }
          if (!fn.endsWith(file_suffix.toQString()))
          {
            fn += file_suffix.toQString();
            OPENMS_LOG_DEBUG << "  Suffix-add: " << file_suffix << "\n";
          }
          fn = QDir::toNativeSeparators(fn);
          output_files_[round][param_index].filenames.push_back(fn);
          if (list_to_single)
          {
            break; // only one iteration required
          }
          if (auto [it, newly_inserted] = filename_output_set.insert(fn); !newly_inserted)
          {
            error_msg = "TOPPAS failed to build correct filenames. Please report this bug, along with your Pipeline\n!";
            OPENMS_LOG_ERROR << error_msg;
            return false;
          }
        }
      } // end for rounds
          
      //std::cerr << "output filenames (" << out_params[i].param_name <<") final: " << ListUtils::concatenate< std::set<QString> >(filename_output_set, ", ") << std::endl;
    } // end for out params (each edge)

    return true;
  }

  void TOPPASToolVertex::smartFileNames_(std::vector<QStringList>& filenames)
  {
    /* TODO:
     * implement this carefully; also take care of what happens after the call
     * of this method in updateCurrentOutputFileNames()
     */

    // special case #1, only one filename in each round (at least 2 rounds), with different directory but same basename
    // --> use LAST directory as new name, e.g. 'subdir' from 'c:\mydir\subdir\samesame.mzML'
    bool passes_constraints = false;
    if (filenames.size() > 1) // more than one round
    {
      passes_constraints = true;
      for (Size i = 1; i < filenames.size(); ++i)
      {
        if ((filenames[i].size() > 1) // one file per round AND unique filename
           || (QFileInfo(filenames[0][0]).fileName() != QFileInfo(filenames[i][0]).fileName()))
        {
          passes_constraints = false;
          break;
        }
      }
    }

    if (passes_constraints) // rename
    {
      for (Size i = 0; i < filenames.size(); ++i)
      {
        QString p = QDir::toNativeSeparators(QFileInfo(filenames[i][0]).canonicalPath());
        if (p.isEmpty())
        {
          continue;
        }
        //std::cout << "PATH: " << p << "\n";
        String tmp = String(p).suffix(String(QString(QDir::separator()))[0]);
        //std::cout << "INTER: " << tmp << "\n";
        if (tmp.size() <= 2 || tmp.has(':'))
        {
          continue; // too small to be reliable; might even be 'c:'
        }
        filenames[i][0] = tmp.toQString();
        //std::cout << "  -->: " << filenames[i][0] << "\n";
      }
      return; // we do not want the next special case on top of this...
    }

    // possibilities for more good naming schemes...
    // special case #2 ...

    return;
  }

  void TOPPASToolVertex::forwardTOPPOutput()
  {
    QProcess* p = qobject_cast<QProcess*>(QObject::sender());
    if (!p)
    {
      return;
    }

    QString out = p->readAllStandardOutput();
    emit toppOutputReady(out);
  }

  void TOPPASToolVertex::toolStartedSlot()
  {
    status_ = TOOL_RUNNING;
    update(boundingRect());
  }

  void TOPPASToolVertex::toolFinishedSlot()
  {
    status_ = TOOL_SUCCESS;
    update(boundingRect());
  }

  void TOPPASToolVertex::toolScheduledSlot()
  {
    status_ = TOOL_SCHEDULED;
    update(boundingRect());
  }

  void TOPPASToolVertex::toolFailedSlot()
  {
    status_ = TOOL_CRASH;
    update(boundingRect());
  }

  void TOPPASToolVertex::toolCrashedSlot()
  {
    status_ = TOOL_CRASH;
    update(boundingRect());
  }

  void TOPPASToolVertex::inEdgeHasChanged()
  {
    // something has changed --> tmp files might be invalid --> reset
    reset(true);
    TOPPASVertex::inEdgeHasChanged();
  }

  void TOPPASToolVertex::outEdgeHasChanged()
  {
    // something has changed --> tmp files might be invalid --> reset
    reset(true);
    TOPPASVertex::outEdgeHasChanged();
  }

  void TOPPASToolVertex::openContainingFolder() const
  {
    QString path = getFullOutputDirectory().toQString();
    GUIHelpers::openFolder(path);
  }

  String TOPPASToolVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = getScene_();
    return QDir::toNativeSeparators(ts->getTempDir() + QDir::separator() + getOutputDir().toQString());
  }

  String TOPPASToolVertex::getOutputDir() const
  {
    TOPPASScene* ts = getScene_();
    String workflow_dir = FileHandler::stripExtension(File::basename(ts->getSaveFileName()));
    if (workflow_dir.empty())
    {
      workflow_dir = "Untitled_workflow";
    }
    String dir = workflow_dir +
                 String(QDir::separator()) +
                 get3CharsNumber_(topo_nr_) + "_" + getName();
    if (!getType().empty())
    {
      dir += "_" + getType();
    }

    return dir;
  }

  void TOPPASToolVertex::createDirs()
  {
    QDir dir;
    if (!dir.mkpath(getFullOutputDirectory().toQString()))
    {
      OPENMS_LOG_ERROR << "TOPPAS: Could not create path " << getFullOutputDirectory() << std::endl;
    }

    // subsdirectories named after the output parameter name
    QStringList files = this->getFileNames();
    foreach(const QString &file, files)
    {
      QString sdir = File::path(file).toQString();
      if (!File::exists(sdir))
      {
        if (!dir.mkpath(sdir))
        {
          OPENMS_LOG_ERROR << "TOPPAS: Could not create path " << String(sdir) << std::endl;
        }
      }
    }
  }

  void TOPPASToolVertex::setTopoNr(UInt nr)
  {
    if (topo_nr_ != nr)
    {
      // topological number changes --> output dir changes --> reset
      reset(true);
      topo_nr_ = nr;
      emit somethingHasChanged();
    }
  }

  void TOPPASToolVertex::reset(bool reset_all_files)
  {
    __DEBUG_BEGIN_METHOD__

    finished_ = false;
    status_ = TOOL_READY;
    output_files_.clear();

    if (reset_all_files)
    {
      QString remove_dir = getFullOutputDirectory().toQString();
      if (File::exists(remove_dir))
      {
        File::removeDirRecursively(remove_dir);
      }
    }

    TOPPASVertex::reset(reset_all_files);

    __DEBUG_END_METHOD__
  }

  bool TOPPASToolVertex::refreshParameters()
  {
    TOPPASScene* ts = getScene_();
    QString old_ini_file = ts->getTempDir() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
    if (!type_.empty())
    {
      old_ini_file += type_.toQString() + "_";
    }
    old_ini_file += File::getUniqueName().toQString() + "_tmp_OLD.ini";
    writeParam_(param_, old_ini_file);

    bool changed = initParam_(old_ini_file);
    QFile::remove(old_ini_file);

    return changed;
  }

  bool TOPPASToolVertex::isToolReady() const
  {
    return tool_ready_;
  }

  void TOPPASToolVertex::writeParam_(const Param& param, const QString& ini_file)
  {
    Param save_param;
    save_param.setValue(name_ + ":1:toppas_dummy", "blub");
    save_param.insert(name_ + ":1:", param);
    save_param.remove(name_ + ":1:toppas_dummy");
    save_param.setSectionDescription(name_ + ":1", "Instance '1' section for '" + name_ + "'");
    ParamXMLFile paramFile;
    paramFile.store(ini_file, save_param);
  }

  void TOPPASToolVertex::toggleBreakpoint()
  {
    breakpoint_set_ = !breakpoint_set_;
  }

}
