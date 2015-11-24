// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QRegExp>

#include <QSvgRenderer>

namespace OpenMS
{


  struct NameComponent
  {
    String prefix, suffix;
    int counter;
    NameComponent()
      : counter(-1)
    {}

    NameComponent(const String& r_prefix, const String& r_suffix)
      : prefix(r_prefix),
      suffix(r_suffix),
      counter(-1)
    {}

    String toString() const
    {
      String s_counter;
      if (counter != -1) s_counter = String(counter).fillLeft('0', 3) + ".";
      return (prefix + s_counter + suffix);
    }

  };

  UInt TOPPASToolVertex::uid_ = 1;

  TOPPASToolVertex::TOPPASToolVertex() :
    TOPPASVertex(),
    name_(),
    type_(),
    param_(),
    status_(TOOL_READY),
    tool_ready_(true),
    breakpoint_set_(false)
  {
    pen_color_ = Qt::black;
    brush_color_ = QColor(245, 245, 245);
    initParam_();
    connect(this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
    connect(this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
    connect(this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
    connect(this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
  }

  TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type) :
    TOPPASVertex(),
    name_(name),
    type_(type),
    param_(),
    tool_ready_(true),
    breakpoint_set_(false)
  {
    pen_color_ = Qt::black;
    brush_color_ = QColor(245, 245, 245);
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
    tool_ready_(rhs.tool_ready_),
    breakpoint_set_(false)
  {
    pen_color_ = Qt::black;
    brush_color_ = QColor(245, 245, 245);
    connect(this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
    connect(this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
    connect(this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
    connect(this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
  }

  TOPPASToolVertex::~TOPPASToolVertex()
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

  bool TOPPASToolVertex::initParam_(const QString& old_ini_file)
  {
    Param tmp_param;
    // this is the only exception for writing directly to the tmpDir, instead of a subdir of tmpDir, as scene()->getTempDir() might not be available yet
    QString ini_file = File::getTempDirectory().toQString() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
    if (type_ != "")
    {
      ini_file += type_.toQString() + "_";
    }
    ini_file += File::getUniqueName().toQString() + "_tmp.ini";
    ini_file = QDir::toNativeSeparators(ini_file);

    QString program = File::findExecutable(name_).toQString();
    QStringList arguments;
    arguments << "-write_ini";
    arguments << ini_file;

    if (type_ != "")
    {
      arguments << "-type";
      arguments << type_.toQString();
    }
    // allow for update using old parameters
    if (old_ini_file != "")
    {
      if (!File::exists(old_ini_file))
      {
        QMessageBox::critical(0, "Error", (String("Could not open '") + old_ini_file + "'!").c_str());
        tool_ready_ = false;
        return false;
      }
      arguments << "-ini";
      arguments << old_ini_file;
    }

    // actually request the INI
    QProcess p;
    p.start(program, arguments);
    if (!p.waitForFinished(-1))
    {
      QMessageBox::critical(0, "Error", (String("Could not execute '") + program + " " + String(arguments.join(" ")) + "'!\n\nMake sure the TOPP tools are present in '" + File::getExecutablePath() + "', that you have permission to write to the temporary file path, and that there is space left in the temporary file path.").c_str());
      tool_ready_ = false;
      return false;
    }
    if (!File::exists(ini_file))
    {
      QMessageBox::critical(0, "Error", (String("Could not open '") + ini_file + "'!").c_str());
      tool_ready_ = false;
      return false;
    }

    ParamXMLFile paramFile;
    paramFile.load(String(ini_file).c_str(), tmp_param);
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
    QFile::remove(ini_file);

    setToolTip(param_.getSectionDescription(name_).toQString());

    return changed;
  }

  void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
  {
    editParam();
  }

  void TOPPASToolVertex::editParam()
  {
    QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
    String default_dir = "";

    // use a copy for editing
    Param edit_param(param_);

    QVector<String> hidden_entries;
    // remove entries that are handled by edges already, user should not see them
    QVector<IOInfo> input_infos;
    getInputParameters(input_infos);
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

    QVector<IOInfo> output_infos;
    getOutputParameters(output_infos);
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
    TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, name_, type_, toolTip(), hidden_entries);
    if (dialog.exec())
    {
      // take new values
      param_.update(edit_param);
      reset(true);
      emit parameterChanged(doesParamChangeInvalidate_());
    }

    qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
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

  void TOPPASToolVertex::getInputParameters(QVector<IOInfo>& input_infos) const
  {
    getParameters_(input_infos, true);
  }

  void TOPPASToolVertex::getOutputParameters(QVector<IOInfo>& output_infos) const
  {
    getParameters_(output_infos, false);
  }

  void TOPPASToolVertex::getParameters_(QVector<IOInfo>& io_infos, bool input_params) const
  {
    String search_tag = input_params ? "input file" : "output file";

    io_infos.clear();

    for (Param::ParamIterator it = param_.begin(); it != param_.end(); ++it)
    {
      if (it->tags.count(search_tag))
      {
        StringList valid_types(it->valid_strings);
        for (Size i = 0; i < valid_types.size(); ++i)
        {
          if (!valid_types[i].hasPrefix("*."))
          {
            std::cerr << "Invalid restriction \"" + valid_types[i] + "\"" + " for parameter \"" + it->name + "\"!" << std::endl;
            break;
          }
          valid_types[i] = valid_types[i].suffix(valid_types[i].size() - 2);
        }

        IOInfo io_info;
        io_info.param_name = it.getName();
        io_info.valid_types = valid_types;
        if (it->value.valueType() == DataValue::STRING_LIST)
        {
          io_info.type = IOInfo::IOT_LIST;
        }
        else if (it->value.valueType() == DataValue::STRING_VALUE)
        {
          io_info.type = IOInfo::IOT_FILE;
        }
        else
        {
          std::cerr << "TOPPAS: Unexpected parameter value!" << std::endl;
        }
        io_infos.push_back(io_info);
      }
    }
    // order in param can change --> sort
    qSort(io_infos);
  }

  void TOPPASToolVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
  {
    //painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform);
    QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
    if (isSelected())
    {
      pen.setWidth(2);
      painter->setBrush(brush_color_.darker(130));
      pen.setColor(Qt::darkBlue);
    }
    else
    {
      painter->setBrush(brush_color_);
    }
    painter->setPen(pen);

    QPainterPath path;
    path.addRect(-70.0, -60.0, 140.0, 120.0);
    painter->drawPath(path);

    pen.setColor(pen_color_);
    painter->setPen(pen);

    QString tmp_str = (type_ == "" ? name_ : name_ + " (" + type_ + ")").toQString();
    for (int i = 0; i < 10; ++i)
    {
      QString prev_str = tmp_str;
      tmp_str = toolnameWithWhitespacesForFancyWordWrapping_(painter, tmp_str);
      if (tmp_str == prev_str)
      {
        break;
      }
    }
    QString draw_str = tmp_str;

    QRectF text_boundings = painter->boundingRect(QRectF(-65, -35, 130, 70), Qt::AlignCenter | Qt::TextWordWrap, draw_str);
    painter->drawText(text_boundings, Qt::AlignCenter | Qt::TextWordWrap, draw_str);

    //topo sort number
    qreal x_pos = -64.0;
    qreal y_pos = -41.0;
    painter->drawText(x_pos, y_pos, QString::number(topo_nr_));

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

    // recycling status
    if (this->allow_output_recycling_)
    {
      painter->setPen(Qt::green);
      QSvgRenderer* svg_renderer = new QSvgRenderer(QString(":/Recycling_symbol.svg"), 0);
      svg_renderer->render(painter, QRectF(-7, -52, 14, 14));
    }

    //breakpoint set?
    if (this->breakpoint_set_)
    {
      QSvgRenderer* svg_renderer = new QSvgRenderer(QString(":/stop_sign.svg"), 0);
      painter->setOpacity(0.35);
      svg_renderer->render(painter, QRectF(-60, -60, 120, 120));
    }
  }

  QString TOPPASToolVertex::toolnameWithWhitespacesForFancyWordWrapping_(QPainter* painter, const QString& str)
  {
    qreal max_width = 130;

    QStringList parts = str.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    QStringList new_parts;

    foreach(QString part, parts)
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
          if (QRegExp("[A-Z]").exactMatch(tmp_str.at(i - 1)))
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

  QPainterPath TOPPASToolVertex::shape() const
  {
    QPainterPath shape;
    shape.addRect(-71.0, -61.0, 142.0, 122.0);
    return shape;
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
    if (!isUpstreamFinished()) return;

    if (finished_)
    {
      LOG_ERROR << "This should not happen. Calling an already finished node '" << this->name_ << "' (#" << this->getTopoNr() << ")!" << std::endl;
      throw Exception::IllegalSelfOperation(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());

    QString ini_file = ts->getTempDir()
                       + QDir::separator()
                       + getOutputDir().toQString()
                       + QDir::separator()
                       + name_.toQString();
    if (type_ != "")
      ini_file += "_" + type_.toQString();
    // do not write the ini yet - we might need to alter it

    RoundPackages pkg;
    String error_msg("");
    bool success = buildRoundPackages(pkg, error_msg);
    if (!success)
    {
      std::cerr << "Could not retrieve input files from upstream nodes...\n";
      emit toolFailed(error_msg.toQString());
      return;
    }

    // all inputs are ready --> GO!
    if (!updateCurrentOutputFileNames(pkg, error_msg)) // based on input, we prepare output names
    {
      emit toolFailed(error_msg.toQString());
      return;
    }

    createDirs();

    //emit toolStarted(); //disabled! Every signal emitted here does only mean the process is queued(!) not that its executed right away

    /// update round status
    round_total_ = (int) pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0; // once round_counter_ reaches round_total_, we are done

    QStringList shared_args;
    if (type_ != "")
      shared_args << "-type" << type_.toQString();

    // get *all* input|output file parameters (regardless if edge exists)
    QVector<IOInfo> in_params, out_params;
    getInputParameters(in_params);
    getOutputParameters(out_params);

    bool ini_round_dependent = false; // indicates if we need a new INI file for each round (usually GenericWrapper issue)

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
          std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
          return;
        }

        String param_name = in_params[param_index].param_name;

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        if (param_name.hasPrefix("ETool:"))
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }
        if (!store_to_ini)
          args << "-" + param_name.toQString();

        QStringList file_list = ite->second.filenames;

        if (store_to_ini)
        {
          if (param_tmp.getValue(param_name).valueType() == DataValue::STRING_LIST)
          {
            param_tmp.setValue(param_name, StringListUtils::fromQStringList(file_list));
          }
          else
          {
            if (file_list.size() > 1)
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            }
            param_tmp.setValue(param_name, String(file_list[0]));
          }
        }
        else
        {
          args << file_list;
        }

      }

      // OUTGOING EDGES
      // ;output names are already prepared by 'updateCurrentOutputFileNames()'
      typedef RoundPackage::iterator EdgeIndexIt;
      for (EdgeIndexIt it_edge  = output_files_[round].begin();
           it_edge != output_files_[round].end();
           ++it_edge)
      {
        int param_index = it_edge->first;
        String param_name = out_params[param_index].param_name;

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        if (param_name.hasPrefix("ETool:"))
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }
        if (!store_to_ini)
          args << "-" + param_name.toQString();

        const QStringList& output_files = output_files_[round][param_index].filenames;

        if (store_to_ini)
        {
          if (param_tmp.getValue(param_name).valueType() == DataValue::STRING_LIST)
          {
            param_tmp.setValue(param_name, StringListUtils::fromQStringList(output_files));
          }
          else
          {
            if (output_files.size() > 1) throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            param_tmp.setValue(param_name, String(output_files[0]));
          }
        }
        else
        {
          args << output_files;
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
        LOG_DEBUG << msg_enqueue << std::endl;
        // show sys-call in logWindow of TOPPAS (or console for non-gui)
        if ((int) param_tmp.getValue("debug") > 0)
        {
          ts->logTOPPOutput(msg_enqueue.toQString());
        }
      }
      toolScheduledSlot();
      ts->enqueueProcess(TOPPASScene::TOPPProcess(p, File::findExecutable(name_).toQString(), args, this));
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

    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());

    //** ERROR handling
    if (es != QProcess::NormalExit)
    {
      emit toolCrashed();
    }
    else if (ec != 0)
    {
      emit toolFailed();
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
          LOG_ERROR << "SOMETHING is very fishy. The vertex is already set to finished, yet there was still a thread spawning..." << std::endl;
          throw Exception::IllegalSelfOperation(__FILE__, __LINE__, __PRETTY_FUNCTION__);
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

    //clean up
    QProcess* p = qobject_cast<QProcess*>(QObject::sender());
    if (p)
    {
      delete p;
    }

    ts->processFinished();

    __DEBUG_END_METHOD__
  }

  bool TOPPASToolVertex::renameOutput_()
  {
    // get all output names
    QStringList files = this->getFileNames();

    std::map<String, NameComponent> name_old_to_new;
    Map<String, int> name_new_count, name_new_idx; // count occurrence (for optional counter infix)

    // a first round to find which filenames are not unique (and require augmentation with a counter)

    foreach(QString file, files)
    {
      QFileInfo fi(file);
      String new_suffix = FileTypes::typeToName(FileHandler::getTypeByContent(file));
      String new_prefix = String(fi.path() + "/" + fi.baseName()) + ".";
      NameComponent nc(new_prefix, new_suffix);
      name_old_to_new[file] = nc;
      ++name_new_count[nc.toString()];
    }
    // for all names which occur more than once, introduce a counter  
    foreach(QString file, files)
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
          // rename file and update record
          QFile file(it->second.filenames[fi]);
          String new_filename = name_old_to_new[it->second.filenames[fi]].toString();
          if (File::exists(new_filename))
          {
            bool success = File::remove(new_filename);
            if (!success)
            {
              std::cerr << "Could not remove '" << new_filename << "'.\n";
              return false;
            }
          }
          bool success = file.rename(new_filename.toQString());
          if (!success)
          {
            std::cerr << "Could not rename " << String(it->second.filenames[fi]) << " to " << new_filename << "\n";
            return false;
          }
          it->second.filenames[fi] = new_filename.toQString();
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
    if (pkg.size() < 1)
    {
      error_msg = "Less than one round received from upstream tools. Something is fishy!\n";
      std::cerr << error_msg;
      return false;
    }

    // look for the input with the most files in round 0 (as this is the maximal number of output files we can produce)
    // we assume the number of files is equal in all rounds...
    // However, we delay using nodes which use 'recycling' of input, as the names will always be the same.
    //          Only if a recycling node gives the most input files we use its names
    int max_size_index = -1;
    int max_size = -1;
    for (int use_recycling = 0; use_recycling < 2; ++use_recycling)
    {
      for (RoundPackageConstIt it  = pkg[0].begin();
           it != pkg[0].end();
           ++it)
      {
        if (use_recycling == 0 && (it->second.edge->getSourceVertex()->isRecyclingEnabled()))
        { // first test all input nodes with disabled recycling
          continue;
        }
        //std::cerr << "Edge: " << (it->second.edge->toString()) << "\n";
        if ((it->second.filenames.size() > max_size) ||   // either just larger 
            // ... or it's from '-in' (which we prefer as naming source).. only for non-recycling -in though
            ((it->second.filenames.size () == max_size) && (it->second.edge->getTargetInParamName() == "in") && (use_recycling == 0)))
        {
          max_size_index = it->first;
          max_size       = it->second.filenames.size();
        }
      }

      if (max_size_index == -1)
      {
        error_msg = "Did not find upstream nodes with un-recycled names. Something is fishy!\n";
        LOG_ERROR << error_msg;
        return false;
      }
    }

    // now we construct output filenames for this node
    // use names from the selected upstream vertex (hoping that this is the maximal number of files we are going to produce)
    std::vector<QStringList> per_round_basenames;
    for (Size i = 0; i < pkg.size(); ++i)
    {
      per_round_basenames.push_back(pkg[i].find(max_size_index)->second.filenames);
      //String s = String(pkg[i].find(max_size_index)->second.filenames.join(" + "));
    }

    // maybe we find something more unique, e.g. last base directory if all filenames are equal
    smartFileNames_(per_round_basenames);

    QRegExp rx_tmp_suffix("_tmp\\d+$"); // regex to remove "_tmp<number>" if its a suffix

    // clear output file list
    output_files_.clear();
    output_files_.resize(pkg.size()); // #rounds

    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    QVector<IOInfo> out_params;
    getOutputParameters(out_params);
    // output names for each outgoing edge
    for (int i = 0; i < out_params.size(); ++i)
    {
      // search for an out edge for this parameter (not required to exist)
      bool found(false);
      int param_index;
      TOPPASEdge* edge_out;
      for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
      {
        param_index = (*it)->getSourceOutParam();
        if (i == param_index) // corresponding out edge found
        {
          edge_out = *it;
          found = true;
          break;
        }
      }
      if (!found)
      {
        continue;
      }

      // create common path of output files
      QString path = ts->getTempDir()
                     + QDir::separator()
                     + getOutputDir().toQString()
                     + QDir::separator()
                     + out_params[param_index].param_name.remove(':').toQString().left(50) // max 50 chars per subdir
                     + QDir::separator();
      if (path.length() > 150)
      {
        LOG_WARN << "Warning: the temporary path '" << String(path) << "' used in TOPPAS has many characters.\n"
                 << "         TOPPAS might not be able to write files properly.\n";
      }

      VertexRoundPackage vrp;
      vrp.edge = edge_out;
      for (Size r = 0; r < per_round_basenames.size(); ++r)
      {
        // store edge for this param for all rounds
        output_files_[r][param_index] = vrp; // index by index of source-out param

        // check if tool consumes list and outputs single file (such as IDMerger or FileMerger)
        if (per_round_basenames[r].size() > 1 && out_params[param_index].type == IOInfo::IOT_FILE)
        {
          QString fn = path + QString(QFileInfo(per_round_basenames[r].first()).fileName()
                            + "_to_"
                            + QFileInfo(per_round_basenames[r].last()).fileName()
                            + "_merged");
          fn = fn.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
          fn += "_tmp" + QString::number(uid_++);
          fn = QDir::toNativeSeparators(fn);
          output_files_[r][param_index].filenames.push_back(fn);
        }
        else // each input file will have a corresponding output file
        {
          foreach(const QString &input_file, per_round_basenames[r])
          {
            QString fn = path + QFileInfo(input_file).fileName(); // out_path + filename
            int tmp_index = rx_tmp_suffix.indexIn(fn);
            if (tmp_index != -1)
            {
              fn = fn.left(tmp_index);
            }
            fn = fn.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
            fn += "_tmp" +  QString("%1").arg(uid_++, 3, 10, QChar('0')); // pad with zeros '00X' etc...
            fn = QDir::toNativeSeparators(fn);
            output_files_[r][param_index].filenames.push_back(fn);
          }
        }
      }
    }
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
        if (p.isEmpty()) continue;
        //std::cout << "PATH: " << p << "\n";
        String tmp = String(p).suffix(String(QString(QDir::separator()))[0]);
        //std::cout << "INTER: " << tmp << "\n";
        if (tmp.size() <= 2 || tmp.has(':')) continue; // too small to be reliable; might even be 'c:'
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

  void TOPPASToolVertex::openContainingFolder()
  {
    QString path = getFullOutputDirectory().toQString();
    GUIHelpers::openFolder(path);
  }

  String TOPPASToolVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    return QDir::toNativeSeparators(ts->getTempDir() + QDir::separator() + getOutputDir().toQString());
  }

  String TOPPASToolVertex::getOutputDir() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    String workflow_dir = File::removeExtension(File::basename(ts->getSaveFileName()));
    if (workflow_dir == "")
    {
      workflow_dir = "Untitled_workflow";
    }
    String dir = String("TOPPAS_tmp") +
                 String(QDir::separator()) +
                 workflow_dir +
                 String(QDir::separator()) +
                 get3CharsNumber_(topo_nr_) + "_" + getName();
    if (getType() != "")
    {
      dir += "_" + getType();
    }

    return dir;
  }

  void TOPPASToolVertex::createDirs()
  {
    QDir dir;
    bool ok = dir.mkpath(getFullOutputDirectory().toQString());

    if (!ok)
    {
      std::cerr << "TOPPAS: Could not create path " << getFullOutputDirectory() << std::endl;
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
          std::cerr << "TOPPAS: Could not create path " << String(sdir) << std::endl;
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
      // reset UID for tmp files
      uid_ = 1;
    }

    TOPPASVertex::reset(reset_all_files);

    __DEBUG_END_METHOD__
  }

  bool TOPPASToolVertex::refreshParameters()
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    QString old_ini_file = ts->getTempDir() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
    if (type_ != "")
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
    save_param.setValue(name_ + ":1:toppas_dummy", DataValue("blub"));
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
