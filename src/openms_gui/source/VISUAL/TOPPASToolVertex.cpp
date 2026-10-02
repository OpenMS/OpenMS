// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <QSvgRenderer>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QProcess>
#include <QtCore/QRegularExpression>
#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QMessageBox>
#include <map>

namespace OpenMS
{
  TOPPASToolVertex::TOPPASToolVertex()
    : TOPPASToolVertex("", "")
  {
  }

  TOPPASToolVertex::TOPPASToolVertex(const std::string& name, const std::string& type) :
    name_(name),
    type_(type)
  {
    brush_color_ = brush_color_.lighter(130); // make TOPP tools more white compared to all other nodes
    initParam_();
    connect(this, &TOPPASToolVertex::toolStarted, this, &TOPPASToolVertex::toolStartedSlot);
    connect(this, &TOPPASToolVertex::toolFinished, this, &TOPPASToolVertex::toolFinishedSlot);
    connect(this, &TOPPASToolVertex::toolFailed, this, &TOPPASToolVertex::toolFailedSlot);
    connect(this, &TOPPASToolVertex::toolCrashed, this, &TOPPASToolVertex::toolCrashedSlot);
  }

  TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs):
      TOPPASVertex(rhs),
      name_(rhs.name_),
      type_(rhs.type_),
      param_(rhs.param_),
      status_(TOOL_READY),
      tool_ready_(rhs.tool_ready_)

  {
  }

  TOPPASToolVertex& TOPPASToolVertex::operator=(const TOPPASToolVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);

    param_ = rhs.param_;
    name_ = rhs.name_;
    type_ = rhs.type_;
    finished_ = false;
    status_ = TOOL_READY;
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
    QString ini_file = toQString(TempFiles::getTemporaryFile());
    QString program = toQString(File::findSiblingTOPPExecutable(name_));
    QStringList arguments;
    arguments << "-write_ini" << ini_file;

    if (!type_.empty())
    {
      arguments << "-type";
      arguments << toQString(type_);
    }
    // allow for update using old parameters
    if (old_ini_file != "")
    {
      if (!File::exists(fromQString(old_ini_file)))
      {
        std::string msg =std::string("Could not open old INI file '") + fromQString(old_ini_file) + "'! File does not exist!";
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
      std::string msg =std::string("Error! Call to '") + fromQString(program) + "' '" + fromQString(arguments.join("' '")) +
          " returned with exit code (" + StringUtils::toStr(p.exitCode()) + "), exit status (" + StringUtils::toStr(p.exitStatus()) + ")." +
          "\noutput:\n" + fromQString(QString(p.readAll())) +
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
    if (!File::exists(fromQString(ini_file)))
    { // it would be weird to get here, since the TOPP tool ran successfully above, so INI file should exist, but nevertheless:
      std::string msg =std::string("Could not open '") + fromQString(ini_file) + "'! It does not exist!";
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
    ParamXMLFile().load(fromQString(ini_file).c_str(), tmp_param);
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
    setToolTip(toQString(std::string(param_.getSectionDescription(name_))));

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

    QVector<std::string> hidden_entries;
    // remove entries that are handled by edges already, user should not see them
    QVector<IOInfo> input_infos = getInputParameters();
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      int index = (*it)->getTargetInParam();
      if (index < 0)
      {
        continue;
      }

      const std::string& name = input_infos[index].param_name;
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

      const std::string& name = output_infos[index].param_name;
      if (edit_param.exists(name))
      {
        hidden_entries.push_back(name);
      }
    }

    // remove entries explained by edges
    for (const std::string &name : hidden_entries)
    {
      edit_param.remove(name);
    }

    // edit_param no longer contains tool description, take it from the node tooltip
    QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
    std::string default_dir;
    TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, name_, type_, fromQString(toolTip()), hidden_entries);
    if (dialog.exec())
    {
      // take new values
      getScene_()->abortPipeline();
      param_.update(edit_param);
      emit parameterChanged(true);
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
    auto add_params = [&](const std::string& search_tag) {
      for (Param::ParamIterator it = param_.begin(); it != param_.end(); ++it)
      {
        if (! it->tags.count(search_tag)) continue; // skip irrelevant parameters

        StringList valid_types(ListUtils::toStringList<std::string>(it->valid_strings));
        for (Size i = 0; i < valid_types.size(); ++i)
        {
          if (! StringUtils::hasPrefix(valid_types[i], "*."))
          {
            std::cerr << "Invalid restriction \"" + valid_types[i] + "\"" + " for parameter \"" + it->name + "\"!" << std::endl;
            break;
          }
          valid_types[i] = StringUtils::suffix(valid_types[i], valid_types[i].size() - 2);
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

    QString draw_str = toQString(type_.empty() ? name_ : name_ + " (" + type_ + ")");
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

  std::string TOPPASToolVertex::getName() const
  {
    return name_;
  }

  const std::string& TOPPASToolVertex::getType() const
  {
    return type_;
  }

  void TOPPASToolVertex::run()
  {
    if (auto* pipeline = getScene_()) { pipeline->resumePipeline(this); }
  }

  void TOPPASToolVertex::emitToolStarted()
  {
    emit toolStarted();
  }


  const Param& TOPPASToolVertex::getParam()
  {
    return param_;
  }

  void TOPPASToolVertex::setParam(const Param& param)
  {
    // Saved workflows may contain only overrides. Keep the current executable's
    // port metadata and defaults, with the same strict checks as PipelineTool.
    Param updated(param_);
    if (! updated.update(param, false, false, true, true, OPENMS_LOG_WARN))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Saved parameters for '" + name_ + "' are incompatible with its current schema.");
    }
    param_ = std::move(updated);
  }

  TOPPASToolVertex::TOOLSTATUS TOPPASToolVertex::getStatus() const
  { return status_; }


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
    QString path = toQString(getFullOutputDirectory());
    GUIHelpers::openFolder(path);
  }

  std::string TOPPASToolVertex::getFullOutputDirectory() const
  {
    const auto files = getFileNames();
    if (! files.empty()) { return File::path(fromQString(files.front())); }
    TOPPASScene* ts = getScene_();
    return fromQString(QDir::toNativeSeparators(ts->getTempDir() + QDir::separator() + toQString(getOutputDir())));
  }

  std::string TOPPASToolVertex::getOutputDir() const
  {
    TOPPASScene* ts = getScene_();
    std::string workflow_dir = File::stemName(ts->getSaveFileName());
    if (workflow_dir.empty())
    {
      workflow_dir = "Untitled_workflow";
    }
    std::string dir = workflow_dir +
                 fromQString(QString(QDir::separator())) +
                 get3CharsNumber_(topo_nr_) + "_" + getName();
    if (!getType().empty())
    {
      dir += "_" + getType();
    }

    return dir;
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

    TOPPASVertex::reset(reset_all_files);

    __DEBUG_END_METHOD__
  }

  bool TOPPASToolVertex::refreshParameters()
  {
    TOPPASScene* ts = getScene_();
    QString old_ini_file = ts->getTempDir() + QDir::separator() + "TOPPAS_" + toQString(name_) + "_";
    if (!type_.empty())
    {
      old_ini_file += toQString(type_) + "_";
    }
    old_ini_file += toQString(File::getUniqueName()) + "_tmp_OLD.ini";
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
    paramFile.store(fromQString(ini_file), save_param);
  }

  void TOPPASToolVertex::toggleBreakpoint()
  {
    breakpoint_set_ = !breakpoint_set_;
  }

}
