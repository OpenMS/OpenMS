// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>

#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <QPainter>
#include <QPainterPath>
#include <QtWidgets/QMessageBox>

#include <QApplication>

namespace OpenMS
{

  TOPPASEdge::TOPPASEdge() :
    QObject(),
    QGraphicsItem(),
    from_(nullptr),
    to_(nullptr),
    hover_pos_(),
    color_(),
    source_out_param_(-1),
    target_in_param_(-1)
  {
    setFlag(QGraphicsItem::ItemIsSelectable, true);
  }

  TOPPASEdge::TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos) :
    QObject(),
    QGraphicsItem(),
    from_(from),
    to_(nullptr),
    hover_pos_(hover_pos),
    color_(),
    source_out_param_(-1),
    target_in_param_(-1)
  {
    setFlag(QGraphicsItem::ItemIsSelectable, true);
  }

  TOPPASEdge::TOPPASEdge(const TOPPASEdge& rhs) :
    QObject(),
    QGraphicsItem(),
    from_(rhs.from_),
    to_(rhs.to_),
    hover_pos_(rhs.hover_pos_),
    color_(rhs.color_),
    source_out_param_(rhs.source_out_param_),
    target_in_param_(rhs.target_in_param_)
  {
    setFlag(QGraphicsItem::ItemIsSelectable, true);
  }

  String TOPPASEdge::toString()
  {
    String s = String("Edge: ") + getSourceOutParamName() + " target-in: " + getTargetInParamName() + "\n";
    return s;
  }
  TOPPASEdge& TOPPASEdge::operator=(const TOPPASEdge& rhs)
  {
    from_ = rhs.from_;
    to_ = rhs.to_;
    hover_pos_ = rhs.hover_pos_;
    color_ = rhs.color_;
    source_out_param_ = rhs.source_out_param_;
    target_in_param_ = rhs.target_in_param_;

    setFlag(QGraphicsItem::ItemIsSelectable, true);

    return *this;
  }

  TOPPASEdge::~TOPPASEdge()
  {
    // notify our childs that we are dying ;)
    emit somethingHasChanged();

    if (from_)
    {
      from_->removeOutEdge(this);
      disconnect(from_, SIGNAL(somethingHasChanged()), this, SLOT(sourceHasChanged()));
    }
    if (to_)
    {
      to_->removeInEdge(this);
      disconnect(this, SIGNAL(somethingHasChanged()), to_, SLOT(inEdgeHasChanged()));
    }
  }

  QRectF TOPPASEdge::boundingRect() const
  {
    qreal min_x = startPos().x() < endPos().x() ? startPos().x() : endPos().x();
    qreal min_y = startPos().y() < endPos().y() ? startPos().y() : endPos().y();
    qreal max_x = startPos().x() > endPos().x() ? startPos().x() : endPos().x();
    qreal max_y = startPos().y() > endPos().y() ? startPos().y() : endPos().y();

    return QRectF(QPointF(min_x - 11.0, min_y - 11.0), QPointF(max_x + 11.0, max_y + 11.0));
  }

  QPainterPath TOPPASEdge::shape() const
  {
    QPainterPath shape_1;
    shape_1.moveTo(startPos() + QPointF(-10, -10));
    shape_1.lineTo(endPos() + QPointF(-10, -10));
    shape_1.lineTo(endPos() + QPointF(10, 10));
    shape_1.lineTo(startPos() + QPointF(10, 10));
    shape_1.closeSubpath();

    QPainterPath shape_2;
    shape_2.moveTo(startPos() + QPointF(-10, 10));
    shape_2.lineTo(endPos() + QPointF(-10, 10));
    shape_2.lineTo(endPos() + QPointF(10, -10));
    shape_2.lineTo(startPos() + QPointF(10, -10));
    shape_2.closeSubpath();

    return shape_1.united(shape_2);
  }

  void TOPPASEdge::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
  {
    painter->setBrush(Qt::white);

    QPen pen(color_);
    if (isSelected())
    {
      pen.setWidth(3);
    }
    else
    {
      pen.setWidth(2);
    }

    TOPPASToolVertex* ttv_source = qobject_cast<TOPPASToolVertex*>(this->getSourceVertex());
    // when copying parameters (using CTRL); only for incomplete edges drawn from tool nodes
    if ((QGuiApplication::keyboardModifiers() & Qt::ControlModifier) && !this->to_ && ttv_source)
    {
      pen.setColor(Qt::darkMagenta);
      pen.setWidth(1);
    }

    painter->setPen(pen);

    // angle of line
    qreal angle = -QLineF(endPos(), startPos()).angle() + 180; // negate since angle() reports counter-clockwise; +180 since painter.rotate() is more intuitive then

    // draw the actual line
    QPainterPath path_line(startPos());
    path_line.lineTo(endPos());
    painter->drawPath(path_line);

    // print names
    qreal text_angle = angle;
    bool invert_text_direction = endPos().x() < startPos().x();
    if (invert_text_direction)
    {
      text_angle += 180;
    }
    QPainterPath path_line_short(borderPoint_(false));
    path_line_short.lineTo(endPos());

    // y offset for printing text; we add -1 to make "_" in param names visible (e.g. in "out_annotation")
    //                             otherwise they are too close to the edge itself
    int y_text = -pen.width() - 1;

    // source name
    QString str = getSourceOutParamName();
    if (!str.isEmpty())
    {
      painter->save();     // hard to avoid multiple calls to save() and restore() since we first translate and then rotate
      QPointF point = path_line_short.pointAtPercent(0.05);
      painter->translate(point);
      painter->rotate(text_angle);
      if (invert_text_direction)
      {
        QFontMetrics fm(painter->fontMetrics());
        int text_width = fm.horizontalAdvance(str);
        painter->drawText(QPoint(-text_width, y_text), str);
      }
      else
      {
        painter->drawText(QPoint(0, y_text), str);
      }
      painter->restore();
    }
    // target name
    const qreal arrow_width = 10; // required multiple times, so defined here to avoid inconsistencies
    str = getTargetInParamName();
    if (!str.isEmpty())
    {
      painter->save();
      // qreal pc = path_line_short.percentAtLength(10);
      QPointF point = path_line_short.pointAtPercent(0.95);
      painter->translate(point);
      painter->rotate(text_angle);
      QFontMetrics fm(painter->fontMetrics());
      int text_width = fm.horizontalAdvance(str);
      int text_height = fm.height();  // shift text below the edge by its own height
      if (invert_text_direction)
      {
        painter->drawText(QPoint(arrow_width, text_height + y_text), str);
      }
      else
      {
        painter->drawText(QPoint(-text_width - arrow_width, text_height + y_text), str);
      }
      painter->restore();
    }

    // draw arrow head
    painter->save();
    painter->translate(endPos());
    painter->rotate(angle);
    QPainterPath path;
    path.moveTo(QPointF(0, 0));
    path.lineTo(QPointF(-arrow_width, 4));
    path.lineTo(QPointF(-arrow_width, -4));
    path.closeSubpath();
    painter->drawPath(path);
    painter->restore();
  }

  QPointF TOPPASEdge::startPos() const
  {
    if (from_)
    {
      return mapFromScene(from_->scenePos());
    }
    else
    {
      return QPointF();
    }
  }

  QPointF TOPPASEdge::endPos() const
  {
    if (!to_)
    {
      // we do not have a target vertex yet
      return (hover_pos_);
    }
    else
    {
      // we have a target node --> line should end at its border
      return (borderPoint_());
    }
  }



  QPointF TOPPASEdge::borderPoint_(bool atTargetVertex) const
  {
    if (!to_ || !from_)
    {
      return QPointF(); // both ends need to be fixed; otherwise we have no input/output slots assigned anyways
    }
    const TOPPASVertex* to = (atTargetVertex ? to_ : from_);
    const TOPPASVertex* from = (!atTargetVertex ? to_ : from_);

    const QPointF fromP = mapFromScene(from->scenePos());
    QRectF target_boundings = mapFromItem(to, to->shape()).boundingRect();

    return GUIHelpers::intersectionPoint(target_boundings, fromP);
  }

  void TOPPASEdge::setHoverPos(const QPointF& pos)
  {
    prepareResize();
    hover_pos_ = pos;
    update();
  }

  void TOPPASEdge::setTargetVertex(TOPPASVertex* tv)
  {
    to_ = tv;
  }

  void TOPPASEdge::setSourceVertex(TOPPASVertex* tv)
  {
    from_ = tv;
  }

  TOPPASVertex* TOPPASEdge::getSourceVertex()
  {
    return from_;
  }

  TOPPASVertex* TOPPASEdge::getTargetVertex()
  {
    return to_;
  }

  void TOPPASEdge::prepareResize()
  {
    prepareGeometryChange();
  }

  void TOPPASEdge::setColor(const QColor& color)
  {
    color_ = color;
  }

  void TOPPASEdge::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
  {
    showIOMappingDialog();
  }

  void TOPPASEdge::showIOMappingDialog()
  {
    TOPPASIOMappingDialog dialog(this);
    if (dialog.exec())
    {
      emit somethingHasChanged();
    }
  }

  TOPPASEdge::EdgeStatus TOPPASEdge::getToolToolStatus_(TOPPASToolVertex* source_tool, int source_param_index, TOPPASToolVertex* target_tool, int target_param_index)
  {
    if (source_param_index < 0)
    {
      return ES_NO_SOURCE_PARAM;
    }
    if (target_param_index < 0)
    {
      return ES_NO_TARGET_PARAM;
    }

    QVector<TOPPASToolVertex::IOInfo> source_output_files = source_tool->getOutputParameters();
    if (source_param_index >= source_output_files.size())
    {
      return ES_TOOL_API_CHANGED;
    }

    QVector<TOPPASToolVertex::IOInfo> target_input_files = target_tool->getInputParameters();
    if (target_param_index >= target_input_files.size())
    {
      return ES_TOOL_API_CHANGED;
    }

    const TOPPASToolVertex::IOInfo& source_param = source_output_files[source_param_index];
    StringList source_param_types = source_param.valid_types;

    const TOPPASToolVertex::IOInfo& target_param = target_input_files[target_param_index];
    StringList target_param_types = target_param.valid_types;

    if (source_param_types.empty() || target_param_types.empty())
    {
      // no type specified --> allow edge
      return ES_VALID;
    }
    else
    {
      // check file type compatibility
      bool found_match = false;
      for (StringList::iterator s_it = source_param_types.begin(); s_it != source_param_types.end(); ++s_it)
      {
        String ext_1 = *s_it;
        ext_1.toLower();
        found_match = false;
        for (StringList::iterator t_it = target_param_types.begin(); t_it != target_param_types.end(); ++t_it)
        {
          String ext_2 = *t_it;
          ext_2.toLower();
          if (ext_1 == ext_2)
          {
            found_match = true;
            break;
          }
        }
        if (found_match)
        {
          break;
        }
      }

      if (!found_match)
      {
        return ES_FILE_EXT_MISMATCH;
      }

      return ES_VALID;
    }
  }

  TOPPASEdge::EdgeStatus TOPPASEdge::getListToolStatus_(TOPPASInputFileListVertex* source_input_list, TOPPASToolVertex* target_tool, int target_param_index)
  {
    QVector<TOPPASToolVertex::IOInfo> target_input_files = target_tool->getInputParameters();
    if (target_param_index >= target_input_files.size())
    {
      return ES_TOOL_API_CHANGED;
    }

    const QStringList& file_names = source_input_list->getFileNames();
    if (file_names.empty())
    {
      // file names are not specified yet
      return ES_NOT_READY_YET;
    }

    if (target_param_index < 0)
    {
      return ES_NO_TARGET_PARAM;
    }

    const TOPPASToolVertex::IOInfo& target_param = target_input_files[target_param_index];
    StringList target_param_types = target_param.valid_types;

    if (target_param_types.empty())
    {
      // no file types specified --> allow
      return ES_VALID;
    }

    // check file type compatibility
    foreach(const QString& q_file_name, file_names)
    {
      bool type_mismatch = true;
      const String& file_name = String(q_file_name);
      String::SizeType extension_start_index = file_name.rfind(".");
      if (extension_start_index != String::npos)
      {
        String extension = file_name.substr(extension_start_index + 1);
        extension.toLower();
        for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
        {
          String other_ext = *it;
          other_ext.toLower();
          if (extension == other_ext || extension == "gz" || extension == "bz2")
          {
            type_mismatch = false;
            break;
          }
        }
      }

      if (type_mismatch)
      {
        return ES_FILE_EXT_MISMATCH;
      }
    }

    // all file types ok
    return ES_VALID;
  }

  TOPPASEdge::EdgeStatus TOPPASEdge::getEdgeStatus()
  {
    TOPPASVertex* source = getSourceVertex();
    TOPPASVertex* target = getTargetVertex();
    TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
    TOPPASMergerVertex* target_merger = qobject_cast<TOPPASMergerVertex*>(target);
    TOPPASSplitterVertex* source_splitter = qobject_cast<TOPPASSplitterVertex*>(source);
    TOPPASSplitterVertex* target_splitter = qobject_cast<TOPPASSplitterVertex*>(target);
    TOPPASInputFileListVertex* source_input_list = qobject_cast<TOPPASInputFileListVertex*>(source);
    TOPPASOutputVertex* target_output = qobject_cast<TOPPASOutputVertex*>(target);
    TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
    TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);

    if (source_tool && source_out_param_ < 0)
    {
      return ES_NO_SOURCE_PARAM;
    }

    if (target_tool && target_in_param_ < 0)
    {
      return ES_NO_TARGET_PARAM;
    }

    if (source_tool)
    {
      if (target_output) // edges to output vertices are always valid
      {
        return ES_VALID;
      }
      if (target_tool)
      {
        return getToolToolStatus_(source_tool, source_out_param_, target_tool, target_in_param_);
      }
    }

    if (source_input_list && target_tool)
    {
      return getListToolStatus_(source_input_list, target_tool, target_in_param_);
    }

    if (source_merger || source_splitter)
    {
      // check compatibility of source with all target nodes of merger (or splitter)
      for (TOPPASVertex::ConstEdgeIterator e_it = source->inEdgesBegin(); e_it != source->inEdgesEnd(); ++e_it)
      {
        TOPPASEdge* merger_in_edge = *e_it;
        TOPPASToolVertex* merger_in_tool = qobject_cast<TOPPASToolVertex*>(merger_in_edge->getSourceVertex());
        TOPPASInputFileListVertex * merger_in_list = qobject_cast<TOPPASInputFileListVertex*>(merger_in_edge->getSourceVertex());

        if (merger_in_tool && target_tool)
        {
          EdgeStatus es = getToolToolStatus_(merger_in_tool, merger_in_edge->getSourceOutParam(), target_tool, target_in_param_);
          if (es != ES_VALID)
          {
            return es;
          }
        }
        else if (merger_in_list)
        {
          if (target_output)
          {
            // [input] -> [merger] -> [output]) makes no sense
            return ES_MERGER_WITHOUT_TOOL;
          }
          else if (target_tool)
          {
            EdgeStatus es = getListToolStatus_(merger_in_list, target_tool, target_in_param_);
            if (es != ES_VALID)
            {
              return es;
            }
          }
        }
      }
      // no incompatible merger target found
      return ES_VALID;
    }

    if (target_merger || target_splitter)
    {
      // check compatibility of source with all target nodes of merger (or splitter)
      for (TOPPASVertex::ConstEdgeIterator e_it = target->outEdgesBegin(); e_it != target->outEdgesEnd(); ++e_it)
      {
        TOPPASEdge* merger_out_edge = *e_it;
        TOPPASToolVertex* merger_out_tool = qobject_cast<TOPPASToolVertex*>(merger_out_edge->getTargetVertex());
        TOPPASOutputFileListVertex* merger_out_output = qobject_cast<TOPPASOutputFileListVertex*>(merger_out_edge->getTargetVertex());

        if (source_tool)
        {
          if (merger_out_tool)
          {
            EdgeStatus es = getToolToolStatus_(source_tool, source_out_param_, merger_out_tool, merger_out_edge->getTargetInParam());
            if (es != ES_VALID)
            {
              return es;
            }
          }
        }
        else if (source_input_list)
        {
          if (merger_out_tool)
          {
            EdgeStatus es = getListToolStatus_(source_input_list, merger_out_tool, merger_out_edge->getTargetInParam());
            if (es != ES_VALID)
            {
              return es;
            }
          }
          else if (merger_out_output)
          {
            // [input] -> [merger] -> [output]) makes no sense
            return ES_MERGER_WITHOUT_TOOL;
          }
        }
      }
      // no incompatible merger target found
      return ES_VALID;
    }

    return ES_UNKNOWN;
  }

  void TOPPASEdge::setSourceOutParam(int out)
  {
    source_out_param_ = out;
  }

  int TOPPASEdge::getSourceOutParam() const
  {
    return source_out_param_;
  }

  void TOPPASEdge::setTargetInParam(int in)
  {
    target_in_param_ = in;
  }

  int TOPPASEdge::getTargetInParam() const
  {
    return target_in_param_;
  }

  QString TOPPASEdge::getTargetInParamName()
  {
    const EdgeStatus& es = getEdgeStatus();
    if (es != ES_TOOL_API_CHANGED)
    {
      TOPPASVertex* target_o = getTargetVertex();
      const TOPPASToolVertex* target = qobject_cast<TOPPASToolVertex*>(target_o);
      if (target && target_in_param_>=0)
      {
         QVector<TOPPASToolVertex::IOInfo> docks = target->getInputParameters();
         const TOPPASToolVertex::IOInfo& param = docks[this->target_in_param_]; 
         return param.param_name.toQString();
      }
    }
    return "";
  }


  QString TOPPASEdge::getSourceOutParamName()
  {
    const EdgeStatus& es = getEdgeStatus();
    if (es != ES_TOOL_API_CHANGED)
    {
      const auto* source_vertex = getSourceVertex();
      const auto* source_tool = qobject_cast<const TOPPASToolVertex*>(source_vertex);
      if (source_tool && source_out_param_>=0)
      {
        return source_tool->getOutputParameters()[this->source_out_param_].param_name.toQString();
      }
    }
    return "";
  }

  void TOPPASEdge::updateColor()
  {
    const EdgeStatus& es = getEdgeStatus();
    if (es == ES_VALID)
    {
      setColor(Qt::darkGreen);
    }
    else if (es == ES_NOT_READY_YET)
    {
      setColor(QColor(255, 165, 0));
    }
    else
    {
      setColor(Qt::red);
    }
    update(boundingRect());
  }

  void TOPPASEdge::sourceHasChanged()
  {
    emit somethingHasChanged();
  }

  void TOPPASEdge::emitChanged()
  {
    emit somethingHasChanged();
  }

  void TOPPASEdge::contextMenuEvent(QGraphicsSceneContextMenuEvent* event)
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    ts->unselectAll();
    setSelected(true);

    QMenu menu;
    menu.addAction("Edit I/O mapping");
    menu.addAction("Remove");

    QAction* selected_action = menu.exec(event->screenPos());
    if (selected_action)
    {
      QString text = selected_action->text();
      if (text == "Edit I/O mapping")
      {
        TOPPASIOMappingDialog dialog(this);
        if (dialog.exec())
        {
          emit somethingHasChanged();
        }
      }
      else if (text == "Remove")
      {
        ts->removeSelected();
      }
      event->accept();
    }
    else
    {
      event->ignore();
    }
  }

} //namespace
