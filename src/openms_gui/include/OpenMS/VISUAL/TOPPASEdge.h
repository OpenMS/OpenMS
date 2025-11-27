// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QGraphicsItem>

namespace OpenMS
{
  class TOPPASVertex;
  class TOPPASToolVertex;
  class TOPPASInputFileListVertex;

  class String;

  /**
      @brief An edge representing a data flow in TOPPAS

      Like all TOPPASVertex classes, TOPPASEdge is a subclass of QGraphicsItem and thus implements
      methods to draw itself and to react on incoming events such as mouse clicks. It holds
      the data needed to represent an edge between two vertices of a TOPPAS workflow.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASEdge :
    public QObject,
    public QGraphicsItem
  {
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
public:

    /// The status of this edge
    enum EdgeStatus
    {
      ES_VALID,
      ES_NO_TARGET_PARAM,
      ES_NO_SOURCE_PARAM,
      ES_FILE_EXT_MISMATCH,
      ES_MERGER_EXT_MISMATCH,
      ES_MERGER_WITHOUT_TOOL,
      ES_NOT_READY_YET,       // no input files given. We cannot know if the types will match.
      ES_TOOL_API_CHANGED,
      ES_UNKNOWN
    };

    /// Standard constructor
    TOPPASEdge();
    /// Constructor
    TOPPASEdge(TOPPASVertex * from, const QPointF & hover_pos);
    /// Copy constructor
    TOPPASEdge(const TOPPASEdge & rhs);
    /// Destructor
    ~TOPPASEdge() override;
    /// Assignment operator
    TOPPASEdge & operator=(const TOPPASEdge & rhs);

    /// for debug output
    String toString();

    /// Returns the bounding rectangle of this item
    QRectF boundingRect() const override;
    /// Returns a more precise shape
    QPainterPath shape() const override;
    /// Paints the item
    void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override;
    /// Returns the start position of this edge
    QPointF startPos() const;
    /// Returns the end position of this edge
    QPointF endPos() const;
    /// Sets the position of the hovering end while edge is being created
    void setHoverPos(const QPointF & pos);
    /// Sets the source vertex of this edge
    void setSourceVertex(TOPPASVertex * tv);
    /// Sets the target vertex of this edge
    void setTargetVertex(TOPPASVertex * tv);
    /// Returns the source vertex
    TOPPASVertex * getSourceVertex();
    /// Returns the target vertex
    TOPPASVertex * getTargetVertex();
    /// Call this before changing the item geometry
    void prepareResize();
    /// Sets the color
    void setColor(const QColor & color);
    /// Returns the status of this edge
    EdgeStatus getEdgeStatus();
    /// Sets the source output parameter index
    void setSourceOutParam(int out);
    /// Returns the source output parameter index
    int getSourceOutParam() const;
    /// Returns the source output parameter name
    QString getSourceOutParamName();
    /// Sets the target input parameter index
    void setTargetInParam(int in);
    /// Returns the target input parameter index
    int getTargetInParam() const;
    /// Returns the target input parameter index
    QString getTargetInParamName();
    /// Updates the edge color
    void updateColor();
    /// Emits the somethingHasChanged() signal
    void emitChanged();
    /// Shows the I/O mapping dialog
    void showIOMappingDialog();

public slots:

    /// Called by the source vertex when it has changed
    void sourceHasChanged();

signals:

    /// Emitted when something has changed
    void somethingHasChanged();

protected:

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent * e) override;
    void contextMenuEvent(QGraphicsSceneContextMenuEvent * event) override;
    //@}

    ///@name helper methods of getEdgeStatus()
    //@{
    EdgeStatus getToolToolStatus_(TOPPASToolVertex * source, int source_param_index, TOPPASToolVertex * target, int target_param_index);
    EdgeStatus getListToolStatus_(TOPPASInputFileListVertex * source, TOPPASToolVertex * target, int target_param_index);
    //@}

    /// point where the current edge touches the source or target (default) vertex
    QPointF borderPoint_(bool atTargetVertex = true) const;

    /// Pointer to the source of this edge
    TOPPASVertex * from_;
    /// Pointer to the target of this edge
    TOPPASVertex * to_;
    /// Position of hovering end while edge is being created
    QPointF hover_pos_;
    /// The color
    QColor color_;
    /// The source output parameter index
    int source_out_param_;
    /// The target input parameter index
    int target_in_param_;
  };
}

