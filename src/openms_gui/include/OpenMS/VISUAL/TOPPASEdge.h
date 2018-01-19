// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASEDGE_H
#define OPENMS_VISUAL_TOPPASEDGE_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtGui/QGraphicsItem>

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
    int getSourceOutParam();
    /// Returns the source output parameter name
    QString getSourceOutParamName();
    /// Sets the target input parameter index
    void setTargetInParam(int in);
    /// Returns the target input parameter index
    int getTargetInParam();
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

    /// Returns the point in the @p list that is nearest to @p origin
    QPointF nearestPoint_(const QPointF & origin, const QList<QPointF> & list) const;
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

#endif
