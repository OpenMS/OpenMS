// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Johannes Junker, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A special vertex that allows to split a list of inputs.

      Tools that produce lists of output files (several files in each processing round, e.g. map alignment tools) cannot directly provide input for tools that only take a single input file in TOPPAS. This "Splitter" node provides the necessary glue, by splitting a list of output files into several rounds of single input files.

      See the @ref OpenMS::TOPPASMergerVertex "Collector" node for the opposite operation.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASSplitterVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASSplitterVertex() = default;
    /// Copy constructor
    TOPPASSplitterVertex(const TOPPASSplitterVertex& rhs);
    /// Destructor
    ~TOPPASSplitterVertex() override = default;
    /// Assignment operator
    TOPPASSplitterVertex& operator=(const TOPPASSplitterVertex& rhs);
    virtual std::unique_ptr<TOPPASVertex> clone() const override;
    /// returns "SplitterVertex"
    String getName() const override;
    /// check if upstream nodes are finished and call downstream nodes
    void run() override;
    // documented in base class
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    void markUnreachable() override;

protected:

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e) override;
    //@}

  };
}

