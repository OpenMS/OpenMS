// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASOutputVertex.h>

namespace OpenMS
{
  /**
      @brief A vertex representing an output folder

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputFolderVertex :
    public TOPPASOutputVertex
  {
    Q_OBJECT

public:
    virtual std::unique_ptr<TOPPASVertex> clone() const override;
    /// returns "OutputFolderVertex"
    String getName() const override;
    // documented in base class
    void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    
    /// Called when the parent node has finished execution
    void run() override;
  };
} //namespace OpenMS
