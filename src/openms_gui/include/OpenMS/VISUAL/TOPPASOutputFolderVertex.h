// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

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
    /// returns "OutputFolderVertex"
    String getName() const override;
    // documented in base class
    void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    void reset(bool reset_all_files = false) override;
    /// opens the folder containing the output data
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent*) override;
    /// Called when the parent node has finished execution
    void run() override;
    /// Returns the full directory (including preceding output path as selected by user)
    String getFullOutputDirectory() const;
    /// Returns the directory where the output files are stored
    String getOutputDir() const;
    /// Creates the output directory for this node
    String createOutputDir() const;
    /// Sets the topological sort number and removes invalidated tmp files
    void setTopoNr(UInt nr) override;
    /// Opens the folders of the output files
    void openContainingFolder() const;
    /// Sets a custom output folder name, which will be integrated into 'getOutputDir()' and 'getFullOutputDirectory()' calls.
    /// @note The string is not checked for validity (avoid characters which are not allowed in directories, e.g. '{')
    void setOutputFolderName(const QString& name);
    /// return the output folder where results are written
    const QString& getOutputFolderName() const;

public slots:
    //documented in base class
    void inEdgeHasChanged() override;

  };
}

