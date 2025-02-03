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

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A vertex representing an input file list

      @ingroup TOPPAS_elements
  */
class OPENMS_GUI_DLLAPI TOPPASInputFileListVertex : public TOPPASVertex
{
  Q_OBJECT

public:
  /// Default constructor
  TOPPASInputFileListVertex() = default;
  /// Constructor
  TOPPASInputFileListVertex(const QStringList& files);
  /// Copy constructor
  TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs) = default;
  /// Destructor
  ~TOPPASInputFileListVertex() override = default;
  /// Assignment operator
  TOPPASInputFileListVertex& operator=(const TOPPASInputFileListVertex& rhs) = default;

  virtual std::unique_ptr<TOPPASVertex> clone() const override;

  /// returns "InputVertex"
  String getName() const override;
  /// Sets the list of files
  void setFilenames(const QStringList & files);
  /// Starts all tools below this node
  void run() override;
  // documented in base class
  void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override;
  // documented in base class
  QRectF boundingRect() const override;
  /// Checks if the given list of file names is valid
  bool fileNamesValid();
  /// Shows the dialog for editing the files
  void showFilesDialog();
  /// Opens the folders of the input files
  void openContainingFolder();
  /// Returns the key (for applying resources from a resource file)
  const QString & getKey();
  /// Sets the key (for applying resources from a resource file)
  void setKey(const QString & key);

public slots:
    /// Called by an outgoing edge when it has changed
    void outEdgeHasChanged() override;

protected:

    /// The key of this input node (for applying resources from a resource file)
    QString key_;

    /// current working dir, i.e. the last position a file was added from
    QString cwd_;

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent * e) override;
    //@}

  };
}

