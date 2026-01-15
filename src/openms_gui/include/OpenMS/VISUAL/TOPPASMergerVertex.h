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
      @brief A special vertex that allows to merge several inputs.

      A special vertex that allows to merge several inputs. Mergers have two modes: The normal,
      round-based merging mode and a "wait & merge all" mode. In round-based mode, a merger
      first takes the first files of each incoming file list and merges them into a list (which
      has as many elements as the merger has incoming edges).

      In "wait & merge all" mode, the merger first waits for all upstream mergers to finish all
      their merging rounds and then merges all collected files from all merging rounds for all
      incoming edges into one single list and calls the next tool with this list of files as input.

      @ingroup TOPPAS_elements
  */
class OPENMS_GUI_DLLAPI TOPPASMergerVertex : public TOPPASVertex
{
  Q_OBJECT

public:
  /// Default constructor
  TOPPASMergerVertex() = default;
  /// Constructor
  TOPPASMergerVertex(bool round_based);
  /// Copy constructor
  TOPPASMergerVertex(const TOPPASMergerVertex& rhs) = default;
  /// Destructor
  ~TOPPASMergerVertex() override = default;
  /// Assignment operator
  TOPPASMergerVertex& operator=(const TOPPASMergerVertex& rhs) = default;

  virtual std::unique_ptr<TOPPASVertex> clone() const override;

  /// returns "MergerVertex"
  String getName() const override;
  /// check if upstream nodes are finished and call downstream nodes
  void run() override;
  /// Determines whether this merger is merging round based or merging all inputs into one list
  bool roundBasedMode() const;
  // documented in base class
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) override;
  // documented in base class
  QRectF boundingRect() const override;
  // documented in base class
  void markUnreachable() override;

public slots:

signals:
    /// Emitted when merging upstream data failed
    void mergeFailed(const QString message);

protected:

    /// Stores whether this merger is merging round based or merging all inputs into one list
    bool round_based_mode_{true};
  };
}

