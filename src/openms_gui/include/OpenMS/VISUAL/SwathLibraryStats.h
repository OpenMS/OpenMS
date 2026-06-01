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

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <QWidget>

namespace Ui
{
  class SwathLibraryStats;
}

namespace OpenMS
{
  
  /// A multi-tabbed widget for the SwathWizard offering setting of parameters, input-file specification and running Swath and more
  class OPENMS_GUI_DLLAPI SwathLibraryStats : public QWidget
  {
    Q_OBJECT

  public:
    explicit SwathLibraryStats(QWidget *parent = nullptr);
    ~SwathLibraryStats();

    /// updates the view with new information immediately
    void update(const TargetedExperiment::SummaryStatistics& stats);

    /// loads the pqp into a TargetedExperiment object while displaying a progress bar and shows the results when ready
    void updateFromFile(const QString& pqp_file);

  private slots:
  
  private:
    Ui::SwathLibraryStats* ui_;
  };

} // ns OpenMS

// this is required to allow Ui_SwathLibraryStats (auto UIC'd from .ui) to have a SwathLibraryStats member
using SwathLibraryStats = OpenMS::SwathLibraryStats;
