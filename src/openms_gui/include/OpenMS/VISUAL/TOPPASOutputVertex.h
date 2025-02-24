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

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
  @brief A vertex representing an output, either folder or files(s)

  @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputVertex : public TOPPASVertex
  {
    Q_OBJECT
  public:
    /// Default C'tor
    TOPPASOutputVertex() = default;
    /// Copy constructor
    TOPPASOutputVertex(const TOPPASOutputVertex& rhs);
    /// Assignment operator
    TOPPASOutputVertex& operator=(const TOPPASOutputVertex& rhs);

    // documented in base class
    void reset(bool reset_all_files = false) override;
    /// opens the folder containing the output data
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent*) override;
    /// Returns the full directory (including preceding output path as selected by user and a trailing '/')
    String getFullOutputDirectory() const;
    /// Returns the directory where the output files are stored (includes a trailing '/')
    String getOutputDir() const;
    /// Creates the output directory for this node (includes a trailing '/')
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

  signals:
    /// Emitted when an output file was written
    void outputFileWritten(const String& file);

    /// Emitted when user has changed the output folder name (i.e. output dir needs to be newly created and packages updates)
    void outputFolderNameChanged();

  public slots:
    // documented in base class
    void inEdgeHasChanged() override;

  protected:
    /// custom output folder name
    QString output_folder_name_;

    // convenience members, not required for operation, but for progress during copying
    int files_written_ = 0; ///< files that were already written
    int files_total_ = 0;   ///< total number of files from upstream
  };

} // namespace OpenMS

