// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
  /**
      @brief A vertex representing an output file list

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputFileListVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASOutputFileListVertex() = default;
    /// Copy constructor
    TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex & rhs);
    /// Destructor
    ~TOPPASOutputFileListVertex() override = default;
    /// Assignment operator
    TOPPASOutputFileListVertex & operator=(const TOPPASOutputFileListVertex & rhs);
    /// returns "OutputVertex"
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

signals:
    /// Emitted when an output file was written
    void outputFileWritten(const String& file);

    /// Emitted when user has changed the output folder name (i.e. output dir needs to be newly created and packages updates)
    void outputFolderNameChanged();

protected:

    // custom output folder name
    QString output_folder_name_;

    static bool copy_(const QString & from, const QString & to); ///< STATIC(!) function which calls QFile::copy(); needs to be static, since we need to pass a function pointer (which does not work on member functions)
    // convenience members, not required for operation, but for progress during copying
    int files_written_ = 0;   ///< files that were already written
    int files_total_ = 0;     ///< total number of files from upstream
  };
}

