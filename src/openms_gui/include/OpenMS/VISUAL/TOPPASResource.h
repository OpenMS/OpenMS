// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtCore/QUrl>
#include <QtCore/QObject>

namespace OpenMS
{
  /**
      @brief Represents a data resource for TOPPAS workflows.

      Currently, the only supported type of resource is local files.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASResource :
    QObject
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASResource(const QString & file);
    /// Constructor from URL
    TOPPASResource(const QUrl & url);
    /// Copy constructor
    TOPPASResource(const TOPPASResource & rhs);
    /// Destructor
    ~TOPPASResource() override;
    /// Assignment operator
    TOPPASResource & operator=(const TOPPASResource & rhs);
    /// Writes this resource to the local file @p file
    void writeToFile(const QString & file_name);
    /// Returns the file name of the local file, or "" if it has not been written yet
    const QString & getLocalFile() const;
    /// Returns the URL of this resource
    const QUrl & getURL() const;
    /// Sets the URL of this resource from @p file
    void fromLocalFile(const QString & file);

    /// Supported schemes
    static QStringList supported_schemes;

protected:

    /// The URL of this resource
    QUrl url_;
    /// The name of the local file
    QString file_name_;
  };
}

