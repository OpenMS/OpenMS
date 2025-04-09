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

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/TOPPASResource.h>

#include <QtCore/QString>
#include <QtCore/QObject>

#include <map>

namespace OpenMS
{
  /**
      @brief A dictionary mapping string keys to lists of TOPPASResource objects

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASResources :
    QObject
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASResources();
    /// Copy constructor
    TOPPASResources(const TOPPASResources & rhs);
    /// Destructor
    ~TOPPASResources() override;
    /// Assignment operator
    TOPPASResources & operator=(const TOPPASResources & rhs);
    /// Adds the (key,resource_list) pair to the dictionary
    void add(const QString & key, const QList<TOPPASResource> & resource_list);
    /// Returns the resource list that @p key is mapped to, or an empty list if @p key does not exist
    const QList<TOPPASResource> & get(const QString & key) const;
    /// Loads the dictionary from file @p file_name
    void load(const QString & file_name);
    /// Writes the dictionary to file @p file_name
    void store(const QString & file_name);
    /// Clears the dictionary
    void clear();

protected:

    /// The dictionary
    std::map<QString, QList<TOPPASResource> > map_;
    /// The empty list
    QList<TOPPASResource> empty_list_;
  };
}

