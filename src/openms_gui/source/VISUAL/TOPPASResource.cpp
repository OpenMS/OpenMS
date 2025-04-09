// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <iostream>
#include <OpenMS/VISUAL/TOPPASResource.h>

namespace OpenMS
{
  QStringList TOPPASResource::supported_schemes = (QStringList() << "file");

  TOPPASResource::TOPPASResource(const QString & file) :
    QObject(),
    url_(),
    file_name_("")
  {
    fromLocalFile(file);
  }

  TOPPASResource::TOPPASResource(const QUrl & url) :
    QObject(),
    url_(),
    file_name_("")
  {
    QString scheme = url.scheme().toLower();
    if (!supported_schemes.contains(scheme))
    {
      std::cerr << "URL scheme not supported!" << std::endl;
    }
    else
    {
      url_ = url;

      if (scheme == "file")
      {
        file_name_ = url.toLocalFile();
      }
    }
  }

  TOPPASResource::TOPPASResource(const TOPPASResource & rhs) :
    QObject(),
    url_(rhs.url_),
    file_name_(rhs.file_name_)
  {
  }

  TOPPASResource::~TOPPASResource() = default;

  TOPPASResource & TOPPASResource::operator=(const TOPPASResource & rhs)
  {
    url_ = rhs.url_;
    file_name_ = rhs.file_name_;

    return *this;
  }

  void TOPPASResource::writeToFile(const QString & file_name)
  {
    // TODO retrieve data and write it to file_name

    file_name_ = file_name;
  }

  const QString & TOPPASResource::getLocalFile() const
  {
    return file_name_;
  }

  const QUrl & TOPPASResource::getURL() const
  {
    return url_;
  }

  void TOPPASResource::fromLocalFile(const QString & file)
  {
    url_ = QUrl::fromLocalFile(file);
    file_name_ = file;
  }

} //namespace
