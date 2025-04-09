// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <iostream>
#include <map>

namespace OpenMS
{
  TOPPASResources::TOPPASResources() :
    QObject(),
    map_(),
    empty_list_()
  {
  }

  TOPPASResources::TOPPASResources(const TOPPASResources& rhs) :
    QObject(),
    map_(rhs.map_),
    empty_list_()
  {
  }

  TOPPASResources::~TOPPASResources() = default;

  TOPPASResources& TOPPASResources::operator=(const TOPPASResources& rhs)
  {
    map_ = rhs.map_;

    return *this;
  }

  void TOPPASResources::load(const QString& file_name)
  {
    Param load_param;
    ParamXMLFile paramFile;
    paramFile.load(String(file_name), load_param);

    for (Param::ParamIterator it = load_param.begin(); it != load_param.end(); ++it)
    {
      StringList substrings;
      String(it.getName()).split(':', substrings);
      if (substrings.size() != 2 ||
          substrings.back() != "url_list" ||
          (it->value).valueType() != ParamValue::STRING_LIST)
      {
        std::cerr << "Invalid file format." << std::endl;
        return;
      }

      QString key = (substrings[0]).toQString();
      StringList url_list = ListUtils::toStringList<std::string>(it->value);
      QList<TOPPASResource> resource_list;
      for (StringList::const_iterator it = url_list.begin(); it != url_list.end(); ++it)
      {
        resource_list << TOPPASResource(QUrl(it->toQString()));
      }

      add(key, resource_list);
    }
  }

  void TOPPASResources::add(const QString& key, const QList<TOPPASResource>& resource_list)
  {
    map_[key] = resource_list;
  }

  void TOPPASResources::store(const QString& file_name)
  {
    Param save_param;

    for (std::map<QString, QList<TOPPASResource> >::const_iterator it = map_.begin(); it != map_.end(); ++it)
    {
      const String& key = String(it->first);
      const QList<TOPPASResource>& resource_list = it->second;
      std::vector<std::string> url_list;
      foreach(const TOPPASResource &res, resource_list)
      {
        url_list.push_back(String(res.getURL().toString().toStdString()));
      }
      save_param.setValue(key + ":url_list", url_list);
    }

    ParamXMLFile paramFile;
    paramFile.store(String(file_name), save_param);
  }

  const QList<TOPPASResource>& TOPPASResources::get(const QString& key) const
  {
    if (map_.find(key) == map_.end())
    {
      return empty_list_;
    }

    return map_.at(key);
  }

  void TOPPASResources::clear()
  {
    map_.clear();
  }

} //namespace
