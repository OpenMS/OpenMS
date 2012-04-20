// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{
  TOPPASResources::TOPPASResources() :
    QObject(),
    map_(),
    empty_list_()
  {
  }

  TOPPASResources::TOPPASResources(const TOPPASResources & rhs) :
    QObject(),
    map_(rhs.map_),
    empty_list_()
  {
  }

  TOPPASResources::~TOPPASResources()
  {
  }

  TOPPASResources & TOPPASResources::operator=(const TOPPASResources & rhs)
  {
    map_ = rhs.map_;

    return *this;
  }

  void TOPPASResources::load(const QString & file_name)
  {
    Param load_param;
    load_param.load(String(file_name));

    for (Param::ParamIterator it = load_param.begin(); it != load_param.end(); ++it)
    {
      StringList substrings;
      it.getName().split(':', substrings);
      if (substrings.size() != 2 ||
          substrings.back() != "url_list" ||
          (it->value).valueType() != DataValue::STRING_LIST)
      {
        std::cerr << "Invalid file format." << std::endl;
        return;
      }

      QString key = (substrings[0]).toQString();
      StringList url_list = (StringList)(it->value);
      QList<TOPPASResource> resource_list;
      for (StringList::ConstIterator it = url_list.begin(); it != url_list.end(); ++it)
      {
        resource_list << TOPPASResource(QUrl(it->toQString()));
      }

      add(key, resource_list);
    }
  }

  void TOPPASResources::add(const QString & key, const QList<TOPPASResource> & resource_list)
  {
    map_[key] = resource_list;
  }

  void TOPPASResources::store(const QString & file_name)
  {
    Param save_param;

    for (Map<QString, QList<TOPPASResource> >::ConstIterator it = map_.begin(); it != map_.end(); ++it)
    {
      const String & key = String(it->first);
      const QList<TOPPASResource> & resource_list = it->second;
      StringList url_list;
      foreach(const TOPPASResource &res, resource_list)
      {
        url_list << String(res.getURL().toString());
      }
      save_param.setValue(key + ":url_list", DataValue(url_list));
    }

    save_param.store(String(file_name));
  }

  const QList<TOPPASResource> & TOPPASResources::get(const QString & key) const
  {
    if (!map_.has(key))
    {
      return empty_list_;
    }

    return map_[key];
  }

  void TOPPASResources::clear()
  {
    map_.clear();
  }

} //namespace
