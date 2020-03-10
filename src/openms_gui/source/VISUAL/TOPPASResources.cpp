// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <iostream>

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

  TOPPASResources::~TOPPASResources()
  {
  }

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
      it.getName().split(':', substrings);
      if (substrings.size() != 2 ||
          substrings.back() != "url_list" ||
          (it->value).valueType() != DataValue::STRING_LIST)
      {
        std::cerr << "Invalid file format." << std::endl;
        return;
      }

      QString key = (substrings[0]).toQString();
      StringList url_list = (it->value);
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

    for (Map<QString, QList<TOPPASResource> >::ConstIterator it = map_.begin(); it != map_.end(); ++it)
    {
      const String& key = String(it->first);
      const QList<TOPPASResource>& resource_list = it->second;
      StringList url_list;
      foreach(const TOPPASResource &res, resource_list)
      {
        url_list.push_back(String(res.getURL().toString()));
      }
      save_param.setValue(key + ":url_list", DataValue(url_list));
    }

    ParamXMLFile paramFile;
    paramFile.store(String(file_name), save_param);
  }

  const QList<TOPPASResource>& TOPPASResources::get(const QString& key) const
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
