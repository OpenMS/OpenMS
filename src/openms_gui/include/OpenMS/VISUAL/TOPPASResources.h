// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_VISUAL_TOPPASRESOURCES_H
#define OPENMS_VISUAL_TOPPASRESOURCES_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/TOPPASResource.h>

#include <QtCore/QString>
#include <QtCore/QObject>

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
    Map<QString, QList<TOPPASResource> > map_;
    /// The empty list
    QList<TOPPASResource> empty_list_;
  };
}

#endif
