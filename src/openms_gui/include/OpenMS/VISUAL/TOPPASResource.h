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

#ifndef OPENMS_VISUAL_TOPPASRESOURCE_H
#define OPENMS_VISUAL_TOPPASRESOURCE_H

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

#endif
