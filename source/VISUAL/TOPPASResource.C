// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Johannes Junker $
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

  TOPPASResource::~TOPPASResource()
  {
  }

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
