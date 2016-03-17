// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/NetworkGetRequest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <QtGui/QTextDocument>
#include <QtNetwork/QNetworkRequest>

#include <iostream>

// #define QUERY_DEBUG
using namespace std;

namespace OpenMS
{

  NetworkGetRequest::NetworkGetRequest(QObject* parent) :
    QObject(parent)
  {
    to_ = 5;
    timeout_.setInterval(1000 * to_);
  }

  NetworkGetRequest::~NetworkGetRequest()
  {
  }

  void NetworkGetRequest::timedOut()
  {
  }

  void NetworkGetRequest::setUrl(QUrl url)
  {
    url_ = url;
  }

  void NetworkGetRequest::run()
  {
    QNetworkRequest request;
    request.setUrl(url_);
    request.setHeader(QNetworkRequest::ContentTypeHeader, "text/plain");
    QNetworkAccessManager* nam = new QNetworkAccessManager();
    connect(nam, SIGNAL(finished(QNetworkReply*)), this, SLOT(replyFinished(QNetworkReply*)));
    QNetworkReply* nr = nam->get(request);

  }

  void NetworkGetRequest::replyFinished(QNetworkReply* reply)
  {
    response_bytes_ = reply->readAll();
    emit done();
  }

  QString NetworkGetRequest::getResponse() const
  {
    return QString(response_bytes_);
  }

  void NetworkGetRequest::endRun_()
  {
    emit done();
  }

  bool NetworkGetRequest::hasError() const
  {
    return error_message_ != "";
  }

  const String& NetworkGetRequest::getErrorMessage() const
  {
    return error_message_;
  }

}
