// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/NetworkGetRequest.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <QtNetwork/QNetworkRequest>
#include <QtGui/QTextDocument>

using namespace std;

namespace OpenMS
{

  NetworkGetRequest::NetworkGetRequest(QObject* parent) :
    QObject(parent), reply_(nullptr)
  {
    manager_ = new QNetworkAccessManager(this);
  }

  NetworkGetRequest::~NetworkGetRequest() = default;

  void NetworkGetRequest::setUrl(const QUrl& url)
  {
    url_ = url;
  }

  void NetworkGetRequest::run()
  {
    if (reply_ == nullptr)
    {
      error_ = QNetworkReply::NoError;
      error_string_ = "";
      QNetworkRequest request;
      request.setUrl(url_);
      request.setHeader(QNetworkRequest::ContentTypeHeader, "text/plain");
      connect(manager_, SIGNAL(finished(QNetworkReply*)), this, SLOT(replyFinished(QNetworkReply*)));
      reply_ = manager_->get(request);
    }
  }

  void NetworkGetRequest::replyFinished(QNetworkReply* reply)
  {
    if (reply_ != nullptr)
    {
      error_ = reply->error();
      error_string_ = error_ != QNetworkReply::NoError ? reply->errorString() : "";
      response_bytes_ = reply->readAll(); // in case of error this will just read the error html from the server
      reply->close();
      reply->deleteLater();;
    }
    emit done();
  }

  void NetworkGetRequest::timeOut()
  {
    if (reply_ != nullptr)
    {
      error_ = QNetworkReply::TimeoutError;
      error_string_ = "TimeoutError: the connection to the remote server timed out";
      reply_->abort();
      reply_->close();
      reply_->deleteLater();
    }
    emit done();
  }

  const QByteArray& NetworkGetRequest::getResponseBinary() const
  {
    return response_bytes_;
  }

  QString NetworkGetRequest::getResponse() const
  {
    return QString(response_bytes_);
  }  

  bool NetworkGetRequest::hasError() const
  {
    return error_ != QNetworkReply::NoError;
  }

  QString NetworkGetRequest::getErrorString() const
  {
    return error_string_;
  }

}
