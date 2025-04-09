// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QUrl>
#include <QtNetwork/QNetworkReply>

namespace OpenMS
{

  class NetworkGetRequest :
    public QObject
  {
    Q_OBJECT

  public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    OPENMS_DLLAPI NetworkGetRequest(QObject* parent = nullptr);

    /// destructor
    OPENMS_DLLAPI ~NetworkGetRequest() override;
    //@}

    // set request parameters
    OPENMS_DLLAPI void setUrl(const QUrl& url);

    /// returns the response
    OPENMS_DLLAPI QString getResponse() const;

    /// returns the response
    OPENMS_DLLAPI const QByteArray& getResponseBinary() const;

    /// returns true if an error occurred during the query
    OPENMS_DLLAPI bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    OPENMS_DLLAPI QString getErrorString() const;

  protected:

    public slots:

    OPENMS_DLLAPI void run();

    OPENMS_DLLAPI void timeOut();

    private slots:

    OPENMS_DLLAPI void replyFinished(QNetworkReply*);

  signals:

    OPENMS_DLLAPI void done();

  private:
    /// assignment operator
    OPENMS_DLLAPI NetworkGetRequest& operator=(const NetworkGetRequest& rhs);
    /// copy constructor
    OPENMS_DLLAPI NetworkGetRequest(const NetworkGetRequest& rhs);

    QByteArray response_bytes_;
    QUrl url_;
    QNetworkAccessManager* manager_;
    QNetworkReply* reply_;
    QNetworkReply::NetworkError error_;
    QString error_string_;
  };
}

