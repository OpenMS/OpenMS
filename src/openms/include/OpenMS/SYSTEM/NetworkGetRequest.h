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

#ifndef OPENMS_SYSTEM_NETWORKGETREQUEST_H
#define OPENMS_SYSTEM_NETWORKGETREQUEST_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtNetwork/QHttpRequestHeader>
#include <QtCore/QTimer>
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
    OPENMS_DLLAPI NetworkGetRequest(QObject* parent = 0);

    /// destructor
    OPENMS_DLLAPI virtual ~NetworkGetRequest();
    //@}

    // set request parameters
    OPENMS_DLLAPI void setUrl(QUrl url);

    /// returns the response
    OPENMS_DLLAPI QString getResponse() const;

    /// returns true if an error occurred during the query
    OPENMS_DLLAPI bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    OPENMS_DLLAPI const String& getErrorMessage() const;

  protected:

    public slots:

    OPENMS_DLLAPI void run();

    private slots:

    OPENMS_DLLAPI void replyFinished(QNetworkReply*);

    OPENMS_DLLAPI void timedOut();

  signals:

    OPENMS_DLLAPI void done();

  private:
    /// assignment operator
    OPENMS_DLLAPI NetworkGetRequest& operator=(const NetworkGetRequest& rhs);
    /// copy constructor
    OPENMS_DLLAPI NetworkGetRequest(const NetworkGetRequest& rhs);

    OPENMS_DLLAPI void endRun_();
    QByteArray response_bytes_;
    QUrl url_;
    String error_message_;
    QTimer timeout_;
    Int to_;
  };
}

#endif
