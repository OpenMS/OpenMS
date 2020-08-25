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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtNetwork/QNetworkAccessManager>
#include <QTimer>
#include <QNetworkReply>


namespace OpenMS
{
  /**
      @brief Class which handles the communication between OpenMS and the Mascot server

      This class provides a communication interface which is able to query the Mascot
      server and reports the identifications provided be the Mascot server

      @htmlinclude OpenMS_MascotRemoteQuery.parameters

  */
  class MascotRemoteQuery :
    public QObject,
    public DefaultParamHandler
  {
    Q_OBJECT

public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    OPENMS_DLLAPI MascotRemoteQuery(QObject* parent = 0);

    /// destructor
    OPENMS_DLLAPI virtual ~MascotRemoteQuery();
    //@}

    /// sets the query spectra, given in MGF file format
    OPENMS_DLLAPI void setQuerySpectra(const String& exp);

    /// returns the Mascot XML response which contains the identifications
    OPENMS_DLLAPI const QByteArray& getMascotXMLResponse() const;

    /// predicate which returns true if an error occurred during the query
    OPENMS_DLLAPI bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    OPENMS_DLLAPI const String& getErrorMessage() const;

    /// returns the search number
    OPENMS_DLLAPI String getSearchIdentifier() const;

protected:

    OPENMS_DLLAPI virtual void updateMembers_();

public slots:

    OPENMS_DLLAPI void run();

private slots:

    /// slot connected to QTimer (timeout_)
    OPENMS_DLLAPI void timedOut();

    /// slot connected to the QNetworkAccessManager::finished signal
    OPENMS_DLLAPI void readResponse(QNetworkReply* reply);

    /// slot connected to signal downloadProgress
    OPENMS_DLLAPI void downloadProgress(qint64 bytes_read, qint64 bytes_total);

    /// slot connected to signal uploadProgress
    OPENMS_DLLAPI void uploadProgress(qint64 bytes_read, qint64 bytes_total);

    /// slot connected to signal gotRedirect
    OPENMS_DLLAPI void followRedirect(QNetworkReply * reply);

signals:

    /// signal when class got a redirect
    OPENMS_DLLAPI void gotRedirect(QNetworkReply * reply);

    /// signal when class is done and results can be collected
    OPENMS_DLLAPI void done();

private:

    /// login to Mascot server
    void login();

    /// execute query (upload file)
    void execQuery();

    /// download result file
    void getResults(QString results_path);

    /// assignment operator
    OPENMS_DLLAPI MascotRemoteQuery& operator=(const MascotRemoteQuery& rhs);

    /// copy constructor
    OPENMS_DLLAPI MascotRemoteQuery(const MascotRemoteQuery& rhs);

    /// finish a run and emit "done"
    OPENMS_DLLAPI void endRun_();

    /**
      @brief Remove host name information from an url, e.g., "http://www.google.de/search" -> "search"

      @param The url that will be manipulated.
    */
    void removeHostName_(QString& url);

    /// helper function to build URL
    QUrl buildUrl_(std::string path);

    /// Write HTTP header to error stream (for debugging)
    OPENMS_DLLAPI void logHeader_(const QNetworkRequest header, const String& what);

    /// Write HTTP header to error stream (for debugging)
    OPENMS_DLLAPI void logHeader_(const QNetworkReply* header, const String& what);

    OPENMS_DLLAPI String getSearchIdentifierFromFilePath(const String& path) const;

    /// parse new response header
    OPENMS_DLLAPI void readResponseHeader(const QNetworkReply* reply);

    QNetworkAccessManager* manager_;

    // Input / Output data
    String query_spectra_;
    QByteArray mascot_xml_;

    // Internal data structures
    QString cookie_;
    String error_message_;
    QTimer timeout_;
    String search_identifier_;

    /// Path on mascot server
    String server_path_;
    /// Hostname of the mascot server
    String host_name_;
    /// Login required
    bool requires_login_;
    /// Use SSL connection
    bool use_ssl_;
    /// boundary string that will be embedded into the HTTP requests
    String boundary_;
    /// Timeout after these many seconds
    Int to_;
  };

}

