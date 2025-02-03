// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QTimer>
#include <QtNetwork/QNetworkAccessManager>
#include <QtNetwork/QNetworkReply>


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

    /// assignment operator
    OPENMS_DLLAPI MascotRemoteQuery& operator=(const MascotRemoteQuery& rhs) = delete;

    /// copy constructor
    OPENMS_DLLAPI MascotRemoteQuery(const MascotRemoteQuery& rhs) = delete;

    /// destructor
    OPENMS_DLLAPI ~MascotRemoteQuery() override;
    //@}

    /// sets the query spectra, given in MGF file format
    OPENMS_DLLAPI void setQuerySpectra(const String& exp);

    /// returns the Mascot XML response which contains the identifications
    OPENMS_DLLAPI const QByteArray& getMascotXMLResponse() const;

    /// returns the Mascot XML response which contains the decoy identifications (note: setExportDecoys must be set to true, otherwise result will be empty)
    OPENMS_DLLAPI const QByteArray& getMascotXMLDecoyResponse() const;

    /// predicate which returns true if an error occurred during the query
    OPENMS_DLLAPI bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    OPENMS_DLLAPI const String& getErrorMessage() const;

    /// returns the search number
    OPENMS_DLLAPI String getSearchIdentifier() const;

    /// request export of decoy summary and decoys (note: internal decoy search must be enabled in the MGF file passed to mascot)
    OPENMS_DLLAPI void setExportDecoys(const bool b);

protected:

    OPENMS_DLLAPI void updateMembers_() override;

public slots:

    OPENMS_DLLAPI void run();

private slots:

    /// slot connected to QTimer (timeout_)
    OPENMS_DLLAPI void timedOut() const;

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
    void getResults(const QString& results_path);

    /// finish a run and emit "done"
    OPENMS_DLLAPI void endRun_();

    /**
      @brief Remove host name information from an url, e.g., "http://www.google.de/search" -> "search"

      @param url The url that will be manipulated.
    */
    void removeHostName_(QString& url);

    /// helper function to build URL
    QUrl buildUrl_(const std::string& path);

    /// Write HTTP header to error stream (for debugging)
    OPENMS_DLLAPI void logHeader_(const QNetworkRequest& header, const String& what);

    /// Write HTTP header to error stream (for debugging)
    OPENMS_DLLAPI void logHeader_(const QNetworkReply* header, const String& what);

    OPENMS_DLLAPI String getSearchIdentifierFromFilePath(const String& path) const;

    /// parse new response header
    OPENMS_DLLAPI void readResponseHeader(const QNetworkReply* reply);

    QNetworkAccessManager* manager_;

    // Input / Output data
    String query_spectra_;
    QByteArray mascot_xml_;
    QByteArray mascot_decoy_xml_;

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

    bool export_decoys_ = false;
  };

}

