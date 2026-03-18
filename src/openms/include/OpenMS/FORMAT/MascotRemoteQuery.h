// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <string>
#include <vector>

typedef void CURL;

namespace OpenMS
{
  /**
      @brief Class which handles the communication between OpenMS and the Mascot server

      This class provides a communication interface which is able to query the Mascot
      server and reports the identifications provided be the Mascot server

      @htmlinclude OpenMS_MascotRemoteQuery.parameters

  */
  class MascotRemoteQuery :
    public DefaultParamHandler
  {
public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    OPENMS_DLLAPI MascotRemoteQuery();

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
    OPENMS_DLLAPI const std::string& getMascotXMLResponse() const;

    /// returns the Mascot XML response which contains the decoy identifications (note: setExportDecoys must be set to true, otherwise result will be empty)
    OPENMS_DLLAPI const std::string& getMascotXMLDecoyResponse() const;

    /// predicate which returns true if an error occurred during the query
    OPENMS_DLLAPI bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    OPENMS_DLLAPI const String& getErrorMessage() const;

    /// returns the search number
    OPENMS_DLLAPI String getSearchIdentifier() const;

    /// request export of decoy summary and decoys (note: internal decoy search must be enabled in the MGF file passed to mascot)
    OPENMS_DLLAPI void setExportDecoys(const bool b);

    /// Execute the full Mascot workflow synchronously (login -> query -> export results).
    OPENMS_DLLAPI void run();

protected:

    OPENMS_DLLAPI void updateMembers_() override;

private:

    /// login to Mascot server
    void login(CURL* curl);

    /// execute query (upload file)
    void execQuery(CURL* curl);

    /// download result file and process response
    std::string getResults(CURL* curl, const std::string& results_path);

    /// perform an HTTP GET and return the response body
    std::string httpGet(CURL* curl, const std::string& path);

    /// perform an HTTP POST with the given body and content type, return the response body
    std::string httpPost(CURL* curl, const std::string& path, const std::string& body, const std::string& content_type);

    /// build full URL from a path
    std::string buildUrl(const std::string& path) const;

    /// Extract search identifier from .dat file path
    OPENMS_DLLAPI String getSearchIdentifierFromFilePath(const String& path) const;

    /// Remove host name information from a url
    std::string removeHostName(const std::string& url) const;

    // Input / Output data
    String query_spectra_;
    std::string mascot_xml_;
    std::string mascot_decoy_xml_;

    // Internal data structures
    String error_message_;
    String search_identifier_;

    /// Path on mascot server
    String server_path_;
    /// Hostname of the mascot server
    String host_name_;
    /// Login required
    bool requires_login_ = false;
    /// Use SSL connection
    bool use_ssl_ = false;
    /// boundary string that will be embedded into the HTTP requests
    String boundary_;
    /// Timeout after these many seconds
    Int to_ = 1500;

    bool export_decoys_ = false;
  };

}
