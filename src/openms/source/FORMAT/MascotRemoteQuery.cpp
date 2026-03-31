// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/CurlInit.h>

#include <curl/curl.h>

#include <regex>

// #define MASCOTREMOTEQUERY_DEBUG
// #define MASCOTREMOTEQUERY_DEBUG_FULL_QUERY

using namespace std;

namespace OpenMS
{
  namespace
  {
    size_t writeCallback(char* ptr, size_t size, size_t nmemb, void* userdata)
    {
      auto* buf = static_cast<string*>(userdata);
      buf->append(ptr, size * nmemb);
      return size * nmemb;
    }
  }

  MascotRemoteQuery::MascotRemoteQuery() :
    DefaultParamHandler("MascotRemoteQuery")
  {
    // server specifications
    defaults_.setValue("hostname", "", "Address of the host where Mascot listens, e.g. 'mascot-server' or '127.0.0.1'");
    defaults_.setValue("host_port", 80, "Port where the Mascot server listens, 80 should be a good guess");
    defaults_.setMinInt("host_port", 0);
    defaults_.setValue("server_path", "mascot", "Path on the server where Mascot server listens, 'mascot' should be a good guess");
    defaults_.setValue("timeout", 1500, "Timeout in seconds, after which the query is declared as failed."
                                        "This is NOT the whole time the search takes, but the time in between two progress steps. Some Mascot servers freeze during this (unstable network etc) and idle forever"
                                        ", the connection is killed. Set this to 0 to disable timeout!");
    defaults_.setMinInt("timeout", 0);
    defaults_.setValue("boundary", "GZWgAaYKjHFeUaLOLEIOMq", "Boundary for the MIME section", {"advanced"});

    // proxy settings
    defaults_.setValue("use_proxy", "false", "Flag which enables the proxy usage for the HTTP requests, please specify at least 'proxy_host' and 'proxy_port'", {"advanced"});
    defaults_.setValidStrings("use_proxy", {"true","false"});
    defaults_.setValue("proxy_host", "", "Host where the proxy server runs on", {"advanced"});
    defaults_.setValue("proxy_port", 0, "Port where the proxy server listens", {"advanced"});
    defaults_.setMinInt("proxy_port", 0);
    defaults_.setValue("proxy_username", "", "Login name for the proxy server, if needed", {"advanced"});
    defaults_.setValue("proxy_password", "", "Login password for the proxy server, if needed", {"advanced"});

    // login for Mascot security
    defaults_.setValue("login", "false", "Flag which should be set 'true' if Mascot security is enabled; also set 'username' and 'password' then.");
    defaults_.setValidStrings("login", {"true","false"});
    defaults_.setValue("username", "", "Name of the user if login is used (Mascot security must be enabled!)");
    defaults_.setValue("password", "", "Password of the user if login is used (Mascot security must be enabled!)");
    defaults_.setValue("use_ssl", "false", "Flag indicating whether you want to send requests to an HTTPS server or not (HTTP). Requires OpenSSL to be installed (see openssl.org)");
    defaults_.setValidStrings("use_ssl", {"true","false"});

    // Mascot export options
    defaults_.setValue("export_params", "_ignoreionsscorebelow=0&_sigthreshold=0.99&_showsubsets=1&show_same_sets=1&report=0&percolate=0&query_master=0", "Adjustable export parameters (passed to Mascot's 'export_dat_2.pl' script). Generally only parameters that control which hits to export are safe to adjust/add. Many settings that govern what types of information to include are required by OpenMS and cannot be changed. Note that setting 'query_master' to 1 may lead to incorrect protein references for peptides.", {"advanced"});
    defaults_.setValue("skip_export", "false", "For use with an external Mascot Percolator: Run the Mascot search, but do not export the results. The output file produced by MascotAdapterOnline will contain only the Mascot search number.", {"advanced"});
    defaults_.setValidStrings("skip_export", {"true","false"});
    defaults_.setValue("batch_size", 50000, "Number of spectra processed in one batch by Mascot (default 50000)", {"advanced"});
    defaults_.setMinInt("batch_size", 0);
    defaultsToParam_();
  }

  MascotRemoteQuery::~MascotRemoteQuery()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    std::cerr << "MascotRemoteQuery::~MascotRemoteQuery()\n";
#endif
  }

  void MascotRemoteQuery::run()
  {
    updateMembers_();
    CurlInit::ensure();

    CURL* curl = curl_easy_init();
    if (!curl)
    {
      error_message_ = "Failed to initialize libcurl";
      return;
    }

    // Enable cookie engine (automatic cookie management)
    curl_easy_setopt(curl, CURLOPT_COOKIEFILE, "");
    // Thread-safe DNS timeout (avoids signal handler issues with OpenMP)
    curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);
    // Follow HTTP redirects automatically
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    // SSL verification
    if (use_ssl_)
    {
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    }
    // Timeout
    if (to_ > 0)
    {
      curl_easy_setopt(curl, CURLOPT_TIMEOUT, static_cast<long>(to_));
    }

    // Proxy
    bool use_proxy(param_.getValue("use_proxy").toBool());
    if (use_proxy)
    {
      String proxy_host(param_.getValue("proxy_host").toString());
      Int proxy_port(param_.getValue("proxy_port"));
      string proxy_url = proxy_host + ":" + String(proxy_port);
      curl_easy_setopt(curl, CURLOPT_PROXY, proxy_url.c_str());
      // SOCKS5 matches the original Qt implementation (QNetworkProxy::Socks5Proxy)
      curl_easy_setopt(curl, CURLOPT_PROXYTYPE, static_cast<long>(CURLPROXY_SOCKS5));

      String proxy_username(param_.getValue("proxy_username").toString());
      if (!proxy_username.empty())
      {
        curl_easy_setopt(curl, CURLOPT_PROXYUSERNAME, proxy_username.c_str());
        String proxy_password(param_.getValue("proxy_password").toString());
        curl_easy_setopt(curl, CURLOPT_PROXYPASSWORD, proxy_password.c_str());
      }
    }

    // Step 1: Optional login
    if (requires_login_)
    {
      login(curl);
      if (hasError())
      {
        curl_easy_cleanup(curl);
        return;
      }
    }

    // Step 2: Execute query
    execQuery(curl);
    if (hasError())
    {
      curl_easy_cleanup(curl);
      return;
    }

    curl_easy_cleanup(curl);
  }

  string MascotRemoteQuery::buildUrl(const string& path) const
  {
    string protocol = use_ssl_ ? "https" : "http";
    Int port = param_.getValue("host_port");
    return protocol + "://" + string(host_name_.c_str()) + ":" + to_string(port) + path;
  }

  string MascotRemoteQuery::httpGet(CURL* curl, const string& path)
  {
    string url = buildUrl(path);
    string response;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPGET, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nhttpGet URL: " << url << "\n";
#endif

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK)
    {
      error_message_ = String("Mascot HTTP GET failed: ") + curl_easy_strerror(res);
      return {};
    }

    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    if (http_code >= 400)
    {
      error_message_ = String("Mascot server returned HTTP error ") + String(http_code) + "\nTry accessing the server\n  " + host_name_ + server_path_ + "\n from your browser and check if it works fine.";
      return {};
    }

    return response;
  }

  string MascotRemoteQuery::httpPost(CURL* curl, const string& path, const string& body, const string& content_type)
  {
    string url = buildUrl(path);
    string response;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, static_cast<long>(body.size()));
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

    struct curl_slist* headers = nullptr;
    string ct_header = "Content-Type: " + content_type;
    headers = curl_slist_append(headers, ct_header.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nhttpPost URL: " << url << "\n";
#endif

    CURLcode res = curl_easy_perform(curl);
    curl_slist_free_all(headers);
    // Reset custom headers
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, nullptr);

    if (res != CURLE_OK)
    {
      error_message_ = String("Mascot HTTP POST failed: ") + curl_easy_strerror(res);
      return {};
    }

    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    if (http_code >= 400)
    {
      error_message_ = String("Mascot server returned HTTP error ") + String(http_code) + "\nTry accessing the server\n  " + host_name_ + server_path_ + "\n from your browser and check if it works fine.";
      return {};
    }

    return response;
  }

  void MascotRemoteQuery::login(CURL* curl)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::login()" << "\n";
#endif

    string boundary = string(boundary_.c_str());
    string body;
    string boundary_line = "--" + boundary + "\r\n";

    auto addField = [&](const string& name, const string& value)
    {
      body += boundary_line;
      body += "Content-Disposition: form-data; name=\"" + name + "\"\r\n";
      body += "\r\n";
      body += value + "\r\n";
    };

    addField("username", string(param_.getValue("username").toString()));
    addField("password", string(param_.getValue("password").toString()));
    addField("submit", "Login");
    addField("referer", "");
    addField("display", "nothing");
    addField("savecookie", "1");
    addField("action", "login");
    addField("userid", "");
    addField("onerrdisplay", "login_prompt");
    body += "--" + boundary + "--\r\n";

    string content_type = "multipart/form-data; boundary=" + boundary;
    string response = httpPost(curl, server_path_ + "/cgi/login.pl", body, content_type);

    if (hasError()) return;

    if (response.find("Logged in successfu") != string::npos)
    {
      OPENMS_LOG_INFO << "Login successful!" << std::endl;
    }
    else if (response.find("Error: You have entered an invalid password") != string::npos)
    {
      error_message_ = "Error: You have entered an invalid password";
    }
    else if (response.find("is not a valid user") != string::npos)
    {
      error_message_ = "Error: Username is not valid";
    }
    else
    {
      error_message_ = "Login failed with unexpected response from Mascot server";
    }
  }

  void MascotRemoteQuery::execQuery(CURL* curl)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::execQuery()" << "\n";
#endif

    string boundary = string(boundary_.c_str());
    string body;
    body += "--" + boundary + "\r\n";
    body += "Content-Disposition: form-data; name=\"QUE\"\r\n";
    body += "\r\n";
    body += string(query_spectra_.c_str());
    body += "\r\n";
    body += "--" + boundary + "--\r\n";

#ifdef MASCOTREMOTEQUERY_DEBUG_FULL_QUERY
    cerr << ">>>> Query (begin):" << "\n"
         << body << "\n"
         << "<<<< Query (end)." << endl;
#endif

    string content_type = "multipart/form-data; boundary=" + boundary;
    string response = httpPost(curl, server_path_ + "/cgi/nph-mascot.exe?1", body, content_type);

    if (hasError()) return;

    // Check for "Click here to see Search Report"
    if (response.find("Click here to see Search Report") == string::npos)
    {
      // Check for Mascot error codes
      regex mascot_error_regex(R"(\[M[0-9]{5}\])");
      smatch match;
      if (regex_search(response, match, mascot_error_regex))
      {
        OPENMS_LOG_ERROR << "Received response with Mascot error message!" << std::endl;
        if (match[0].str() == "[M00380]")
        {
          error_message_ = "You must enter an email address and user name when using the Matrix Science public web site [M00380].";
        }
        else
        {
          OPENMS_LOG_ERROR << "Error code: " << match[0].str() << std::endl;
          error_message_ = response;
        }
      }
      else
      {
        error_message_ = "Error: Unexpected response from Mascot server after query submission";
      }
      return;
    }

    // Extract .dat file path
    regex rx(R"(file=(.+/\d+/\w+\.dat))");
    smatch dat_match;
    if (!regex_search(response, dat_match, rx))
    {
      error_message_ = "Error: Could not extract .dat file path from Mascot response";
      return;
    }
    string dat_file_path = dat_match[1].str();
    search_identifier_ = getSearchIdentifierFromFilePath(dat_file_path);

    if (param_.exists("skip_export") && (param_.getValue("skip_export") == "true"))
    {
      return;
    }

    // Build export URL path
    string required_params = "&do_export=1&export_format=XML&generate_file=1&group_family=1&peptide_master=1&protein_master=1&search_master=1&show_unassigned=1&show_mods=1&show_header=1&show_params=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&pep_expect=1&pep_var_mod=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1";
    string adjustable_params = string(param_.getValue("export_params").toString());

    string results_path = string(server_path_.c_str()) + "/cgi/export_dat_2.pl?file=" + dat_file_path + required_params + "&" + adjustable_params;

    // Get results (with continuation link handling for Mascot 2.4)
    string xml_response = getResults(curl, results_path);
    if (hasError()) return;

    // Check if this is a decoy response
    if (xml_response.find("<decoy>") != string::npos)
    {
      mascot_decoy_xml_ = std::move(xml_response);
      return;
    }

    mascot_xml_ = std::move(xml_response);

    // Now retrieve decoys if requested
    if (export_decoys_)
    {
      string decoy_path = string(server_path_.c_str()) + "/cgi/export_dat_2.pl?file=" + dat_file_path + required_params + "&" + adjustable_params + "&_show_decoy_report=1&show_decoy=1";

      string decoy_response = getResults(curl, decoy_path);
      if (hasError()) return;
      mascot_decoy_xml_ = std::move(decoy_response);
    }
  }

  string MascotRemoteQuery::getResults(CURL* curl, const string& results_path)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::getResults()" << "\n";
    cerr << " from path " << results_path << "\n";
#endif

    string response = httpGet(curl, results_path);
    if (hasError()) return {};

    // Handle Mascot 2.4 continuation links (these are in-body HTML, NOT HTTP redirects)
    while (response.find("Finished after") != string::npos &&
           response.find("<a id=\"continuation-link\"") != string::npos)
    {
      regex rx(R"xx(<a id="continuation-link" href="(.*?)")xx");
      smatch match;
      if (!regex_search(response, match, rx))
      {
        break;
      }
      string continuation_path = removeHostName(match[1].str());

#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Following continuation link: " << continuation_path << "\n";
#endif

      response = httpGet(curl, continuation_path);
      if (hasError()) return {};
    }

    // Check for empty response
    if (response.empty())
    {
      error_message_ = "Error: Reply from mascot server is empty! Possible server overload - see the Mascot Admin!";
      return {};
    }

    // Check for Mascot error codes in final response
    regex mascot_error_regex(R"(\[M[0-9]{5}\])");
    smatch error_match;
    if (regex_search(response, error_match, mascot_error_regex))
    {
      OPENMS_LOG_ERROR << "Received response with Mascot error message!" << std::endl;
      if (error_match[0].str() == "[M00380]")
      {
        error_message_ = "You must enter an email address and user name when using the Matrix Science public web site [M00380].";
      }
      else
      {
        OPENMS_LOG_ERROR << "Error code: " << error_match[0].str() << std::endl;
        error_message_ = response;
      }
      return {};
    }

    return response;
  }

  string MascotRemoteQuery::removeHostName(const string& url) const
  {
    string result = url;
    if (result.substr(0, 7) == "http://")
      result = result.substr(7);
    else if (result.substr(0, 8) == "https://")
      result = result.substr(8);

    string hostname = string(host_name_.c_str());
    if (result.substr(0, hostname.size()) == hostname)
    {
      result = result.substr(hostname.size());
    }

    // Check for port suffix
    Int port = param_.getValue("host_port");
    string port_suffix = ":" + to_string(port);
    if (result.substr(0, port_suffix.size()) == port_suffix)
    {
      result = result.substr(port_suffix.size());
    }

    // ensure path starts with /
    if (result.empty() || result[0] != '/') result = "/" + result;

    return result;
  }

  void MascotRemoteQuery::setQuerySpectra(const String& exp)
  {
    query_spectra_ = exp;
  }

  const string& MascotRemoteQuery::getMascotXMLResponse() const
  {
    return mascot_xml_;
  }

  const string& MascotRemoteQuery::getMascotXMLDecoyResponse() const
  {
    return mascot_decoy_xml_;
  }

  void MascotRemoteQuery::setExportDecoys(const bool b)
  {
    export_decoys_ = b;
  }

  bool MascotRemoteQuery::hasError() const
  {
    return !error_message_.empty();
  }

  const String& MascotRemoteQuery::getErrorMessage() const
  {
    return error_message_;
  }

  String MascotRemoteQuery::getSearchIdentifier() const
  {
    return search_identifier_;
  }

  void MascotRemoteQuery::updateMembers_()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "MascotRemoteQuery::updateMembers_()" << "\n";
#endif
    server_path_ = param_.getValue("server_path").toString();
    if (!server_path_.empty())
    {
      server_path_ = "/" + server_path_;
    }
    host_name_ = param_.getValue("hostname").toString();

    use_ssl_ = param_.getValue("use_ssl").toBool();

    boundary_ = param_.getValue("boundary").toString();
    mascot_xml_.clear();
    mascot_decoy_xml_.clear();
    error_message_.clear();
    search_identifier_.clear();

    to_ = param_.getValue("timeout");
    requires_login_ = param_.getValue("login").toBool();
  }

  String MascotRemoteQuery::getSearchIdentifierFromFilePath(const String& path) const
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    std::cerr << "MascotRemoteQuery::getSearchIdentifierFromFilePath " << path << std::endl;
#endif
    size_t pos = path.find_last_of("/\\");
    String tmp = path.substr(pos + 1);
    pos = tmp.find_last_of(".");
    tmp = tmp.substr(1, pos - 1);
#ifdef MASCOTREMOTEQUERY_DEBUG
    std::cerr << "MascotRemoteQuery::getSearchIdentifierFromFilePath will return" << tmp << std::endl;
#endif
    return tmp;
  }
}
