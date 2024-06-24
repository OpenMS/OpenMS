// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QtGui/QTextDocument>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QNetworkProxy>
#include <QtNetwork/QSslSocket>

#include <QRegularExpression>

// #define MASCOTREMOTEQUERY_DEBUG
// #define MASCOTREMOTEQUERY_DEBUG_FULL_QUERY

using namespace std;

namespace OpenMS
{

  MascotRemoteQuery::MascotRemoteQuery(QObject* parent) :
    QObject(parent),
    DefaultParamHandler("MascotRemoteQuery"),
    manager_(nullptr)
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
    defaults_.setValue("skip_export", "false", "For use with an external Mascot Percolator (via GenericWrapper): Run the Mascot search, but do not export the results. The output file produced by MascotAdapterOnline will contain only the Mascot search number.", {"advanced"});
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
    if (manager_) {delete manager_;}
  }

  void MascotRemoteQuery::timedOut() const
  {
    OPENMS_LOG_FATAL_ERROR << "Mascot request timed out after " << to_ << " seconds! See 'timeout' parameter for details!" << std::endl;
  }

  void MascotRemoteQuery::run()
  {
    // Due to the asynchronous nature of Qt network requests (and the resulting use
    // of signals and slots), the information flow in this class is not very
    // clear. After the initial call to "run", the steps are roughly as follows:
    //
    // 1. optional: log into Mascot server (function "login")
    // 2. send query (function "execQuery")
    // 3. read query result, prepare exporting (function "readResponse")
    // 4. send export request (function "getResults")
    // 5. Mascot 2.4: read result, check for redirect (function "readResponse")
    // 6. Mascot 2.4: request redirected (caching) page (function "getResults")
    // (5. and 6. can happen multiple times - keep following redirects)
    // 7. Mascot 2.4: read result, check if caching is done (function "readResponse")
    // 8. Mascot 2.4: request results again (function "getResults")
    // 9. read results, which should now contain the XML (function "readResponse")
    //

    updateMembers_();

    // Make sure we do not mess with the asynchronous nature of the call and
    // start a second one while the first one is still running.
    if (manager_)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Error: Please call run() only once per MascotRemoteQuery.");
    }
    manager_ = new QNetworkAccessManager(this);

    if (!use_ssl_)
    {
      manager_->connectToHost(host_name_.c_str(), (UInt)param_.getValue("host_port"));
    }
    else
    {
#ifndef QT_NO_SSL
      manager_->connectToHostEncrypted(host_name_.c_str(), (UInt)param_.getValue("host_port"));
#else
      // should not happen since it is checked during parameter reading. Kept for safety.
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Error: Usage of SSL encryption requested but the linked QT library was not compiled with SSL support. Please recompile QT.");
#endif
    }

    connect(this, SIGNAL(gotRedirect(QNetworkReply *)), this, SLOT(followRedirect(QNetworkReply *)));
    connect(&timeout_, SIGNAL(timeout()), this, SLOT(timedOut()));
    connect(manager_, SIGNAL(finished(QNetworkReply*)), this, SLOT(readResponse(QNetworkReply*)));

    if (param_.getValue("login").toBool())
    {
      login();
    }
    else
    {
      execQuery();
    }
  }

  void MascotRemoteQuery::login()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::login()" << "\n";
#endif

    QUrl url = buildUrl_(server_path_ + "/cgi/login.pl");
    QNetworkRequest request(url);

    QByteArray boundary = boundary_.toQString().toUtf8();
    request.setHeader(QNetworkRequest::ContentTypeHeader, "multipart/form-data, boundary=" + boundary);

    // header
    request.setRawHeader("Host", host_name_.c_str());
    request.setRawHeader("Cache-Control", "no-cache");
    request.setRawHeader("Accept", "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");

    // content
    QByteArray loginbytes;
    QByteArray boundary_string("--" + boundary + "\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"username\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(String(param_.getValue("username").toString()).toQString().toUtf8());
    loginbytes.append("\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"password\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(String(param_.getValue("password").toString()).toQString().toUtf8());
    loginbytes.append("\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"submit\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("Login\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"referer\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"display\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("nothing\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"savecookie\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("1\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"action\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("login\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"userid\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"onerrdisplay\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append("login_prompt\r\n");
    loginbytes.append("--" + boundary + "--\r\n");

    request.setHeader(QNetworkRequest::ContentLengthHeader, loginbytes.length());
    QNetworkReply* pReply = manager_->post(request, loginbytes);
    connect(pReply, SIGNAL(uploadProgress(qint64, qint64)), this, SLOT(uploadProgress(qint64, qint64)));

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(request, "send");
    cerr << "  login sent reply " << pReply << std::endl;
#endif

  }

  void MascotRemoteQuery::getResults(const QString& results_path)
  {
    // Tidy up again and run another request...
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::getResults()" << "\n";
    cerr << " from path " << results_path.toStdString() << "\n";
#endif

    QUrl url = buildUrl_(results_path.toStdString());
    QNetworkRequest request(url);

#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << " server_path_ " << server_path_ << std::endl;
    cerr << " URL : " << url.toDisplayString().toStdString() << std::endl;
#endif

    // header
    request.setRawHeader("Host", host_name_.c_str());
    // request.setRawHeader("Host", String("http://www.chuchitisch.ch").c_str());
    request.setRawHeader("Accept", "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
    request.setRawHeader("Keep-Alive", "300");
    request.setRawHeader("Connection", "keep-alive");
    if (cookie_ != "")
    {
      request.setRawHeader(QByteArray::fromStdString("Cookie"), QByteArray::fromStdString(cookie_.toStdString()));
    }

    QNetworkReply* pReply = manager_->get(request);
    connect(pReply, SIGNAL(uploadProgress(qint64, qint64)), this, SLOT(uploadProgress(qint64, qint64)));

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(request, "request results");
    cerr << "  getResults expects reply " << pReply << std::endl;
#endif

  }

  void MascotRemoteQuery::followRedirect(QNetworkReply * r)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::followRedirect(QNetworkReply *)" << "\n";
    logHeader_(r, "follow redirect");

    // print full response
#ifdef MASCOTREMOTEQUERY_DEBUG_FULL_QUERY
    QByteArray new_bytes = r->readAll();
    QString str = QString::fromUtf8(new_bytes.data(), new_bytes.size());

    cerr << "Response of query: " << "\n";
    cerr << "-----------------------------------" << "\n";
    cerr << QString(new_bytes.constData()).toStdString() << "\n";
    cerr << "-----------------------------------" << "\n";
    cerr << str.toStdString() << "\n";
#endif
#endif

    // parse the location mascot wants us to go
    QString location = r->header(QNetworkRequest::LocationHeader).toString();
    removeHostName_(location);

    QUrl url = buildUrl_(location.toStdString());
    QNetworkRequest request(url);
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "Location: " << location.toStdString() << "\n";
#endif

    // header
    request.setRawHeader("Host", host_name_.c_str());
    request.setRawHeader("Accept", "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
    request.setRawHeader("Keep-Alive", "300");
    request.setRawHeader("Connection", "keep-alive");
    if (cookie_ != "")
    {
      request.setRawHeader(QByteArray::fromStdString("Cookie"), QByteArray::fromStdString(cookie_.toStdString()));
    }

#ifdef MASCOTREMOTEQUERY_DEBUG
    QNetworkReply* pReply = manager_->get(request);
    logHeader_(request, "following redirect");
    cerr << "  followRedirect expects reply " << pReply << std::endl;
#else
    manager_->get(request);
#endif
  }

  void MascotRemoteQuery::execQuery()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::execQuery()" << "\n";
#endif

    QUrl url = buildUrl_(server_path_ + "/cgi/nph-mascot.exe?1");
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << " server_path_ " << server_path_ << std::endl;
    cerr << " URL : " << url.toDisplayString().toStdString() << std::endl;
#endif
    QNetworkRequest request(url);

    QByteArray boundary = boundary_.toQString().toUtf8();
    request.setHeader(QNetworkRequest::ContentTypeHeader, "multipart/form-data, boundary=" + boundary);

    // header
    request.setRawHeader("Host", host_name_.c_str());
    request.setRawHeader("Cache-Control", "no-cache");
    request.setRawHeader("Accept", "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
    request.setRawHeader("Accept", "text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*");

    if (cookie_ != "")
    {
      request.setRawHeader(QByteArray::fromStdString("Cookie"), QByteArray::fromStdString(cookie_.toStdString()));
    }

    QByteArray querybytes;
    querybytes.append("--" + boundary + "--\n");
    querybytes.append("Content-Disposition: ");
    querybytes.append("form-data; name=\"QUE\"\n");
    querybytes.append("\n");
    querybytes.append(query_spectra_.c_str());
    querybytes.append("--" + boundary + "--\n");

    querybytes.replace("\n", "\r\n");

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(request, "request");
#ifdef MASCOTREMOTEQUERY_DEBUG_FULL_QUERY
    cerr << ">>>> Query (begin):" << "\n"
         << querybytes.constData() << "\n"
         << "<<<< Query (end)." << endl;
#endif
#endif

    if (to_ > 0)
    {
      timeout_.start();
    }

    request.setHeader(QNetworkRequest::ContentLengthHeader, querybytes.length());
    QNetworkReply* pReply = manager_->post(request, querybytes);
    connect(pReply, SIGNAL(uploadProgress(qint64, qint64)), this, SLOT(uploadProgress(qint64, qint64)));

#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "  execQuery expects reply " << pReply << std::endl;
#endif
  }

#ifdef MASCOTREMOTEQUERY_DEBUG
  void MascotRemoteQuery::downloadProgress(qint64 bytes_read, qint64 bytes_total)
#else
  void MascotRemoteQuery::downloadProgress(qint64 /*bytes_read*/, qint64 /*bytes_total*/)
#endif
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::downloadProgress(): " << bytes_read << " bytes of " << bytes_total << " read." << "\n";
    /*
    if (http_->state() == QHttp::Reading)
    {
      QByteArray new_bytes = http_->readAll();
      cerr << "Current response: " << "\n";
      cerr << QString(new_bytes.constData()).toStdString() << "\n";
    }
     */
#endif
  }

#ifdef MASCOTREMOTEQUERY_DEBUG
  void MascotRemoteQuery::uploadProgress(qint64 bytes_read, qint64 bytes_total)
#else
  void MascotRemoteQuery::uploadProgress(qint64 /* bytes_read */, qint64 /* bytes_total */)
#endif
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::uploadProgress(): " << bytes_read << " bytes of " << bytes_total << " sent." << "\n";
#endif
  }

  void MascotRemoteQuery::readResponseHeader(const QNetworkReply* reply)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::readResponseHeader(const QNetworkReply* reply)" << "\n";
    logHeader_(reply, "read");
#endif

    int status = reply->attribute( QNetworkRequest::HttpStatusCodeAttribute ).toInt();

    if (status >= 400)
    {
      error_message_ = String("MascotRemoteQuery: The server returned an error status code '") + status + "': " + reply->attribute( QNetworkRequest::HttpReasonPhraseAttribute ).toString() + "\nTry accessing the server\n  " + host_name_ + server_path_ + "\n from your browser and check if it works fine.";
      endRun_();
    }

    // Get session and username and so on...
    if (reply->header(QNetworkRequest::SetCookieHeader).isValid())
    {
      // QVariant qcookie_ = reply->header(QNetworkRequest::SetCookieHeader);

      QByteArray tmp = QByteArray::fromStdString(String("Set-Cookie"));
      QString response = reply->rawHeader(tmp);

      QRegularExpression rx("MASCOT_SESSION=(\\w+);\\spath");
      QString mascot_session(rx.match(response).captured(1));

      rx.setPattern("MASCOT_USERNAME=(\\w+);\\spath");
      QString mascot_username(rx.match(response).captured(1));

      rx.setPattern("MASCOT_USERID=(\\d+);\\spath");
      QString mascot_user_ID(rx.match(response).captured(1));

      //Put the cookie together...
      cookie_ = "userName=; userEmail=; MASCOT_SESSION=";
      cookie_.append(mascot_session);
      cookie_.append("; MASCOT_USERNAME=");
      cookie_.append(mascot_username);
      cookie_.append("; MASCOT_USERID=");
      cookie_.append(mascot_user_ID);

#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Cookie created: " << cookie_.toStdString() << "\n";
      cerr << "Cookie created from: " << response.toStdString() << "\n";
#endif
    }
  }

  void MascotRemoteQuery::endRun_()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "void MascotRemoteQuery::endRun_()" << std::endl;
#endif

    // for consumers of this class
    emit done();
  }

  void MascotRemoteQuery::readResponse(QNetworkReply* reply)
  {
    static QString dat_file_path_;
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "\nvoid MascotRemoteQuery::readResponse(const QNetworkReply* reply): ";
    if (reply->error())
    {
      cerr << "'" << reply->errorString().toStdString() << "'" << "\n";
    }
    else
    {
      cerr << "\n";
    }
    cerr << "  -- got query " << reply << endl;
#endif

    readResponseHeader(reply); // first read and parse the header (also log)

    timeout_.stop();

    if (reply->error())
    {
      error_message_ = String("Mascot Server replied: '") + String(reply->errorString().toStdString()) + "'";
      cerr << "   ending run with " + String("Mascot Server replied: '") + String(reply->errorString().toStdString()) + "'\n";
      endRun_();
      return;
    }

    QByteArray new_bytes = reply->readAll();
#ifdef MASCOTREMOTEQUERY_DEBUG_FULL_QUERY
    QString str = QString::fromUtf8(new_bytes.data(), new_bytes.size());

    cerr << "Response of query: " << "\n";
    cerr << "-----------------------------------" << "\n";
    cerr << QString(new_bytes.constData()).toStdString() << "\n";
    cerr << "-----------------------------------" << "\n";
    cerr << str.toStdString() << "\n";
#endif


    int status = reply->attribute( QNetworkRequest::HttpStatusCodeAttribute ).toInt();
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "HTTP status: " << status << "\n";
    cerr << "Response size: " << QString(new_bytes).trimmed().size() << "\n";
#endif

    // Catch "ghost queries": strange responses that are completely empty and
    // have a HTTP status of zero (and we never actually sent out)
    if (QString(new_bytes).trimmed().size() == 0 && status == 0)
    {
      return;
    }

    if (QString(new_bytes).trimmed().size() == 0 && status != 303) // status != 303: no redirect
    {
      error_message_ = "Error: Reply from mascot server is empty! Possible server overload - see the Mascot Admin!";
      endRun_();
      return;
    }

    // Successful login? fire off the search
    if (new_bytes.contains("Logged in successfu")) // Do not use the whole string. Currently Mascot writes 'successfully', but that might change...
    {
      OPENMS_LOG_INFO << "Login successful!" << std::endl;
      execQuery();
    }
    else if (new_bytes.contains("Error: You have entered an invalid password"))
    {
      error_message_ = "Error: You have entered an invalid password";
      endRun_();
    }
    else if (new_bytes.contains("is not a valid user"))
    {
      error_message_ = "Error: Username is not valid";
      endRun_();
    }
    else if (new_bytes.contains("Click here to see Search Report"))
    {
      //Extract filename from this...
      //e.g. data might look like
      // <A HREF="../cgi/master_results.pl?file=../data/20100728/F018032.dat">Click here to see Search Report</A>
      QString response(new_bytes);

      QRegularExpression rx(R"(file=(.+/\d+/\w+\.dat))");
      rx.setPatternOptions(QRegularExpression::InvertedGreedinessOption);
      dat_file_path_ = rx.match(response).captured(1);
      search_identifier_ = getSearchIdentifierFromFilePath(dat_file_path_);

      if (param_.exists("skip_export") &&
          (param_.getValue("skip_export") == "true"))
      {
        endRun_();
        return;
      }

      QString results_path("");
      results_path.append(server_path_.toQString());
      results_path.append("/cgi/export_dat_2.pl?file=");
      results_path.append(dat_file_path_);

#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Results path to export: " << results_path.toStdString() << "\n";
#endif

      // see http://www.matrixscience.com/help/export_help.html for parameter documentation
      String required_params = "&do_export=1&export_format=XML&generate_file=1&group_family=1&peptide_master=1&protein_master=1&search_master=1&show_unassigned=1&show_mods=1&show_header=1&show_params=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&pep_expect=1&pep_var_mod=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1";
    // TODO: check that export_params don't contain _show_decoy_report=1. This would add <decoy> at the top of the XML document and we can't easily distinguish the two files (target/decoy) from the search results later.
      String adjustable_params = param_.getValue("export_params").toString();

      results_path.append(required_params.toQString() + "&" +
                          adjustable_params.toQString());
      // results_path.append("&show_same_sets=1&show_unassigned=1&show_queries=1&do_export=1&export_format=XML&pep_rank=1&_sigthreshold=0.99&_showsubsets=1&show_header=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&show_mods=1&pep_var_mod=1&protein_master=1&search_master=1&show_params=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1&pep_expect=1&peptide_master=1&generate_file=1&group_family=1");

      // Finished search, fire off results retrieval
      getResults(results_path);
    }
    else if (status == 303)
    {
#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Retrieved redirect \n";
#endif
      emit gotRedirect(reply);
    }
    // mascot 2.4 export page done
    else if (new_bytes.contains("Finished after") &&  new_bytes.contains("<a id=\"continuation-link\""))
    {
      // parsing mascot 2.4 continuation download link
      // <p>Finished after 0 s. <a id="continuation-link" ..
      QString response(new_bytes);
      QRegularExpression rx("<a id=\"continuation-link\" href=\"(.*)\"");
      rx.setPatternOptions(QRegularExpression::InvertedGreedinessOption);
      QString results_path = rx.match(response).captured(1);

      // fill result link again and download
      removeHostName_(results_path);
      //Finished search, fire off results retrieval
      getResults(results_path);
    }
    else
    {
      // check whether Mascot responded using an error code e.g. [M00440], pipe through results else
      QString response_text = new_bytes;
      QRegularExpression mascot_error_regex(R"(\[M[0-9][0-9][0-9][0-9][0-9]\])");
      QRegularExpressionMatch match;
      if (response_text.contains(mascot_error_regex, &match))
      {
        OPENMS_LOG_ERROR << "Received response with Mascot error message!" << std::endl;
        if (match.captured() == "[M00380]")
        {
          // we know this error, so we give a much shorter and readable error message for the user
          error_message_ = "You must enter an email address and user name when using the Matrix Science public web site [M00380].";
          OPENMS_LOG_ERROR << error_message_ << std::endl;
        }
        else
        {
          OPENMS_LOG_ERROR << "Error code: " << match.captured().toStdString() << std::endl;
          error_message_ = response_text;
        }
        endRun_();
      }
      else
      {
        // seems to be fine so grab the xml
#ifdef MASCOTREMOTEQUERY_DEBUG
        cerr << "Get the XML File" << "\n";
#endif
        if (new_bytes.contains("<decoy>"))
        {
          mascot_decoy_xml_ = new_bytes;          
          endRun_();
        }
        else
        {
          mascot_xml_ = new_bytes;

          // now retrieve decos
          if (export_decoys_)
          {
            QString results_path("");
            results_path.append(server_path_.toQString());
            results_path.append("/cgi/export_dat_2.pl?file=");
            results_path.append(dat_file_path_); // export again from dat file (now decoys)        
            // see http://www.matrixscience.com/help/export_help.html for parameter documentation
            String required_params = "&do_export=1&export_format=XML&generate_file=1&group_family=1&peptide_master=1&protein_master=1&search_master=1&show_unassigned=1&show_mods=1&show_header=1&show_params=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&pep_expect=1&pep_var_mod=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1";
            // TODO: check that export_params don't contain _show_decoy_report=1. This would add <decoy> at the top of the XML document and we can't easily distinguish the two files (target/decoy) from the search results later.
            String adjustable_params = param_.getValue("export_params").toString();
            results_path.append(required_params.toQString() + "&" +
                            adjustable_params.toQString() + "&_show_decoy_report=1&show_decoy=1"); // request 

            getResults(results_path); 
          }
          else
          {
            endRun_();
          }
        }
      }
    }
  }

  void MascotRemoteQuery::setQuerySpectra(const String& exp)
  {
    query_spectra_ = exp;
  }

  const QByteArray& MascotRemoteQuery::getMascotXMLResponse() const
  {
    return mascot_xml_;
  }

  const QByteArray& MascotRemoteQuery::getMascotXMLDecoyResponse() const
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
    //MascotRemoteQuery_test
    if (!server_path_.empty())
    {
      server_path_ = "/" + server_path_;
    }
    host_name_ = param_.getValue("hostname").toString();
    
    use_ssl_ = param_.getValue("use_ssl").toBool();
#ifndef QT_NO_SSL
    if (use_ssl_ && !QSslSocket::supportsSsl())
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Error: Usage of SSL encryption requested but the OpenSSL library was not found at runtime. Please install OpenSSL system-wide.");
    }
#else
    if (use_ssl_)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Error: Usage of SSL encryption requested but the linked QT library was not compiled with SSL support. Please recompile QT.");
    }
#endif
    
    boundary_ = param_.getValue("boundary").toString();
    cookie_ = "";
    mascot_xml_ = "";
    mascot_decoy_xml_ = "";

    to_ = param_.getValue("timeout");
    timeout_.setInterval(1000 * to_);
    requires_login_ = param_.getValue("login").toBool();

    bool use_proxy(param_.getValue("use_proxy").toBool());
    if (use_proxy)
    {
      QNetworkProxy proxy;
      proxy.setType(QNetworkProxy::Socks5Proxy);
      String proxy_host(param_.getValue("proxy_host").toString());
      proxy.setHostName(proxy_host.toQString());
      String proxy_port(param_.getValue("proxy_port").toString());
      proxy.setPort(proxy_port.toInt());

      String proxy_password(param_.getValue("proxy_password").toString());
      proxy.setPassword(proxy_password.toQString());

      String proxy_username(param_.getValue("proxy_username").toString());
      if (!proxy_username.empty())
      {
        proxy.setUser(proxy_username.toQString());
      }

      QNetworkProxy::setApplicationProxy(proxy);
    }


    return;
  }

  void MascotRemoteQuery::removeHostName_(QString& url)
  {
    if (url.startsWith("http://"))
      url.remove("http://");
    else if (url.startsWith("https://"))
      url.remove("https://");

    if (!url.startsWith(host_name_.toQString()))
    {
      OPENMS_LOG_ERROR << "Invalid location returned by mascot! Abort." << std::endl;
      endRun_();
      return;
    }
    
    // makes sure that only the first occurrence in the String is replaced,
    // in case of an equal server_name_
    url.replace(url.indexOf(host_name_.toQString()),
                          host_name_.toQString().size(), QString(""));

    // ensure path starts with /
    if (url[0] != '/') url.prepend('/');
  }

  QUrl MascotRemoteQuery::buildUrl_(const std::string& path)
  {
    String protocol;
    if (use_ssl_)
    {
      protocol = "https";
    }
    else
    {
      protocol = "http";
    }
    return QUrl(String(protocol + "://" + host_name_ + path).c_str());
  }

  void MascotRemoteQuery::logHeader_(const QNetworkReply* header, const String& what)
  {
    QList<QByteArray> header_list = header->rawHeaderList();
    cerr << ">>>> Header to " << what << " (begin):\n";
    foreach (QByteArray rawHeader, header_list)
    {
      cerr << "    " << rawHeader.toStdString() << " : " << header->rawHeader(rawHeader).toStdString() << std::endl;
    }
    cerr << "<<<< Header to " << what << " (end)." << endl;
  }

  void MascotRemoteQuery::logHeader_(const QNetworkRequest& header, const String& what)
  {
    QList<QByteArray> header_list = header.rawHeaderList();
    cerr << ">>>> Header to " << what << " (begin):\n";
    foreach (QByteArray rawHeader, header_list)
    {
      cerr << "    " << rawHeader.toStdString() << " : " << header.rawHeader(rawHeader).toStdString() << std::endl;
    }
    cerr << "<<<< Header to " << what << " (end)." << endl;
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

