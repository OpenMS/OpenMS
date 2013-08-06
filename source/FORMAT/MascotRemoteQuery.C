// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QtGui/QTextDocument>
#include <iostream>

// #define MASCOTREMOTEQUERY_DEBUG

using namespace std;

namespace OpenMS
{

  MascotRemoteQuery::MascotRemoteQuery(QObject* parent) :
    QObject(parent),
    DefaultParamHandler("MascotRemoteQuery")
  {
    http_ = new QHttp();

    // server specifications
    defaults_.setValue("hostname", "", "Address of the host where Mascot listens, e.g. 'mascot-server' or '127.0.0.1'");
    defaults_.setValue("host_port", 80, "Port where the Mascot server listens, 80 should be a good guess");
    defaults_.setMinInt("host_port", 0);
    defaults_.setValue("server_path", "mascot", "Path on the server where Mascot server listens, 'mascot' should be a good guess");
    defaults_.setValue("timeout", 1500, "Timeout in seconds, after which the query is declared as failed."
                                        "This is NOT the whole time the search takes, but the time in between two progress steps. Some Mascot servers freeze during this (unstable network etc) and idle forever"
                                        ", the connection is killed. Set this to 0 to disable timeout!");
    defaults_.setMinInt("timeout", 0);
    defaults_.setValue("boundary", "GZWgAaYKjHFeUaLOLEIOMq", "Boundary for the MIME section", StringList::create("advanced"));

    // proxy settings
    defaults_.setValue("use_proxy", "false", "Flag which enables the proxy usage for the http requests, please specify at least 'proxy_host' and 'proxy_port'", StringList::create("advanced"));
    defaults_.setValidStrings("use_proxy", StringList::create("true,false"));
    defaults_.setValue("proxy_host", "", "Host where the proxy server runs on", StringList::create("advanced"));
    defaults_.setValue("proxy_port", 0, "Port where the proxy server listens", StringList::create("advanced"));
    defaults_.setMinInt("proxy_port", 0);
    defaults_.setValue("proxy_username", "", "Login name for the proxy server, if needed", StringList::create("advanced"));
    defaults_.setValue("proxy_password", "", "Login password for the proxy server, if needed", StringList::create("advanced"));

    // login for Mascot security
    defaults_.setValue("login", "false", "Flag which should be set 'true' if Mascot security is enabled; also set 'username' and 'password' then.");
    defaults_.setValidStrings("login", StringList::create("true,false"));
    defaults_.setValue("username", "", "Name of the user if login is used (Mascot security must be enabled!)");
    defaults_.setValue("password", "", "Password of the user if login is used (Mascot security must be enabled!)");
    defaults_.setValue("use_ssl", "false", "Flag indicating wether the server uses https or not.");
    defaults_.setValidStrings("use_ssl", StringList::create("true,false"));

    // Mascot specific options
    defaults_.setValue("query_master", "false", "If this option is set to true, query peptides will be returned with non-truncated lists, however, protein references of peptides will not be correct.", StringList::create("advanced"));
    defaults_.setValidStrings("query_master", StringList::create("true,false"));
    defaults_.setValue("max_hits", 0, "Maximum number of hits to be exported from the server to XML (0 = AUTO)");
    defaults_.setMinInt("max_hits", 0);
    defaultsToParam_();
  }

  MascotRemoteQuery::~MascotRemoteQuery()
  {
    if (http_->state() != QHttp::Unconnected)
    {
#ifdef MASCOTREMOTEQUERY_DEBUG
      std::cerr << "Aborting open connection!\n";
#endif
      http_->abort(); // hardcore close connection (otherwise server might have too many dangling requests)
    }
    delete http_;
  }

  void MascotRemoteQuery::timedOut()
  {
    LOG_FATAL_ERROR << "Mascot request timed out after " << to_ << " seconds! See 'timeout' parameter for details!" << std::endl;
    http_->abort(); // one might try to resend the job here instead...
  }

  void MascotRemoteQuery::run()
  {
    // Due to the asynchronous nature of QHttp::request (and the resulting use
    // of signals and slots), the information flow in this class is not very
    // clear. After the initial call to "run", the steps are roughly as follows:
    // 1. optional: log into Mascot server (function "login")
    // 2. send query (function "execQuery")
    // 3. read query result, prepare exporting (function "httpDone")
    // 4. send export request (function "getResults")
    // 5. Mascot 2.4: read result, check for redirect (function "httpDone")
    // 6. Mascot 2.4: request redirected (caching) page (function "getResults")
    // (5. and 6. can happen multiple times - keep following redirects)
    // 7. Mascot 2.4: read result, check if caching's done (function "httpDone")
    // 8. Mascot 2.4: request results again (function "getResults")
    // 9. read results, which should now contain the XML (function "httpDone")

    updateMembers_();
    connect(http_, SIGNAL(requestFinished(int, bool)), this, SLOT(httpRequestFinished(int, bool)));
    connect(http_, SIGNAL(requestStarted(int)), this, SLOT(httpRequestStarted(int)));
    connect(http_, SIGNAL(done(bool)), this, SLOT(httpDone(bool)));
    connect(http_, SIGNAL(stateChanged(int)), this, SLOT(httpStateChanged(int)));
    connect(http_, SIGNAL(readyRead(const QHttpResponseHeader &)), this, SLOT(readyReadSlot(const QHttpResponseHeader &)));
    connect(http_, SIGNAL(responseHeaderReceived(const QHttpResponseHeader &)), this, SLOT(readResponseHeader(const QHttpResponseHeader &)));
    connect(this, SIGNAL(gotRedirect(const QHttpResponseHeader &)), this, SLOT(followRedirect(const QHttpResponseHeader &)));
    connect(&timeout_, SIGNAL(timeout()), this, SLOT(timedOut()));

    // get progress information
    connect(http_, SIGNAL(dataReadProgress(int, int)), this, SLOT(httpDataReadProgress(int, int)));

    if (param_.getValue("login").toBool()) login();
    else execQuery();
  }

  void MascotRemoteQuery::login()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::login()" << "\n";
#endif

    // header
    QHttpRequestHeader header;
    QString boundary = boundary_.toQString();
    header.setRequest("POST", (server_path_ + "/cgi/login.pl").c_str());
    header.setValue("Host", host_name_.c_str());
    header.setValue("Content-Type", "multipart/form-data, boundary=" + boundary);
    header.setValue("Cache-Control", "no-cache");
    header.setValue("Accept", "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");

    // content
    QByteArray loginbytes;
    QString boundary_string("--" + boundary + "\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"username\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(((String)param_.getValue("username")).c_str());
    loginbytes.append("\r\n");
    loginbytes.append(boundary_string);
    loginbytes.append("Content-Disposition: ");
    loginbytes.append("form-data; name=\"password\"\r\n");
    loginbytes.append("\r\n");
    loginbytes.append(((String)param_.getValue("password")).c_str());
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

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(header, "send");
#endif

    header.setContentLength(loginbytes.length());
    http_->request(header, loginbytes);

  }

  void MascotRemoteQuery::getResults(QString results_path)
  {
    //Tidy up again and run another request...
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::getResults()" << "\n";
#endif

    QHttpRequestHeader header;
    header.setRequest("GET", results_path);
    header.setValue("Host", host_name_.toQString());
    header.setValue("Accept", "text/xml,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
    header.setValue("Keep-Alive", "300");
    header.setValue("Connection", "keep-alive");
    if (cookie_ != "")
    {
      header.setValue("Cookie", cookie_);
    }

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(header, "request results");
#endif

    http_->request(header);
  }

  void MascotRemoteQuery::followRedirect(const QHttpResponseHeader& resp)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::followRedirect(const QHttpResponseHeader &)" << "\n";
    cerr << "Header containing redirect: \n";
    cerr << resp.toString().toStdString();
    cerr << "END HEADER" << "\n";
    cerr << "Location header: " << resp.value("Location").toStdString() << "\n";
#endif
    // parse the location mascot wants us to go
    QString location = resp.value("Location");
    removeHostName_(location);

    QHttpRequestHeader header;
    header.setRequest("GET", location);
    header.setValue("Host", host_name_.toQString());
    header.setValue("Accept", "text/xml,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
    header.setValue("Keep-Alive", "300");
    header.setValue("Connection", "keep-alive");
    if (cookie_ != "")
    {
      header.setValue("Cookie", cookie_);
    }

#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << ">>> Header to follow redirect: " << "\n";
    cerr << header.toString().toStdString() << "\n";
    cerr << "ended: " << "\n";
#endif

    http_->request(header);
  }

  void MascotRemoteQuery::execQuery()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::execQuery()" << "\n";
#endif

    QHttpRequestHeader header;
    QString boundary(boundary_.toQString());

    header.setRequest("POST", (server_path_ + "/cgi/nph-mascot.exe?1").c_str());
    header.setValue("Host", host_name_.toQString());
    header.setValue("Content-Type", ("multipart/form-data, boundary=" + boundary));
    header.setValue("Cache-Control", "no-cache");
    if (cookie_ != "")
    {
      header.setValue("Cookie", cookie_);
    }
    header.setValue("Accept", "text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*");

    QByteArray querybytes;
    querybytes.append("--" + boundary + "--\n");
    querybytes.append("Content-Disposition: ");
    querybytes.append("form-data; name=\"QUE\"\n");
    querybytes.append("\n");
    querybytes.append(query_spectra_.c_str());
    querybytes.append("--" + boundary + "--\n");

    querybytes.replace("\n", "\r\n");

    header.setContentLength(querybytes.length());

#ifdef MASCOTREMOTEQUERY_DEBUG
    logHeader_(header, "request");
    cerr << ">>>> Query (begin):" << "\n"
         << querybytes.constData() << "\n"
         << "<<<< Query (end)." << endl;
#endif

    if (to_ > 0)
      timeout_.start();
    http_->request(header, querybytes);
  }

  void MascotRemoteQuery::httpRequestFinished(int requestId, bool error)
  {
    if (error)
    {
      cerr << "MascotRemoteQuery: An error occurred (requestId=" << requestId << "): " << http_->errorString().toStdString() << " (QT Error Code: " << int(http_->error()) << ")\n";
    }
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "Request Finished Id: " << requestId << "\n";
    cerr << "Error: " << error << "(" << http_->errorString().toStdString() << ")" << "\n";
#endif

  }

#ifdef MASCOTREMOTEQUERY_DEBUG
  void MascotRemoteQuery::httpDataReadProgress(int bytes_read, int bytes_total)
#else
  void MascotRemoteQuery::httpDataReadProgress(int /*bytes_read*/, int /*bytes_total*/)
#endif
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::httpDataReadProgress(): " << bytes_read << " bytes of " << bytes_total << " read." << "\n";
    /*
    if (http_->state() == QHttp::Reading)
    {
      QByteArray new_bytes = http_->readAll();
      cerr << "Current response: " << "\n";
      QTextDocument doc;
      doc.setHtml(new_bytes.constData());
      cerr << doc.toPlainText().toStdString() << "\n";
    }
     */
#endif
  }

#ifdef MASCOTREMOTEQUERY_DEBUG
  void MascotRemoteQuery::httpDataSendProgress(int bytes_sent, int bytes_total)
#else
  void MascotRemoteQuery::httpDataSendProgress(int /*bytes_sent*/, int /*bytes_total*/)
#endif
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::httpDataSendProgress(): " << bytes_sent << " bytes of " << bytes_total << " sent." << "\n";
#endif
  }

  void MascotRemoteQuery::httpRequestStarted(int
#ifdef MASCOTREMOTEQUERY_DEBUG
                                             requestId
#endif
                                             )
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cout << "Request started: " << requestId << "\n";
#endif
  }

  void MascotRemoteQuery::httpStateChanged(int
#ifdef MASCOTREMOTEQUERY_DEBUG
                                           state
#endif
                                           )
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    switch (state)
    {
    case QHttp::Closing:
      cout << "State change to QHttp::Closing\n";
      break;

    case QHttp::Unconnected:
      cout << "State change to QHttp::Unconnected\n";
      break;

    case QHttp::HostLookup:
      cout << "State change to QHttp::HostLookup\n";
      break;

    case QHttp::Sending:
      cout << "State change to QHttp::Sending\n";
      break;

    case QHttp::Reading:
      cout << "State change to QHttp::Reading\n";
      break;

    case QHttp::Connected:
      cout << "State change to QHttp::Connected\n";
      break;

    case QHttp::Connecting:
      cout << "State change to QHttp::Connecting\n";
      break;
    }
#endif
  }

  void MascotRemoteQuery::readyReadSlot(const QHttpResponseHeader& /* resp */)
  {
    //if (http_->bytesAvailable() < 1000) std::cerr << "new bytes: " << http_->bytesAvailable() << " from " << resp.toString() << " with code " <<  resp.statusCode() << " and httpstat: " << http_->state() << "\n";
    if (to_ > 0)
      timeout_.start(); // reset timeout
  }

  void MascotRemoteQuery::readResponseHeader(const QHttpResponseHeader& response_header)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::readResponseHeader(const QHttpResponseHeader &responseHeader)" << "\n";
    logHeader_(response_header, "read");
#endif

    if (response_header.statusCode() >= 400)
    {
      error_message_ = String("MascotRemoteQuery: The server returned an error status code '") + response_header.statusCode() + "': " + response_header.reasonPhrase() + "\nTry accessing the server\n  " + host_name_ + server_path_ + "\n from your browser and check if it works fine.";
      endRun_();
    }


    //Get session and username and so on...
    if (response_header.hasKey("Set-Cookie"))
    {
      QString response = response_header.toString();

      QRegExp rx("MASCOT_SESSION=(\\w+);\\spath");
      rx.indexIn(response);
      QString mascot_session(rx.cap(1));

      rx.setPattern("MASCOT_USERNAME=(\\w+);\\spath");
      rx.indexIn(response);
      QString mascot_username(rx.cap(1));

      rx.setPattern("MASCOT_USERID=(\\d+);\\spath");
      rx.indexIn(response);
      QString mascot_user_ID(rx.cap(1));

      //Put the cookie together...
      cookie_ = "userName=; userEmail=; MASCOT_SESSION=";
      cookie_.append(mascot_session);
      cookie_.append("; MASCOT_USERNAME=");
      cookie_.append(mascot_username);
      cookie_.append("; MASCOT_USERID=");
      cookie_.append(mascot_user_ID);

#ifdef MASCOTREMOTEQUERY_DEBUG
      cout << "Cookie created:" << cookie_.toStdString() << "\n";
#endif
    }
  }

  void MascotRemoteQuery::endRun_()
  {
    if (http_->state() != QHttp::Unconnected)
    {
      http_->clearPendingRequests();
      http_->close();
    }
    emit done();
  }

  void MascotRemoteQuery::httpDone(bool error)
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "void MascotRemoteQuery::httpDone(bool error): ";
    if (error)
    {
      cerr << "'" << http_->errorString().toStdString() << "'" << "\n";
    }
    else
    {
      cerr << "\n";
    }
#endif

    timeout_.stop();

    if (error)
    {
      error_message_ = String("Mascot Server replied: '") + String(http_->errorString().toStdString()) + "'";
      endRun_();
      return;
    }

    QByteArray new_bytes = http_->readAll();
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "Response of query: " << "\n";
    QTextDocument doc;
    doc.setHtml(new_bytes.constData());
    cerr << doc.toPlainText().toStdString() << "\n";
#endif

    if (QString(new_bytes).trimmed().size() == 0 && !(http_->lastResponse().isValid() && http_->lastResponse().statusCode() == 303))
    {
      error_message_ = "Error: Reply from mascot server is empty! Possible server overload - see the Mascot Admin!";
      endRun_();
      return;
    }

    //Successful login? fire off the search
    if (new_bytes.contains("Logged in successfu")) // Do not use the whole string. Currently Mascot writes 'successfuly', but that might change...
    {
      //Successful login? fire off the search
      LOG_INFO << "Login successful!" << std::endl;
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

      QRegExp rx("file=(.+/\\d+/\\w+\\.dat)");
      rx.setMinimal(true);
      rx.indexIn(response);

      QString results_path("");
      results_path.append(server_path_.toQString());
      results_path.append("/cgi/export_dat_2.pl?file=");
      results_path.append(rx.cap(1));
      results_path.append(QString("&report=") + QString::number(max_hits_));

#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Results path to export: " << results_path.toStdString() << "\n";
#endif
      //results_path_.append("&show_same_sets=1&show_unassigned=1&show_queries=1&do_export=1&export_format=XML&pep_rank=1&_sigthreshold=0.99&_showsubsets=1&show_header=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&show_mods=1&pep_var_mod=1&protein_master=1&prot_score=1&search_master=1&show_header=1&show_params=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1&pep_expect=1&peptide_master=1"); // contains duplicate options
      results_path.append("&show_same_sets=1&show_unassigned=1&show_queries=1&do_export=1&export_format=XML&pep_rank=1&_sigthreshold=0.99&_showsubsets=1&show_header=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&show_mods=1&pep_var_mod=1&protein_master=1&search_master=1&show_params=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1&pep_expect=1&peptide_master=1&generate_file=1&group_family=1");

      if (param_.getValue("query_master").toBool())
      {
        results_path.append("&query_master=1");
      }
      else
      {
        results_path.append("&query_master=0");
      }
      //Finished search, fire off results retrieval
      getResults(results_path);
    }
    else if (http_->lastResponse().statusCode() == 303)
    {
#ifdef MASCOTREMOTEQUERY_DEBUG
      cerr << "Retrieved redirect \n";
#endif
      emit gotRedirect(http_->lastResponse());
    }
    // mascot 2.4 export page done
    else if (new_bytes.contains("Finished after") &&  new_bytes.contains("<a id=\"continuation-link\""))
    {
      // parsing mascot 2.4 continuation download link
      // <p>Finished after 0 s. <a id="continuation-link" ..
      QString response(new_bytes);
      QRegExp rx("<a id=\"continuation-link\" href=\"(.*)\"");
      rx.setMinimal(true);
      rx.indexIn(response);
      QString results_path = rx.cap(1);

      // fill result link again and download
      removeHostName_(results_path);
      //Finished search, fire off results retrieval
      getResults(results_path);
    }
    else
    {
      // check whether Mascot responded using an error code e.g. [M00440], pipe through results else
      QString response_text = new_bytes;
      QRegExp mascot_error_regex("\\[M[0-9][0-9][0-9][0-9][0-9]\\]");
      if (response_text.contains(mascot_error_regex))
      {
        LOG_ERROR << "Received response with Mascot error message!" << std::endl;
        if (mascot_error_regex.cap() == "[M00380]")
        {
          // we know this error, so we give a much shorter and readable error message for the user
          LOG_ERROR << "You must enter an email address and user name when using the Matrix Science public web site [M00380]." << std::endl;
          error_message_ = "You must enter an email address and user name when using the Matrix Science public web site [M00380].";
        }
        else
        {
          LOG_ERROR << "Error code: " << mascot_error_regex.cap().toStdString() << std::endl;
          QTextDocument doc;
          doc.setHtml(response_text);
          error_message_ = doc.toPlainText().toStdString();
        }
        endRun_();
      }
      else
      {
        // seems to be fine so grab the xml
#ifdef MASCOTREMOTEQUERY_DEBUG
        cerr << "Get the XML File" << "\n";
#endif
        mascot_xml_ = new_bytes;
        endRun_();
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

  bool MascotRemoteQuery::hasError() const
  {
    return error_message_ != "";
  }

  const String& MascotRemoteQuery::getErrorMessage() const
  {
    return error_message_;
  }

  void MascotRemoteQuery::updateMembers_()
  {
#ifdef MASCOTREMOTEQUERY_DEBUG
    cerr << "MascotRemoteQuery::updateMembers_()" << "\n";
#endif
    server_path_ = param_.getValue("server_path");
    //MascotRemoteQuery_test
    if (server_path_ != "")
      server_path_ = "/" + server_path_;

    host_name_ = param_.getValue("hostname");
    use_ssl_ = param_.getValue("use_ssl").toBool();

    // clear all content from this class
    delete http_;
    http_ = new QHttp(this);
    http_->setHost(host_name_.c_str(),
                   (use_ssl_ ? QHttp::ConnectionModeHttps : QHttp::ConnectionModeHttp),
                   (UInt)param_.getValue("host_port"));

    boundary_ = param_.getValue("boundary");
    cookie_ = "";
    mascot_xml_ = "";

    to_ = param_.getValue("timeout");
    timeout_.setInterval(1000 * to_);
    requires_login_ = param_.getValue("login").toBool();
    max_hits_ = param_.getValue("max_hits");

    bool use_proxy(param_.getValue("use_proxy").toBool());
    if (use_proxy)
    {
      String proxy_host(param_.getValue("proxy_host"));
      String proxy_port(param_.getValue("proxy_port"));
      String proxy_username(param_.getValue("proxy_username"));
      String proxy_password(param_.getValue("proxy_password"));
      if (proxy_username != "")
      {
        http_->setProxy(proxy_host.c_str(), proxy_port.toInt(), proxy_username.c_str(), proxy_password.c_str());
      }
      else
      {
        http_->setProxy(proxy_host.c_str(), proxy_port.toInt());
      }
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
      LOG_ERROR << "Invalid location returned by mascot! Abort." << std::endl;
      endRun_();
      return;
    }
    url.remove(host_name_.toQString());

    // ensure path starts with /
    if (url[0] != '/') url.prepend('/');
  }

  void MascotRemoteQuery::logHeader_(const QHttpHeader& header,
                                     const String& what)
  {
    cerr << ">>>> Header to " << what << " (begin):\n"
         << header.toString().toStdString()
         << "<<<< Header to " << what << " (end)." << endl;
  }

}
