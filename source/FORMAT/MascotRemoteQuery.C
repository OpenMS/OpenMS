// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QtGui/QTextDocument>
#include <iostream>

#define MASCOTREMOTEQUERY_DEBUG
#undef  MASCOTREMOTEQUERY_DEBUG

using namespace std;

namespace OpenMS
{

	MascotRemoteQuery::MascotRemoteQuery(QObject* parent)
		:	QObject(parent),
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
		defaults_.setValue("use_proxy", "false", "Flag which enabled the proxy usage for the http requests, please specify at least 'proxy_host' and 'proxy_port'", StringList::create("advanced"));
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

		// Mascot specific options
		defaults_.setValue("query_master", "false", "If this option is set to true, query peptides will be returned with non-truncated lists, however, protein references of peptides will not be correct.", StringList::create("advanced"));
		defaults_.setValidStrings("query_master", StringList::create("true,false"));

		defaultsToParam_();
	}

	MascotRemoteQuery::~MascotRemoteQuery()
	{
    if (http_->state() != QHttp::Unconnected) 
    {
      std::cerr << "Aborting open connection!\n";
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
		updateMembers_();
    connect(http_, SIGNAL(requestFinished(int, bool)), this, SLOT(httpRequestFinished(int,bool)));
    connect(http_, SIGNAL(requestStarted(int)), this, SLOT(httpRequestStarted(int)));
    connect(http_, SIGNAL(done(bool)), this, SLOT(httpDone(bool)));
    connect(http_, SIGNAL(stateChanged(int)), this, SLOT(httpStateChanged(int)));
    connect(http_, SIGNAL(readyRead ( const QHttpResponseHeader & )), this, SLOT(readyReadSlot ( const QHttpResponseHeader & )) );
    connect(http_, SIGNAL(responseHeaderReceived(const QHttpResponseHeader &)), this, SLOT(readResponseHeader(const QHttpResponseHeader &)));
    connect(this, SIGNAL(loginDone()), this, SLOT(execQuery()));
    connect(this, SIGNAL(queryDone()), this, SLOT(getResults()));
    connect(&timeout_, SIGNAL(timeout()), this, SLOT(timedOut()));

    if (param_.getValue("login").toBool())
		{
			login();
		}
		else
		{
			emit loginDone();
		}
		return;
	}


void MascotRemoteQuery::login() 
{
#ifdef MASCOTREMOTEQUERY_DEBUG	
	cerr << "void MascotRemoteQuery::login()" << "\n";
#endif

	// header
	QHttpRequestHeader header;
	QString boundary(((String)param_.getValue("boundary")).c_str());
	header.setRequest("POST", ("/" + (String)param_.getValue("server_path") + "/cgi/login.pl").c_str());
	header.setValue("Host", ((String)param_.getValue("hostname")).c_str());
	header.setValue("Content-Type", "multipart/form-data, boundary=" + boundary);
 	header.setValue("Cache-Control", "no-cache");
	header.setValue("Accept","text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");	

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
	cerr << ">>>> Header to send: " << "\n";
	cerr << header.toString().toStdString() << "\n";
	cerr << "ended" << "\n";
#endif

	header.setContentLength(loginbytes.length());
 	http_->request(header, loginbytes);

}

void MascotRemoteQuery::getResults() 
{
	//Tidy up again and run another request...
#ifdef MASCOTREMOTEQUERY_DEBUG	
	cerr << "void MascotRemoteQuery::getResults()" << "\n";
#endif
	
	QHttpRequestHeader header;
	header.setRequest("GET", results_path_);
	header.setValue("Host", ((String)param_.getValue("hostname")).c_str());
	header.setValue("Accept","text/xml,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
	header.setValue("Keep-Alive","300");
	header.setValue("Connection","keep-alive");
	if (cookie_ != "")
	{
		header.setValue("Cookie", cookie_);
	}

#ifdef MASCOTREMOTEQUERY_DEBUG
	cerr << ">>> Header to request results: " << "\n";
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
	
	header.setRequest("POST", ("/" + (String)param_.getValue("server_path") + "/cgi/nph-mascot.exe?1").c_str());
	header.setValue("Host", ((String)param_.getValue("hostname")).c_str());
 	header.setValue("Content-Type", ("multipart/form-data, boundary=" + (String)param_.getValue("boundary")).c_str());
 	header.setValue("Cache-Control", "no-cache");
	if (cookie_ != "")
	{
		header.setValue("Cookie", cookie_);	
	}
	header.setValue("Accept","text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*");
 	
	QByteArray querybytes;
	querybytes.append(query_spectra_.c_str());
	querybytes.replace("\n", "\r\n");
 	header.setContentLength(querybytes.length());
#ifdef MASCOTREMOTEQUERY_DEBUG
	cerr << ">>>> Header to request:" << "\n";
	cerr << header.toString().toStdString() << "\n";
	cerr << "ended: " << "\n";

	cerr << ">>>> Query:" << "\n";
	cerr << querybytes.constData() << "\n";
	cerr << "ended: " << "\n";
#endif
  
  if (to_ > 0) timeout_.start();
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
	cout<<"Request started: "<<requestId<<"\n";
#endif
}

void MascotRemoteQuery::httpStateChanged(int 
#ifdef MASCOTREMOTEQUERY_DEBUG
state
#endif
)
{
#ifdef MASCOTREMOTEQUERY_DEBUG
	cout<<"State change: "<<state<<"\n";
#endif
}

void MascotRemoteQuery::readyReadSlot ( const QHttpResponseHeader & /* resp */)
{
  //if (http_->bytesAvailable() < 1000) std::cerr << "new bytes: " << http_->bytesAvailable() << " from " << resp.toString() << " with code " <<  resp.statusCode() << " and httpstat: " << http_->state() << "\n";
  if (to_ > 0) timeout_.start(); // reset timeout
}


void MascotRemoteQuery::readResponseHeader(const QHttpResponseHeader& response_header) 
{
#ifdef MASCOTREMOTEQUERY_DEBUG
	cerr << "void MascotRemoteQuery::readResponseHeader(const QHttpResponseHeader &responseHeader)" << "\n";

	cerr << ">>>>> Header to read: " << "\n";
	cerr <<  response_header.toString().toStdString() << "\n";
	cerr << "ended" << "\n";
#endif

	if (response_header.statusCode() >= 400)
	{
		error_message_ = String("MascotRemoteQuery: The server returned an error status code '") + response_header.statusCode() + "': " + response_header.reasonPhrase() + "\nTry accessing the server\n  " + (String)param_.getValue("hostname") + "/" + (String)param_.getValue("server_path") + "\n from your browser and check if it works fine.";
    endRun();
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
cout<<"Cookie created:"<<cookie_.toStdString()<<"\n";
#endif
	}
}	

void MascotRemoteQuery::endRun()
{
  if (http_->state() != QHttp::Unconnected) http_->close();
  emit done();
}


void MascotRemoteQuery::httpDone(bool error)
{
	if (error)
	{
    error_message_ = String("Mascot Server replied: '") + String(http_->errorString().toStdString()) + "'";
    endRun();
	}

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
  

	QByteArray new_bytes = http_->readAll();
#ifdef MASCOTREMOTEQUERY_DEBUG
	cerr << "Response of query: " << "\n";
	QTextDocument doc;
	doc.setHtml(new_bytes.constData());
	cerr << doc.toPlainText().toStdString() << "\n";
#endif
	
  if (QString(new_bytes).trimmed().size() == 0 )
  {
    error_message_ = "Error: Reply from mascot server is empty! Possible server overload - see the Mascot Admin!";
    endRun();
  }

	//Successful login? fire off the search
	if (new_bytes.contains("Logged in successfuly")) // this is NOT a typo. Mascot writes 'successfuly' that way!
	{
		emit loginDone();
	} 
	else if (new_bytes.contains("Error: You have entered an invalid password"))
	{
		error_message_ = "Error: You have entered an invalid password";
    endRun();
	}
	else if (new_bytes.contains("is not a valid user"))
	{
		error_message_ = "Error: Username is not valid";
    endRun();
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

		results_path_ = "";
    results_path_.append("/mascot/cgi/export_dat_2.pl?file=");        
    results_path_.append(rx.cap(1));
    
#ifdef MASCOTREMOTEQUERY_DEBUG
		cerr << "Results path to export: " << results_path_.toStdString() << "\n";
#endif
		results_path_.append("&show_same_sets=1&show_unassigned=1&show_queries=1&do_export=1&export_format=XML&pep_rank=1&_sigthreshold=0.99&_showsubsets=1&show_header=1&prot_score=1&pep_exp_z=1&pep_score=1&pep_seq=1&pep_homol=1&pep_ident=1&show_mods=1&pep_var_mod=1&protein_master=1&prot_score=1&search_master=1&show_header=1&show_params=1&pep_scan_title=1&query_qualifiers=1&query_peaks=1&query_raw=1&query_title=1&pep_expect=1&peptide_master=1");
		
		if (param_.getValue("query_master").toBool())
		{
			results_path_.append("&query_master=1");
		}
		else
		{
			results_path_.append("&query_master=0");
		}
			//Finished search, fire off results retrieval
		emit queryDone();
	} 
	else 
	{	
		// check whether Mascot responded using an error code e.g. [M00440], pipe through results else
		QString response_text = new_bytes;
		QRegExp mascot_error_regex("\\[M[0-9][0-9][0-9][0-9][0-9]\\]");
		if (response_text.contains(mascot_error_regex))
		{
			QTextDocument doc;
			doc.setHtml(response_text);
			error_message_ = doc.toPlainText().toStdString();
      endRun();
		}
		else
		{
			// seems to be fine so grab the xml
			#ifdef MASCOTREMOTEQUERY_DEBUG
			cerr << "Get the XML File" << "\n";
			#endif
			mascot_xml_ = new_bytes;
      endRun();
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
		// clear all content from this class
		delete http_;
		http_ = new QHttp(this);
		http_->setHost(((String)param_.getValue("hostname")).c_str(), (UInt)param_.getValue("host_port"));

		cookie_ = "";
		mascot_xml_ = "";
		results_path_ = "";

    to_ = param_.getValue("timeout");
    timeout_.setInterval(1000 * to_);

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
}
