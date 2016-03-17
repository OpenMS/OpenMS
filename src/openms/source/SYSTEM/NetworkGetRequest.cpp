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

#include <OpenMS/SYSTEM/NetworkGetRequest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <QtGui/QTextDocument>
#include <iostream>

// #define QUERY_DEBUG
using namespace std;

namespace OpenMS
{

  NetworkGetRequest::NetworkGetRequest(QObject* parent) :
    QObject(parent)
  {
    http_ = new QHttp();
    to_ = 5;
    timeout_.setInterval(1000 * to_);
  }

  NetworkGetRequest::~NetworkGetRequest()
  {
    if (http_->state() != QHttp::Unconnected)
    {
#ifdef QUERY_DEBUG
      std::cerr << "Aborting open connection!\n";
#endif
      http_->abort(); // hardcore close connection (otherwise server might have too many dangling requests)
    }
    delete http_;
  }

  void NetworkGetRequest::timedOut()
  {
    http_->abort();
  }

  void NetworkGetRequest::setHost(QString host)
  {
    host_ = host;
  }

  void NetworkGetRequest::setPath(QString path)
  {
    path_ = path;
  }

  void NetworkGetRequest::run()
  {
    // Due to the asynchronous nature of QHttp::request (and the resulting use
    // of signals and slots), the information flow in this class is not very
    // clear. After the initial call to "run", the steps are:
    // 1. send query
    // 2. read query result, prepare exporting (function "httpDone")

    http_->setHost(host_.c_str(), 80);

    connect(http_, SIGNAL(requestFinished(int, bool)), this, SLOT(httpRequestFinished(int, bool)));
    //connect(http_, SIGNAL(requestStarted(int)), this, SLOT(httpRequestStarted(int)));
    connect(http_, SIGNAL(done(bool)), this, SLOT(httpDone(bool)));
    connect(http_, SIGNAL(stateChanged(int)), this, SLOT(httpStateChanged(int)));
    connect(http_, SIGNAL(readyRead(const QHttpResponseHeader &)), this, SLOT(readyReadSlot(const QHttpResponseHeader &)));
    connect(http_, SIGNAL(responseHeaderReceived(const QHttpResponseHeader &)), this, SLOT(readResponseHeader(const QHttpResponseHeader &)));
    connect(&timeout_, SIGNAL(timeout()), this, SLOT(timedOut()));

    query();
  }

  void NetworkGetRequest::query()
  {
#ifdef QUERY_DEBUG
    cerr << "void NetworkGetRequest::getResults()" << "\n";
#endif
    QHttpRequestHeader header;
    header.setRequest("GET", path_);
    header.setRequest("Host", host_);
    header.setValue("Accept", "text/plain");
    header.setValue("Keep-Alive", "300");
    header.setValue("Connection", "keep-alive");

#ifdef QUERY_DEBUG
    logHeader_(header, "request results");
#endif

    http_->request(header);
  }

  void NetworkGetRequest::httpRequestFinished(int requestId, bool error)
  {
    if (error)
    {
      cerr << "NetworkGetRequest: An error occurred (requestId=" << requestId << "): " << http_->errorString().toStdString() << " (QT Error Code: " << int(http_->error()) << ")\n";
    }
#ifdef QUERY_DEBUG
    cerr << "Request Finished Id: " << requestId << "\n";
    cerr << "Error: " << error << "(" << http_->errorString().toStdString() << ")" << "\n";
#endif

  }

  void NetworkGetRequest::httpStateChanged(int
#ifdef QUERY_DEBUG
    state
#endif
    )
  {
#ifdef QUERY_DEBUG
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

  void NetworkGetRequest::readyReadSlot(const QHttpResponseHeader& /* resp */)
  {
    if (to_ > 0)
      timeout_.start(); // reset timeout
  }

  void NetworkGetRequest::readResponseHeader(const QHttpResponseHeader& response_header)
  {
    if (response_header.statusCode() >= 400)
    {
      error_message_ = String("NetworkGetRequest: The server returned an error status code '") + response_header.statusCode() + "': " + response_header.reasonPhrase() + "\n";
      endRun_();
    }
  }

  QString NetworkGetRequest::getResponse() const
  {
    return QString(response_bytes_);
  }

  void NetworkGetRequest::endRun_()
  {
    if (http_->state() != QHttp::Unconnected)
    {
      http_->clearPendingRequests();
      http_->close();
    }
    emit done();
  }

  void NetworkGetRequest::httpDone(bool error)
  {
#ifdef QUERY_DEBUG
    cerr << "void NetworkGetRequest::httpDone(bool error): ";
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
      error_message_ = String("Server replied: '") + String(http_->errorString().toStdString()) + "'";
      endRun_();
      return;
    }

    QByteArray new_bytes = http_->readAll();
#ifdef QUERY_DEBUG
    cerr << "Response of query: " << "\n";
    QTextDocument doc;
    doc.setHtml(new_bytes.constData());
    cerr << doc.toPlainText().toStdString() << "\n";
#endif

    response_bytes_ = new_bytes;
    endRun_();
  }

  bool NetworkGetRequest::hasError() const
  {
    return error_message_ != "";
  }

  const String& NetworkGetRequest::getErrorMessage() const
  {
    return error_message_;
  }

}
