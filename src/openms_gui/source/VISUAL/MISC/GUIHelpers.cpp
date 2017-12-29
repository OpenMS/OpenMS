// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QDesktopServices>
#include <QDir>
#include <QUrl>
#include <QMessageBox>
#include <QString>
#include <QStringList>
#include <QProcess>

namespace OpenMS
{

  void GUIHelpers::openFolder(const QString& folder)
  {
#if defined(__APPLE__)
    QProcess* p = new QProcess();
    p->setProcessChannelMode(QProcess::ForwardedChannels);
    QStringList app_args;
    app_args.append(folder);
    p->start("/usr/bin/open", app_args);
    if (!p->waitForStarted())
    {
      // execution failed
      QMessageBox::warning(0, "Open Folder Error", "The folder '" + folder + "' could not be opened!");
      LOG_ERROR << "Failed to open folder '" << folder.toStdString() << "'" << std::endl;
      LOG_ERROR << p->errorString().toStdString() << std::endl;
    }
#else
    if (!QDir(folder).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + folder, QUrl::TolerantMode))))
    {
      QMessageBox::warning(nullptr, "Open Folder Error", "The folder '" + folder + "' could not be opened!");
    }
#endif
  }

  void GUIHelpers::startTOPPView(const QStringList& args)
  {
    QProcess* p = new QProcess();
    p->setProcessChannelMode(QProcess::ForwardedChannels);
#if defined(__APPLE__)
    // check if we can find the TOPPView.app
    QString app_path = (File::getExecutablePath() + "../../../TOPPView.app").toQString();

    if (File::exists(app_path))
    {
      // we found the app
      QStringList app_args;
      app_args.append("-a");
      app_args.append(app_path);
      app_args.append("--args");
      app_args.append(args);
      p->start("/usr/bin/open", app_args);
    }
    else
    {
      // we could not find the app, try it the Linux way
      QString toppview_executable = (File::findExecutable("TOPPView")).toQString();
      p->start(toppview_executable, args);
    }
#else
    // LINUX+WIN
    QString toppview_executable = (File::findExecutable("TOPPView")).toQString();
    p->start(toppview_executable, args);
#endif


    if (!p->waitForStarted())
    {
      // execution failed
      LOG_ERROR << p->errorString().toStdString() << std::endl;
  #if defined(Q_WS_MAC)
      LOG_ERROR << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
  #endif
    }
  }

  void GUIHelpers::openURL(const QString& target)
  {
    QUrl url_target;

    // add protocol handler if none is given
    if (!(target.startsWith("http://") || target.startsWith("https://")))
    {
      // we expect all unqualified urls to be file urls
      try
      {
        String local_url = File::findDoc(target);
        url_target = QUrl::fromLocalFile(local_url.toQString());
      }
      catch (Exception::FileNotFound&)
      {
        // we fall back to the web url
        url_target = QUrl(QString("http://www.openms.de/current_doxygen/%1").arg(target), QUrl::TolerantMode);
      }
    }
    else
    {
      url_target = QUrl(target, QUrl::TolerantMode);
    }

    if (!QDesktopServices::openUrl(url_target))
    {
      QMessageBox::warning(nullptr,
                           QObject::tr("Error"),
                           QObject::tr("Unable to open\n") + target + 
                           QObject::tr("\n\nPossible reason: security settings or misconfigured Operating System"));
    }
  }

} //namespace OpenMS
