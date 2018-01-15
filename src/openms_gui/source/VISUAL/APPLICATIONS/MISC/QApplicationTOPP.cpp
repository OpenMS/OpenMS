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

#include <cstdio>
#include <cstdlib>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>

//Qt
#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>
#include <QMessageBox>
#include <QFile>
#include <QFileOpenEvent>


namespace OpenMS
{

  QApplicationTOPP::QApplicationTOPP(int& argc, char** argv) :
    QApplication(argc, argv)
  {
    // register GUI ProgressLogger that can be used in GUI tools
    Factory<ProgressLogger::ProgressLoggerImpl>::registerProduct(GUIProgressLoggerImpl::getProductName(), &GUIProgressLoggerImpl::create);

    // set plastique style unless windows / mac style is available
    if (QStyleFactory::keys().contains("windowsxp", Qt::CaseInsensitive))
    {
      this->setStyle("windowsxp");
    }
    else if (QStyleFactory::keys().contains("macintosh", Qt::CaseInsensitive))
    {
      this->setStyle("macintosh");
    }
    else if (QStyleFactory::keys().contains("plastique", Qt::CaseInsensitive))
    {
      this->setStyle("plastique");
    }

    // customize look and feel via Qt style sheets
    String filename = File::find("GUISTYLE/qtStyleSheet.qss");
    QFile fh(filename.toQString());
    fh.open(QFile::ReadOnly);
    QString style_string = QLatin1String(fh.readAll());
    //std::cerr << "Stylesheet content: " << style_string.toStdString() << "\n\n\n";
    this->setStyleSheet(style_string);
  }

  QApplicationTOPP::~QApplicationTOPP()
  {
  }

  /*
    @brief: Catch exceptions in Qt GUI applications, preventing ungraceful exit

    Re-implementing QApplication::notify() to catch exception thrown in event handlers (which is most likely OpenMS code).
  */
  bool QApplicationTOPP::notify(QObject* rec, QEvent* ev)
  {
    // this is called quite often (whenever a signal is fired), so mind performance!
    try
    {
      return QApplication::notify(rec, ev);
    }
    catch (Exception::BaseException& e)
    {
      String msg = String("Caught exception: '") + e.getName() + "' with message '" + e.getMessage() + "'";
      LOG_ERROR << msg << "\n";
      QMessageBox::warning(nullptr, QString("Unexpected error occurred"), msg.toQString());
      return false;
      // we could also exit() here... but no for now
    }

    return false; // never reached, so return value does not matter
  }

  bool QApplicationTOPP::event(QEvent* event)
  {
    switch (event->type())
    {
    case QEvent::FileOpen:
      emit fileOpen(static_cast<QFileOpenEvent*>(event)->file());
      return true;

    default:
      return QApplication::event(event);
    }
  }

}
