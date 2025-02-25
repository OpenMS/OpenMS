// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/VISUAL/GUIProgressLoggerImpl.h>

//Qt
#include <QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QStyleFactory>
#include <QtWidgets/QPushButton>
#include <QMessageBox>
#include <QFile>
#include <QFileOpenEvent>
#include <QLibraryInfo>


namespace OpenMS
{

  QApplicationTOPP::QApplicationTOPP(int& argc, char** argv) :
    QApplication(argc, argv)
  {
    // inject the GUIProgressLoggerImpl to be used by OpenMS lib via an extern variable
    make_gui_progress_logger = 
      []() -> ProgressLogger::ProgressLoggerImpl* { return new GUIProgressLoggerImpl(); };

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

  QApplicationTOPP::~QApplicationTOPP() = default;

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
      String msg = String("Caught exception: '") + e.getName() + "' with message '" + e.what() + "'";
      OPENMS_LOG_ERROR << msg << "\n";
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

  void QApplicationTOPP::showAboutDialog(QWidget* parent, const QString& toolname)
  {
    // dialog and grid layout
    QDialog* dlg = new QDialog(parent);
    QGridLayout* grid = new QGridLayout(dlg);
    dlg->setWindowTitle("About " + toolname);

    // image
    QLabel* label = new QLabel(dlg);
    label->setPixmap(QPixmap(":/TOPP_about.png"));
    grid->addWidget(label, 0, 0);

    // text
    QString text = QString("<BR>"
                           "<FONT size=+3>%1</font><BR>"
                           "<BR>"
                           "Version %2 %3"
                           "<BR>"
                           "OpenMS and TOPP is free software available under the<BR>"
                           "BSD 3-Clause License (BSD-new)<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "<BR>"
                           "Any published work based on TOPP and OpenMS shall cite:<BR>%4")
    .arg(toolname)
    .arg(VersionInfo::getVersion().toQString())
    .arg( // if we have a revision, embed it also into the shown version number
      VersionInfo::getRevision().empty() ? "" : QString(" (") + VersionInfo::getRevision().toQString() + ")")
    .arg((TOPPBase::cite_openms.title + "<BR>" + TOPPBase::cite_openms.when_where + "<BR>doi:" + TOPPBase::cite_openms.doi).c_str());

    label = new QLabel(text, dlg);

    grid->addWidget(label, 0, 1, Qt::AlignTop | Qt::AlignLeft);

    // close button
    QPushButton* button = new QPushButton("Close", dlg);
    grid->addWidget(button, 1, 1, Qt::AlignBottom | Qt::AlignRight);
    connect(button, SIGNAL(clicked()), dlg, SLOT(close()));

    // execute
    dlg->exec();
  }

}
