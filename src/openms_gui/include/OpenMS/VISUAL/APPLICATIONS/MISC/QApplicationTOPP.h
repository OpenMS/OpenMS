// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//Qt
#include <QtWidgets/QApplication>

namespace OpenMS
{
  /**
    @brief Extension to the QApplication for running TOPPs GUI tools.

    Basically re-implements notify of QApplication to prevent ungraceful exit.
  */
  class OPENMS_GUI_DLLAPI QApplicationTOPP :
    public QApplication
  {

    Q_OBJECT

public:
    /// Constructor (no NOT remove the "&" from argc, since Qt will segfault on some platforms otherwise!)
    QApplicationTOPP(int& argc, char** argv);

    /// Destructor
    ~QApplicationTOPP() override;

    /**
      @brief: Catch exceptions in Qt GUI applications, preventing ungraceful exit

      Re-implementing QApplication::notify() to catch exception thrown in event
      handlers (which is most likely OpenMS code).
    */
    bool notify(QObject* rec, QEvent* ev) override;

    /**
      Reimplemented from QApplication, to handle QEvent::FileOpen to enable handling of odoc event on MacOSX
    */
    bool event(QEvent*) override;

    /**
      @brief Show the About-Dialog with License and Citation for all GUI tools

      @param parent Parent widget (usually 'this')
      @param toolname name of the tool (used as heading)
    */
    static void showAboutDialog(QWidget* parent, const QString& toolname);


signals:
    void fileOpen(QString file);

  };

}

