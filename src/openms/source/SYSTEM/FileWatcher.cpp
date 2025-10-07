// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/FileWatcher.h>
#include <QtCore/QTimer>

using namespace std;

namespace OpenMS
{
  FileWatcher::FileWatcher(QObject * parent) :
    QFileSystemWatcher(parent),
    timers_(),
    delay_in_seconds_(1.0)
  {
    // Connect the slot for monitoring file changes
    connect(this, &FileWatcher::fileChanged, [this](const String& s) { monitorFileChanged_(s.toQString()); });
  }

  FileWatcher::~FileWatcher() = default;

  void FileWatcher::monitorFileChanged_(const QString & name)
  {
    //cout << "File changed: " << String(name) << endl;
    // Look up if there is already a timer for this file
    QTimer * timer = nullptr;
    for (map<QString, QString>::const_iterator it = timers_.begin(); it != timers_.end(); ++it)
    {
      if (it->second == name)     //we found the timer name and id
      {
        //cout << " - Found timer name: " << String(it->second) << endl;
        //search for the timer instance with the corresponding Id
        timer = findChild<QTimer *>(it->first);
      }
    }

    //timer does not exist => create and start a new one
    if (!timer)
    {
      //static timer counter
      static int timer_id = 0;
      //cout << " - no timer found => creating a new one with name: ";
      timer = new QTimer(this);
      timer->setInterval((int)(1000.0 * delay_in_seconds_));
      timer->setSingleShot(true);
      timer->setObjectName(QString::number(++timer_id));
      connect(timer, SIGNAL(timeout()), this, SLOT(timerTriggered_()));
      timer->start();
      timers_[QString::number(timer_id)] = name;
      //cout << timer_id << endl;

    }
    //timer exists => restart it as the file changed another time
    else
    {
      //cout << " - timer found => resetting" << endl;
      timer->start();
    }
  }

  void FileWatcher::timerTriggered_()
  {
    //cout << "Timer activated" << endl;
    //get the timer instance
    QTimer * timer = qobject_cast<QTimer *>(sender());
    //emit the final for the file corresponding to the timer name
    //cout << " - timer name: " << String(timer->objectName()) << endl;
    //cout << " - timer file: " << String(timers_[timer->objectName()]) << endl;
    emit fileChanged(String(timers_[timer->objectName()]));
    //erase the timer name from the list
    timers_.erase(timer->objectName());
  }

  //OPENMS_DLLAPI FileWatcher myFileWatcher_instance; // required, such that the moc file get generated during building OpenMS.dll, not later during OpenMS_GUI.dll as DLL flags are wrong then

} // namespace OpenMS
