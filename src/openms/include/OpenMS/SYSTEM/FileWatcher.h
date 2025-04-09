// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <map>

//Qt
#include <QtCore/QFileSystemWatcher>

//STL
#include <map>

namespace OpenMS
{
  class String;

  /**
      @brief Watcher that monitors file changes.

      This class can be used similar to QFileSystemWatcher.
      Additionally it offers a delayed fileChanged signal.

      This behaviour is required for the following reason:
      Normally QFileSystemWatcher emits a signal every time a file is changed.
      This causes several signals for large files (one for each flush of the buffer).

      @ingroup System
  */
  class OPENMS_DLLAPI FileWatcher :
    public QFileSystemWatcher       //find out why ICC requires public instead of protected
  {
    Q_OBJECT

public:
    /// Constructor
    FileWatcher(QObject * parent = nullptr);

    /// Destructor
    ~FileWatcher() override;

    ///Sets the delay in seconds (default: 1s)
    inline void setDelayInSeconds(double delay)
    {
      delay_in_seconds_ = delay;
    }

    ///Adds a file to the watcher
    inline void addFile(const String & path)
    {
      QFileSystemWatcher::addPath(path.toQString());
    }

    ///removes a file from the watcher
    inline void removeFile(const String & path)
    {
      QFileSystemWatcher::removePath(path.toQString());
    }

signals:
    ///Delayed file change signal
    void fileChanged(const String &);

protected slots:
    /// Slot that is connected to the fileChanged signal in order to track the changes
    void monitorFileChanged_(const QString & name);
    /// Slot that is called when the delay is over
    void timerTriggered_();

protected:
    /// A map that links timer name and file
    std::map<QString, QString> timers_;
    /// Delay (seconds)
    double delay_in_seconds_;
  };

  // OPENMS_DLLAPI extern FileWatcher myFileWatcher_instance;
}

