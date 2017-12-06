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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_SYSTEM_FILEWATCHER_H
#define OPENMS_SYSTEM_FILEWATCHER_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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
  class FileWatcher :
    public QFileSystemWatcher       //find out why ICC requires public instead of protected
  {
    Q_OBJECT

public:
    ///Constructor
    OPENMS_DLLAPI FileWatcher(QObject * parent = nullptr);

    ///Destructor
    OPENMS_DLLAPI ~FileWatcher() override;

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
    OPENMS_DLLAPI void fileChanged(const String &);

protected slots:
    /// Slot that is connected to the fileChanged signal in order to track the changes
    OPENMS_DLLAPI void monitorFileChanged_(const QString & name);
    /// Slot that is called when the delay is over
    OPENMS_DLLAPI void timerTriggered_();

protected:
    /// A map that links timer name and file
    Map<QString, QString> timers_;
    /// Delay (seconds)
    double delay_in_seconds_;
  };

  // OPENMS_DLLAPI extern FileWatcher myFileWatcher_instance;
}

#endif // OPENMS_SYSTEM_FILE_H
