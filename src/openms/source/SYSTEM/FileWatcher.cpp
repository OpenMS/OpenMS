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
    //Connect the slot for monitoring file changes
    connect(this, SIGNAL(fileChanged(const QString &)), this, SLOT(monitorFileChanged_(const QString &)));
  }

  FileWatcher::~FileWatcher()
  {
  }

  void FileWatcher::monitorFileChanged_(const QString & name)
  {
    //cout << "File changed: " << String(name) << endl;
    //Look up if there is already a timer for this file
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
