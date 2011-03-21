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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/FileWatcher.h>
#include <QtCore/QTimer>

using namespace std;

namespace OpenMS 
{
	FileWatcher::FileWatcher(QObject* parent)
		: QFileSystemWatcher(parent),
			timers_(),
			delay_in_seconds_(1.0)
	{
		//Connect the slot for monitoring file changes
		connect(this,SIGNAL(fileChanged(const QString&)),this,SLOT(monitorFileChanged_(const QString&)));
	}

	FileWatcher::~FileWatcher()
	{
	}
	
	void FileWatcher::monitorFileChanged_(const QString& name)
	{
		//static timer counter
		static int timer_id = 0;

    //cout << "File changed: " << String(name) << endl;
		//Look up if there is already a timer for this file
		QTimer* timer = 0;	  
		for (map<QString,QString>::const_iterator it=timers_.begin(); it!=timers_.end(); ++it)
		{
			if (it->second==name) //we found the timer name and id
			{
				//cout << " - Found timer name: " << String(it->second) << endl;
				//search for the timer instance with the corresponding Id
				timer = findChild<QTimer*>(it->first);				
			}
		}
		
		//timer does not exist => create and start a new one
		if (!timer) 
		{
			//cout << " - no timer found => creating a new one with name: ";
			timer = new QTimer(this);
			timer->setInterval((int)(1000.0*delay_in_seconds_));
			timer->setSingleShot(true);
			timer->setObjectName(QString::number(++timer_id));
			connect(timer,SIGNAL(timeout()),this,SLOT(timerTriggered_()));
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
  	QTimer* timer = qobject_cast<QTimer*>(sender());
  	//emit the final for the file corresponding to the timer name
		//cout << " - timer name: " << String(timer->objectName()) << endl;
		//cout << " - timer file: " << String(timers_[timer->objectName()]) << endl;
  	emit fileChanged(String(timers_[timer->objectName()]));
  	//erase the timer name from the list
  	timers_.erase(timer->objectName());
  }

  //OPENMS_DLLAPI FileWatcher myFileWatcher_instance; // required, such that the moc file get generated during building OpenMS.dll, not later during OpenMS_GUI.dll as DLL flags are wrong then

} // namespace OpenMS
