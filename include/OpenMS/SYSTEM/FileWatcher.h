// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
	class OPENMS_DLLAPI FileWatcher
		: public QFileSystemWatcher //find out why ICC requires public instead of protected
	{
    Q_OBJECT
    
    public:
    	///Constructor
    	FileWatcher(QObject* parent=0);
			
			///Destructor
			virtual ~FileWatcher();
			    	
    	///Sets the delay in seconds (default: 1s)
    	inline void setDelayInSeconds(DoubleReal delay)
    	{
    		delay_in_seconds_ = delay;
    	}
    	
    	///Adds a file to the watcher
    	inline void addFile(const String& path)
    	{
    		QFileSystemWatcher::addPath(path.toQString());
    	}

    	///removes a file from the watcher
    	inline void removeFile(const String& path)
    	{
    		QFileSystemWatcher::removePath(path.toQString());
    	}
    	
    signals:
    	///Delayed file change signal
    	void fileChanged(const String&);
    
    protected slots:
			/// Slot that is connected to the fileChanged signal in order to track the changes
    	void monitorFileChanged_(const QString& name);
			/// Slot that is called when the delay is over
    	void timerTriggered_();

    protected:
    	/// A map that links timer name and file
    	Map<QString,QString> timers_;
			/// Delay (seconds)
    	DoubleReal delay_in_seconds_;
	};

}

#endif // OPENMS_SYSTEM_FILE_H
