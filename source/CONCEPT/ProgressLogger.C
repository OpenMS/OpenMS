// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QString>
#include <QtGui/QProgressDialog>

#include <iostream>

using namespace std;

namespace OpenMS
{
	int ProgressLogger::recursion_depth_ = 0;

	ProgressLogger::ProgressLogger()
		:	type_(NONE),
			begin_(0),
			end_(0),
			value_(0),
			dlg_(0),
			stop_watch_(),
			last_invoke_()
	{
	}

	ProgressLogger::~ProgressLogger()
	{
		delete(dlg_);
	}

	void ProgressLogger::setLogType(LogType type) const
	{
		type_ = type;
	}

	ProgressLogger::LogType ProgressLogger::getLogType() const
	{
		return type_;
	}


	void ProgressLogger::startProgress(SignedSize begin, SignedSize end, const String& label) const
	{
		OPENMS_PRECONDITION(begin <= end, "ProgressLogger::init : invalid range!");
		last_invoke_ = time (NULL);
		
		switch (type_)
		{
			case CMD:
				begin_ = begin;
				end_ = end;
				if ( recursion_depth_ ) cout << '\n';
				cout << string(2*recursion_depth_,' ') << "Progress of '" << label << "':" << endl;
				stop_watch_.reset();
				stop_watch_.start();
				break;
			case GUI:
				begin_ = begin;
				end_ = end;
				if(!dlg_) dlg_ = new QProgressDialog(label.c_str(), QString(), int(begin), int(end));
				dlg_->setWindowTitle(label.c_str());
				dlg_->setWindowModality(Qt::WindowModal);
				dlg_->show();
				break;
			case NONE:
				break;
		};
		++recursion_depth_;
		return;
	}

	void ProgressLogger::setProgress(SignedSize value) const
	{
		// update only if at least 1 second has passed
		if (last_invoke_ == time (NULL)) return;
		
		last_invoke_ = time (NULL);
	
		switch (type_)
		{
			case CMD:
				if (begin_==end_)
				{
					cout << '.' << flush;
				}
				else if (value<begin_ || value >end_)
				{
					cout << "ProgressLogger: Invalid progress value '" << value
							 << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << endl;
				}
				else
				{
					cout << '\r' << string(2*recursion_depth_,' ') << QString::number(Real(value -begin_) / Real(end_ - begin_) * 100.0,'f',2).toStdString()  << " %               ";
					cout << flush;
				}
				break;
			case GUI:
				if (value<begin_ || value >end_)
				{
					cout << "ProgressLogger: Invalid progress value '" << value << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << endl;
				}
				else
				{
					if (dlg_)
					{
						dlg_->setValue((int)value);
					}
					else
					{
						cout << "ProgressLogger warning: 'setValue' called before 'startProgress'!" << endl;
					}
				}
				break;
			case NONE:
				break;
		};
	}

	void ProgressLogger::endProgress() const
	{
		if (recursion_depth_) --recursion_depth_;
		switch (type_)
		{
			case CMD:
				stop_watch_.stop();
				if (begin_==end_)
				{
					if ( recursion_depth_ ) cout << '\n';
          cout << endl << string(2*recursion_depth_,' ') << "-- done [took " << String::number(stop_watch_.getCPUTime(),3) << " s(CPU), " << String::number(stop_watch_.getClockTime(),3) << " s(Wall)] -- " << endl;
				}
				else
				{
					cout << '\r' << string(2*recursion_depth_,' ') << "-- done [took " << String::number(stop_watch_.getCPUTime(),3) << " s(CPU), " << String::number(stop_watch_.getClockTime(),3) << " s(Wall)] -- " << endl;
				}
				break;
			case GUI:
				if (dlg_)
				{
					dlg_->setValue((int)end_);
				}
				else
				{
					cout << "ProgressLogger warning: 'endProgress' called before 'startProgress'!" << endl;
				}
				break;
			case NONE:
				break;
		};
	}

}//namespace OpenMS
