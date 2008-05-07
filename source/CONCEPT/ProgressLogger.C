// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <QtCore/QString>
#include <QtGui/QProgressDialog>

#include <iostream>

using namespace std;

namespace OpenMS
{
	ProgressLogger::ProgressLogger()
		:	type_(NONE),
			begin_(0),
			end_(0),
			value_(0),
			dlg_(0)
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
	
	void ProgressLogger::startProgress(UInt begin, UInt end, const String& label) const
	{
		OPENMS_PRECONDITION(begin <= end, "ProgressLogger::init : invalid range!");
		
		switch (type_)
		{
			case CMD:
				begin_ = begin;
				end_ = end;
				cout << "Progress of '" << label << "':" << endl;		
				break;
			case GUI:
				begin_ = begin;
				end_ = end;
				if(!dlg_) dlg_ = new QProgressDialog();
				dlg_->setRange(begin,end);
				dlg_->setLabelText(label.c_str());
				dlg_->setWindowTitle(label.c_str());
				dlg_->show();
				break;
			case NONE:
				break;
		};	
	}
	
	void ProgressLogger::setProgress(UInt value) const
	{		
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
					cout << '\r' <<QString::number(Real(value -begin_) / Real(end_ - begin_) * 100.0,'f',2).toStdString()  << " %               ";
				}
				break;
			case GUI:
				if (begin_==end_)
				{
					dlg_->setValue(value);
				}
				else if (value<begin_ || value >end_)
				{
					cout << "ProgressLogger: Invalid progress value '" << value 
							 << "'. Should be between '" << begin_ << "' and '" << end_ << "'!" << endl;
				}
				else
				{
					dlg_->setValue(value);
				}	
				break;
			case NONE:
				break;
		};	
	}
	
	void ProgressLogger::endProgress() const
	{
		switch (type_)
		{
			case CMD:
				if (begin_==end_)
				{
					cout << endl << " -- done --  " << endl;
				}
				else
				{
					cout << "\r -- done --          " << endl;
				}
				break;
			case GUI:
				dlg_->setValue(end_);
				break;
			case NONE:
				break;
		};	
	}

}//namespace OpenMS
