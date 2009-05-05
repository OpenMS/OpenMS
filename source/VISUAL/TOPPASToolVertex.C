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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>

namespace OpenMS
{
	UInt TOPPASToolVertex::instance_counter = 1;

	TOPPASToolVertex::TOPPASToolVertex()
		:	TOPPASVertex(),
			param_(),
			instance_nr_(instance_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(250,200,0);
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type)
		: TOPPASVertex(name, type),
			param_(),
			instance_nr_(instance_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(250,200,0);
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs)
		:	TOPPASVertex(rhs),
			param_(rhs.param_),
			instance_nr_(instance_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(250,200,0);
	}

	TOPPASToolVertex::~TOPPASToolVertex()
	{
	
	}
	
	TOPPASToolVertex& TOPPASToolVertex::operator= (const TOPPASToolVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		
		param_ = rhs.param_;
		
		return *this;
	}
	
	void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
		String default_dir = "";
		String ini_file = name_ + "__";
		if (type_ != "")
		{
			ini_file += type_ + "__";
		}
		ini_file += QString::number(instance_nr_) + ".ini";
		
		String call = name_ + " -write_ini " + ini_file + " -log " + ini_file + ".log";
		if (type_ != "")
		{
			call += " -type " + type_;
		}
		
		if (system(call.c_str()) != 0)
		{
			QMessageBox::critical(parent_widget,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
		}
		else if(!File::exists(ini_file))
		{
			QMessageBox::critical(parent_widget,"Error",(String("Could not open '")+ini_file+"'!").c_str());
		}
		else
		{
			TOPPASToolConfigDialog dialog(parent_widget, ini_file, default_dir, name_, type_, instance_nr_);
			if (dialog.exec())
			{
				// ...
			}
		}
	}
}

