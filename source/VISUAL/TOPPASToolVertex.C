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
#include <OpenMS/VISUAL/TOPPASScene.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>

namespace OpenMS
{
	UInt TOPPASToolVertex::id_counter = 1;

	TOPPASToolVertex::TOPPASToolVertex()
		:	TOPPASVertex(),
			name_(),
			type_(),
			param_(),
			id_(id_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type)
		: TOPPASVertex(),
			name_(name),
			type_(type),
			param_(),
			id_(id_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs)
		:	TOPPASVertex(rhs),
			name_(rhs.name_),
			type_(rhs.type_),
			param_(rhs.param_),
			id_(id_counter++)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
	}

	TOPPASToolVertex::~TOPPASToolVertex()
	{
	
	}
	
	TOPPASToolVertex& TOPPASToolVertex::operator= (const TOPPASToolVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		
		param_ = rhs.param_;
		name_ = rhs.name_;
		type_ = rhs.type_;
		
		return *this;
	}
	
	void TOPPASToolVertex::initParam_()
	{
		Param tmp_param;
		String ini_file = "TOPPAS_" + name_ + "_";
		if (type_ != "")
		{
			ini_file += type_ + "_";
		}
		ini_file += QString::number(id_) + ".tmp.ini";
		
		String call = name_ + " -write_ini " + ini_file + " -log " + ini_file + ".log";
		if (type_ != "")
		{
			call += " -type " + type_;
		}
		
		if (system(call.c_str()) != 0)
		{
			QMessageBox::critical(0,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
			return;
		}
		else if(!File::exists(ini_file))
		{
			QMessageBox::critical(0,"Error",(String("Could not open '")+ini_file+"'!").c_str());
			return;
		}
		
		tmp_param.load((ini_file).c_str());
		param_=tmp_param.copy(name_+":1:",true);
		param_.remove("log");
		param_.remove("no_progress");
		param_.remove("debug");
	}
	
	void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
		String default_dir = "";
		
		TOPPASToolConfigDialog dialog(parent_widget, param_, default_dir, name_, type_);
		dialog.exec();
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
	}
	
	void TOPPASToolVertex::getInputFiles(QVector<IOInfo>& input_infos)
	{
		getFiles_(input_infos,true);
	}
	
	void TOPPASToolVertex::getOutputFiles(QVector<IOInfo>& output_infos)
	{
		getFiles_(output_infos,false);
	}
	
	void TOPPASToolVertex::getFiles_(QVector<IOInfo>& io_infos, bool input_files)
	{
		String search_tag = input_files ? "input file" : "output file";
		
		io_infos.clear();
		
		for (Param::ParamIterator it = param_.begin(); it != param_.end(); ++it)
		{
			if (it->tags.count(search_tag))
			{
				StringList valid_types;
				
				const String& desc = it->description;
				String::SizeType index = desc.find("valid formats",0);
				if (index != String::npos)
				{
					String::SizeType types_start_pos = desc.find("'",index) + 1;
					String::SizeType types_length = desc.find("'",types_start_pos) - types_start_pos;
					String types_string = desc.substr(types_start_pos, types_length);
					if (types_string.find(",",0) == String::npos)
					{
						valid_types.push_back(types_string.trim());
					}
					else
					{
						types_string.split(',', valid_types);
					}
				}
				
				IOInfo io_info;
				io_info.param_name = it->name;
				io_info.valid_types = valid_types;
				if (it->value.valueType() == DataValue::STRING_LIST)
				{
					io_info.type = IOInfo::IOT_LIST;
				}
				else if (it->value.valueType() == DataValue::STRING_VALUE)
				{
					io_info.type = IOInfo::IOT_FILE;
				}
				else
				{
					std::cerr << "this should not happen\n" << std::endl;
				}
				io_infos.push_back(io_info);
			}
		}
	}

	void TOPPASToolVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
		
		QPainterPath path;
		path.addRect(-70.0, -60.0, 140.0, 120.0);		
 		painter->drawPath(path);
 		
		if (type_ == "")
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), name_.toQString());
		}
		else
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), -(int)(text_boundings.height()/3.0), name_.toQString());
			text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, type_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), +(int)(text_boundings.height()/1.33), type_.toQString());
		}
	}
	
	QRectF TOPPASToolVertex::boundingRect() const
	{
		return QRectF(-71,-61,142,122);
	}
	
	QPainterPath TOPPASToolVertex::shape () const
	{
		QPainterPath shape;
		shape.addRect(-71.0, -61.0, 142.0, 122.0);				
		return shape;
	}
	
	const String& TOPPASToolVertex::getName()
	{
		return name_;
	}
	
	const String& TOPPASToolVertex::getType()
	{
		return type_;
	}
}

