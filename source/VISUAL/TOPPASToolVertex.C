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
#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>

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
			id_(id_counter++),
			finished_(false)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type, const String& tmp_path)
		: TOPPASVertex(),
			name_(name),
			type_(type),
			tmp_path_(tmp_path),
			param_(),
			id_(id_counter++),
			finished_(false)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs)
		:	TOPPASVertex(rhs),
			name_(rhs.name_),
			type_(rhs.type_),
			tmp_path_(rhs.tmp_path_),
			param_(rhs.param_),
			id_(rhs.id_),
			finished_(rhs.finished_)
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
		tmp_path_ = rhs.tmp_path_;
		id_ = rhs.id_;
		finished_ = rhs.finished_;
		
		return *this;
	}
	
	void TOPPASToolVertex::initParam_()
	{
		Param tmp_param;
		String ini_file = tmp_path_ + "TOPPAS_" + name_ + "_";
		if (type_ != "")
		{
			ini_file += type_ + "_";
		}
		ini_file += File::getUniqueName() + "_tmp.ini";
		
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
		
		// handled by TOPPAS anyway:
		param_.remove("type");
	}
	
	void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
		String default_dir = "";
		
		// remove entries that are handled by edges already, user should not see them
		QVector<Param::ParamEntry> hidden_entries;
		QVector<IOInfo> input_infos;
		getInputParameters(input_infos);
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int index = (*it)->getTargetInParam();
			if (index < 0)
			{
				continue;
			}
			
			const String& name = input_infos[index].param_name;
			if (param_.exists(name))
			{
				hidden_entries.push_back(param_.getEntry(name));
				param_.remove(name);
			}
		}
		
		QVector<IOInfo> output_infos;
		getOutputParameters(output_infos);
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			int index = (*it)->getSourceOutParam();
			if (index < 0)
			{
				continue;
			}
			
			const String& name = output_infos[index].param_name;
			if (param_.exists(name))
			{
				hidden_entries.push_back(param_.getEntry(name));
				param_.remove(name);
			}
		}
		
		TOPPASToolConfigDialog dialog(parent_widget, param_, default_dir, name_, type_);
		dialog.exec();
		
		// restore the removed entries
		foreach (const Param::ParamEntry& pe, hidden_entries)
		{
			StringList tags;
			for (std::set<String>::const_iterator it = pe.tags.begin(); it != pe.tags.end(); ++it)
			{
				tags.push_back(*it);
			}
			param_.setValue(pe.name, pe.value, pe.description, tags);
		}
		
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
	}
	
	void TOPPASToolVertex::getInputParameters(QVector<IOInfo>& input_infos)
	{
		getParameters_(input_infos,true);
	}
	
	void TOPPASToolVertex::getOutputParameters(QVector<IOInfo>& output_infos)
	{
		getParameters_(output_infos,false);
	}
	
	void TOPPASToolVertex::getParameters_(QVector<IOInfo>& io_infos, bool input_params)
	{
		String search_tag = input_params ? "input file" : "output file";
		
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
	
	void TOPPASToolVertex::compute()
	{
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv && !(tv->isFinished()))
			{
				tv->compute();
			}
		}
		
		// all inputs ready now --> write param to ini file and run tool
		String ini_file_name = tmp_path_+"TOPPAS_" + name_ + "_" + type_ + "_" + File::getUniqueName() + ".ini";
		param_.store(ini_file_name);
		
		String call = name_+" -ini "+ini_file_name+" -log "+ini_file_name+".log";
		if (type_ != "")
		{
			call += " -type "+type_;
		}
		
		// add all input file parameters
		QVector<IOInfo> in_params;
		getInputParameters(in_params);
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int param_index = (*it)->getTargetInParam();
			if (param_index < 0)
			{
				std::cerr << "blub" << std::endl;
				continue;
			}
			call += " -" + in_params[param_index].param_name;
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv)
			{
				int out_param_index = (*it)->getSourceOutParam();
				if (out_param_index < 0)
				{
					std::cerr << "blub" << std::endl;
					continue;
				}
				call += " " + String(tv->output_file_names_[out_param_index].join(" "));
				continue;
			}
			TOPPASInputFileVertex* ifv = qobject_cast<TOPPASInputFileVertex*>((*it)->getSourceVertex());
			if (ifv)
			{
				call += " " + String(ifv->getFilename());
				continue;
			}
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
			if (iflv)
			{
				call += " " + String(iflv->getFilenames().join(" "));
				continue;
			}
		}
		
		// add all output file parameters
		QVector<IOInfo> out_params;
		getOutputParameters(out_params);
		
		output_file_names_.clear();
		output_file_names_.resize(out_params.size()); // TODO: better solution?
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			int param_index = (*it)->getSourceOutParam();
			if (param_index < 0)
			{
				std::cerr << "blub" << std::endl;
				continue;
			}
			call += " -" + out_params[param_index].param_name;
			
			if (out_params[param_index].param_name == "out" && out_params[param_index].type == IOInfo::IOT_LIST)
			{
				// TODO think about how to determine number of output arguments for output params with list type!
				// for now, use this workaround for "out" (as many as for "in")... what about the others?!
				for (int i = 0; i < in_params.size(); ++i)
				{
					output_file_names_[param_index].push_back((tmp_path_ + name_ + "_" + type_ + "_" + File::getUniqueName() + ".out").toQString());
				}
				call += " " + String(output_file_names_[param_index].join(" "));
			}
			else if (out_params[param_index].type == IOInfo::IOT_FILE)
			{
				output_file_names_[param_index].push_back((tmp_path_ + name_ + "_" + type_ + "_" + File::getUniqueName() + ".out").toQString());
				call += " " + String(output_file_names_[param_index].join(" "));
			}
			else // IOT_LIST
			{
				std::cerr << "I don't know what to do in this situation, yet." << std::endl;
			}
		}
		
		if(system(call.c_str())!=0)
		{
			QMessageBox::critical(0,"Error",("Something went wrong!\n\nCheck the log file " + ini_file_name + ".log for possible reasons. If it does not exist, make sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
		}
		
		finished_ = true;
	}
	
	bool TOPPASToolVertex::isFinished()
	{
		return finished_;
	}
	
	void TOPPASToolVertex::setFinished(bool b)
	{
		finished_ = b;
	}
	
}

