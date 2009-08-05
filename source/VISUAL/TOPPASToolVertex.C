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
#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>
#include <QtCore/QFile>

namespace OpenMS
{
	TOPPASToolVertex::TOPPASToolVertex()
		:	TOPPASVertex(),
			name_(),
			type_(),
			param_(),
			finished_(false),
			started_here_(false)
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
			finished_(false),
			started_here_(false)
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
			finished_(rhs.finished_),
			started_here_(rhs.started_here_)
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
		finished_ = rhs.finished_;
		started_here_ = rhs.started_here_;
		
		return *this;
	}
	
	void TOPPASToolVertex::initParam_()
	{
		Param tmp_param;
		ini_file_ = tmp_path_ + "TOPPAS_" + name_ + "_";
		if (type_ != "")
		{
			ini_file_ += type_ + "_";
		}
		ini_file_ += File::getUniqueName() + "_tmp.ini";
		
		String call = name_ + " -write_ini " + ini_file_ + " -log " + ini_file_ + ".log";
		if (type_ != "")
		{
			call += " -type " + type_;
		}
		
		if (system(call.c_str()) != 0)
		{
			QMessageBox::critical(0,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
			return;
		}
		else if(!File::exists(ini_file_))
		{
			QMessageBox::critical(0,"Error",(String("Could not open '")+ini_file_+"'!").c_str());
			return;
		}
		
		tmp_param.load((ini_file_).c_str());
		param_=tmp_param.copy(name_+":1:",true);
		param_.remove("log");
		param_.remove("no_progress");
		param_.remove("debug");
		
		// handled by TOPPAS anyway:
		param_.remove("type");
		
		// remove tmp ini file
		QFile::remove(ini_file_.toQString());
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
		
		//emit somethingHasChanged();
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
					std::cerr << "Unexpected parameter value!" << std::endl;
				}
				io_infos.push_back(io_info);
			}
		}
		// order in param can change --> sort
		qSort(io_infos);
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
	
	void TOPPASToolVertex::runRecursively()
	{
		if (started_here_)
		{
			// make sure pipelines are not run multiple times
			return;
		}
		
		bool we_depend_on_other_tools = false;
		// recursive execution of all parent nodes that are tools
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv)
			{
				we_depend_on_other_tools = true;
				tv->runRecursively();
			}
		}
		if (!we_depend_on_other_tools)
		{
			// start actual pipeline execution here
			started_here_ = true;
			runToolIfInputReady();
		}
	}
	
	void TOPPASToolVertex::runToolIfInputReady()
	{
		//check if everything ready
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv && !tv->isFinished())
			{
				// some tool that we depend on has not finished execution yet --> do not start yet
				return;
			}
		}
		
		// all inputs are ready --> GO!
		param_.store(ini_file_);
		
		QStringList args;
		args	<< "-ini"
					<< ini_file_.toQString()
					<< "-no_progress";
		if (type_ != "")
		{
			args << "-type" << type_.toQString();
		}
		
		// add all input file parameters
		QVector<IOInfo> in_params;
		getInputParameters(in_params);
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int param_index = (*it)->getTargetInParam();
			if (param_index < 0)
			{
				std::cerr << "Input parameter index out of bounds!" << std::endl;
				break;
			}
			args << "-"+(in_params[param_index].param_name).toQString();
			
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv)
			{
				int out_param_index = (*it)->getSourceOutParam();
				if (out_param_index < 0)
				{
					std::cerr << "Output parameter index out of bounds!" << std::endl;
					break;
				}
				args << tv->output_file_names_[out_param_index];
				continue;
			}
			
			TOPPASInputFileVertex* ifv = qobject_cast<TOPPASInputFileVertex*>((*it)->getSourceVertex());
			if (ifv)
			{
				args << ifv->getFilename();
				continue;
			}
			
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
			if (iflv)
			{
				args << iflv->getFilenames();
				continue;
			}
		}
		
		// add all output file parameters
		QVector<IOInfo> out_params;
		getOutputParameters(out_params);
		
		for (int i = 0; i < out_params.size(); ++i)
		{
			// search for an out edge for this parameter
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				int param_index = (*it)->getSourceOutParam();
				
				if (i == param_index)
				{
					args	<< "-"+(out_params[param_index].param_name).toQString()
								<< output_file_names_[param_index]; // can be single file name or list
					break; // (regardless of the number of out edges, every argument must appear only once)
				}
			}
		}
		
		//start process
		QProcess* p = new QProcess();
		p->setProcessChannelMode(QProcess::MergedChannels);
		connect(p,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(executionFinished(int,QProcess::ExitStatus)));
		
		//start process
		p->start(name_.toQString(), args);
		p->waitForStarted();
		emit toolStarted();
	}
	
	void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
	{
		// remove tmp ini file
		QFile::remove(ini_file_.toQString());
		
		if (es != QProcess::NormalExit)
		{
			std::cerr << "TOPP tool crashed!" << std::endl;
			emit toolCrashed();
			return;
		}
		
		if (ec != 0)
		{
			std::cerr << "TOPP tool execution failed! (Exit code: " << ec << ")" << std::endl;
			emit toolFailed();
			return;
		}
		
		finished_ = true;
		// notify all childs that we are finished
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getTargetVertex());
			if (tv)
			{
				tv->runToolIfInputReady();
				continue;
			}
			TOPPASOutputFileVertex* ofv = qobject_cast<TOPPASOutputFileVertex*>((*it)->getTargetVertex());
			if (ofv)
			{
				ofv->finished();
				continue;
			}
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>((*it)->getTargetVertex());
			if (oflv)
			{
				oflv->finished();
				continue;
			}
		}
		
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (p)
		{
			delete p;
		}
		
		emit toolFinished();
	}
	
	bool TOPPASToolVertex::isFinished()
	{
		return finished_;
	}
	
	void TOPPASToolVertex::setFinished(bool b)
	{
		finished_ = b;
	}
	
	const Param& TOPPASToolVertex::getParam()
	{
		return param_;
	}
	
	void TOPPASToolVertex::setParam(const Param& param)
	{
		param_ = param;
	}
	
	const QVector<QStringList>& TOPPASToolVertex::getOutputFileNames()
	{
		return output_file_names_;
	}
	
	void TOPPASToolVertex::updateOutputFileNames()
	{
		// recurse
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv)
			{
				tv->updateOutputFileNames();
			}
		}
		
		// now, parents are up-to-date --> update ourself:
		
		/*	first, determine (expected) number of elements of
				output parameters with list type (see HACK below) */
		QVector<IOInfo> in_params;
		int in_list_element_count = 0;
		getInputParameters(in_params);
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int param_index = (*it)->getTargetInParam();
			if (param_index < 0)
			{
				std::cerr << "Input parameter index out of bounds!" << std::endl;
				break;
			}
			
			if (in_params[param_index].param_name == "in")
			{
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
					in_list_element_count = (iflv->getFilenames()).count();
					break;
				}
				
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
				if (tv)
				{
					int out_param_index = (*it)->getSourceOutParam();
					if (out_param_index < 0)
					{
						std::cerr << "Output parameter index out of bounds!" << std::endl;
						break;
					}
					in_list_element_count = tv->output_file_names_[out_param_index].count();
					break;
				}
			}
		}
		
		// now, update the output file names:
		QVector<IOInfo> out_params;
		getOutputParameters(out_params);
		
		output_file_names_.clear();
		output_file_names_.resize(out_params.size());
		
		for (int i = 0; i < out_params.size(); ++i)
		{
			// search for an out edge for this parameter
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				int param_index = (*it)->getSourceOutParam();
				
				if (i == param_index)
				{
					// corresponding out edge found
					if (out_params[param_index].type == IOInfo::IOT_FILE)
					{
						output_file_names_[param_index].push_back((tmp_path_ + "TOPPAS_" + name_ + "_" + type_ + "_" + File::getUniqueName() + "__tmp.out").toQString());
					}
					else if (out_params[param_index].type == IOInfo::IOT_LIST)
					{
						/*	HACK: for output parameters with list type, we expect
								the number of elements to be the same as the number of
								elements of the "in" input list (thus, we also assume that there is
								an input parameter "in" with list type) */
						for (int j = 0; j < in_list_element_count; ++j)
						{
							output_file_names_[param_index].push_back((tmp_path_ + "TOPPAS_" + name_ + "_" + type_ + "_" + File::getUniqueName() + "__tmp.out").toQString());
						}
					}
					
					break; // we are done, don't push_back further file names for this parameter
				}
			}
		}
	}
	
	void TOPPASToolVertex::setStartedHere(bool b)
	{
		started_here_ = b;
	}

}

