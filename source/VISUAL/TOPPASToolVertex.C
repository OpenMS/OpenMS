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
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtGui/QImage>

namespace OpenMS
{	
	TOPPASToolVertex::TOPPASToolVertex()
		:	TOPPASVertex(),
			name_(),
			type_(),
			param_(),
			finished_(false),
			progress_color_(Qt::gray),
			iteration_nr_(0),
			input_list_length_(1)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type, const String& tmp_path)
		: TOPPASVertex(),
			name_(name),
			type_(type),
			tmp_path_(tmp_path),
			param_(),
			finished_(false),
			progress_color_(Qt::gray),
			iteration_nr_(0),
			input_list_length_(1)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		initParam_();
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs)
		:	TOPPASVertex(rhs),
			name_(rhs.name_),
			type_(rhs.type_),
			tmp_path_(rhs.tmp_path_),
			param_(rhs.param_),
			finished_(rhs.finished_),
			progress_color_(rhs.progress_color_),
			iteration_nr_(rhs.iteration_nr_),
			input_list_length_(rhs.input_list_length_)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
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
		progress_color_ = rhs.progress_color_;
		iteration_nr_ = rhs.iteration_nr_;
		input_list_length_ = rhs.input_list_length_;
		
		return *this;
	}
	
	void TOPPASToolVertex::initParam_()
	{
		Param tmp_param;
		QString ini_file = QDir::tempPath() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
		if (type_ != "")
		{
			ini_file += type_.toQString() + "_";
		}
		ini_file += File::getUniqueName().toQString() + "_tmp.ini";
		
		String call = name_ + " -write_ini " + ini_file;
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
		
		tmp_param.load(String(ini_file).c_str());
		param_=tmp_param.copy(name_+":1:",true);
		param_.remove("log");
		param_.remove("no_progress");
		param_.remove("debug");
		
		// handled by TOPPAS anyway:
		param_.remove("type");
		
		// remove tmp ini file
		QFile::remove(ini_file);
	}
	
	void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		editParam();
	}
	
	void TOPPASToolVertex::editParam()
	{
		QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
		String default_dir = "";
		
		// use a copy for editing
		Param edit_param(param_);
		
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
			if (edit_param.exists(name))
			{
				hidden_entries.push_back(edit_param.getEntry(name));
				edit_param.remove(name);
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
			if (edit_param.exists(name))
			{
				hidden_entries.push_back(edit_param.getEntry(name));
				edit_param.remove(name);
			}
		}
		
		TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, name_, type_);
		if (dialog.exec())
		{
			param_ = edit_param;
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
			
			progress_color_ = Qt::gray;
			emit somethingHasChanged();
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
		//painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform);
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
			pen.setColor(Qt::darkBlue);
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
		
		QPainterPath path;
		path.addRect(-70.0, -60.0, 140.0, 120.0);		
 		painter->drawPath(path);
 		
 		pen.setColor(pen_color_);
 		painter->setPen(pen);
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
		
		// progress light
		painter->setPen(Qt::black);
		painter->setBrush(progress_color_);
		painter->drawEllipse(45,-52, 14, 14);
		
		//topo sort number
		qreal x_pos = -62.0;
		qreal y_pos = 48.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
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
		updateCurrentOutputFileNames();
		createDirs();
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QString ini_file = ts->getOutDir()
							+QDir::separator()
							+getOutputDir().toQString()
							+QDir::separator()
							+name_.toQString();
		if (type_ != "")
		{
			ini_file += "_"+type_.toQString();
		}
		ini_file += ".ini";
		
		Param save_param;
		save_param.setValue(name_+":1:toppas_dummy", DataValue("blub"));
		save_param.insert(name_+":1:", param_);
		save_param.remove(name_+":1:toppas_dummy");
		save_param.setSectionDescription(name_+":1", "Instance '1' section for '"+name_+"'");
		save_param.store(ini_file);
		
		QStringList shared_args;
		shared_args	<< "-ini"
								<< ini_file
								<< "-no_progress";
		if (type_ != "")
		{
			shared_args << "-type" << type_.toQString();
		}
		
		ts->setPipelineRunning(true);
		emit toolStarted();
		
		/* Decide whether or not to iterate over the single files of the incoming lists
		 * depending on the type of the input parameter "-in" (file or list). */
		num_iterations_ = in_parameter_has_list_type_ ? 1 : input_list_length_;
		iteration_nr_ = 0; // needed in executionFinished()
		
		for (int i = 0; i < num_iterations_; ++i)
		{
			QStringList args = shared_args;
			
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
					const QStringList& source_out_files = tv->current_output_files_[out_param_index];
					
					if (in_parameter_has_list_type_)
					{
						args << source_out_files;
					}
					else
					{
						args << source_out_files[i];
					}
					continue;
				}
				
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
				if (mv)
				{
					const QStringList& input_files = mv->getCurrentOutputList();

					if (in_parameter_has_list_type_)
					{
						args << input_files;
					}
					else
					{
						args << input_files[i];
					}
					continue;
				}

				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
					const QStringList& input_files = iflv->getFilenames();
					if (in_parameter_has_list_type_)
					{
						args << input_files;
					}
					else
					{
						args << input_files[i];
					}
					continue;
				}
			}
			
			// add all output file parameters
			QVector<IOInfo> out_params;
			getOutputParameters(out_params);
			
			for (int j = 0; j < out_params.size(); ++j)
			{
				// search for an out edge for this parameter
				for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
				{
					int param_index = (*it)->getSourceOutParam();
					
					if (j == param_index)
					{
						args << "-"+(out_params[param_index].param_name).toQString();
						const QStringList& output_files = current_output_files_[param_index];
						if (in_parameter_has_list_type_)
						{
							args << output_files;
						}
						else
						{
							args << output_files[i];
						}
						
						break; // (regardless of the number of out edges, every argument must appear only once)
					}
				}
			}
			
			//start process
			QProcess* p = new QProcess();
			p->setProcessChannelMode(QProcess::MergedChannels);
			connect(p,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(executionFinished(int,QProcess::ExitStatus)));
			connect(p,SIGNAL(readyReadStandardOutput()),this,SLOT(forwardTOPPOutput()));
			connect(ts,SIGNAL(terminateCurrentPipeline()),p,SLOT(kill()));
			
			//start process
			p->start(name_.toQString(), args);
			p->waitForStarted();
		}
	}
	
	void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
	{
		iteration_nr_++;
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		
		if (es != QProcess::NormalExit)
		{
			ts->setPipelineRunning(false);
			emit toolCrashed();
			//clean up
			QProcess* p = qobject_cast<QProcess*>(QObject::sender());
			if (p)
			{
				delete p;
			}
			return;
		}
		
		if (ec != 0)
		{
			ts->setPipelineRunning(false);
			emit toolFailed();
			//clean up
			QProcess* p = qobject_cast<QProcess*>(QObject::sender());
			if (p)
			{
				delete p;
			}
			return;
		}
		
		if (iteration_nr_ == num_iterations_) // all iterations performed --> proceed in pipeline
		{
			finished_ = true;
			emit toolFinished();
			
			if (all_written_output_files_.size() != current_output_files_.size())
			{
				all_written_output_files_.resize(current_output_files_.size());
			}
			for (int i = 0; i < current_output_files_.size(); ++i)
			{
				all_written_output_files_[i] << current_output_files_[i];
			}
			
			// notify all childs that we are finished, proceed in pipeline
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getTargetVertex());
				if (tv)
				{
					tv->runToolIfInputReady();
					continue;
				}
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getTargetVertex());
				if (mv)
				{
					mv->forwardPipelineExecution();
					continue;
				}
				TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>((*it)->getTargetVertex());
				if (oflv)
				{
					oflv->finish();
					continue;
				}
			}
		}
		
		/*	Normally, the subtree finished signal is propagated upstream from finished output nodes.
				However, if there is a blocking "merge all" node in the way, this will not happen
				--> additionally check here if subtree is finished ("merge all" nodes will return true
				in every case)	*/
		checkIfSubtreeFinished();
		
		/*	Because "merge all" nodes wait for the signal that all round-based mergers above them have completed
				their merging rounds, we additionally check this here (in case there is no merger that sends this signal
				downwards	*/ 
		checkIfAllUpstreamMergersFinished();
		
		//clean up
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (p)
		{
			delete p;
		}
	}
	
	bool TOPPASToolVertex::isFinished()
	{
		return finished_;
	}
	
	const Param& TOPPASToolVertex::getParam()
	{
		return param_;
	}
	
	void TOPPASToolVertex::setParam(const Param& param)
	{
		param_ = param;
	}
	
	const QVector<QStringList>& TOPPASToolVertex::getCurrentOutputFileNames()
	{
		return current_output_files_;
	}
	
	const QVector<QStringList>& TOPPASToolVertex::getAllWrittenOutputFileNames()
	{
		return all_written_output_files_;
	}
	
	void TOPPASToolVertex::updateCurrentOutputFileNames()
	{
		QVector<IOInfo> in_params;
		input_list_length_ = 1; // stays like that if -in param is not a list
		getInputParameters(in_params);
		QStringList input_file_basenames;
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		
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
				in_parameter_has_list_type_ = (in_params[param_index].type == IOInfo::IOT_LIST);
				
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
					const QStringList& input_files = iflv->getFilenames();
					input_list_length_ = input_files.count();
					foreach (const QString& str, input_files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
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
					const QStringList& input_files = tv->current_output_files_[out_param_index];
					input_list_length_ = input_files.count();
					foreach (const QString& str, input_files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
					break;
				}

				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
				if (mv)
				{
					const QStringList& files = mv->getCurrentOutputList();
					input_list_length_ = files.count();

					foreach (const QString& str, files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
					break;
				}
			}
		}
		
		// now, update the output file names:
		QVector<IOInfo> out_params;
		getOutputParameters(out_params);
		
		current_output_files_.clear();
		current_output_files_.resize(out_params.size());
		
		for (int i = 0; i < out_params.size(); ++i)
		{
			// search for an out edge for this parameter
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				int param_index = (*it)->getSourceOutParam();
				
				if (i == param_index) // corresponding out edge found
				{
					// check if tool consumes list and outputs single file (such as IDMerger or FileMerger)
					if (in_parameter_has_list_type_ && out_params[param_index].type == IOInfo::IOT_FILE)
					{
						QString f = ts->getOutDir()
							+QDir::separator()
							+getOutputDir().toQString()
							+QDir::separator()
							+out_params[param_index].param_name.toQString()
							+QDir::separator()
							+input_file_basenames.first()
							+"_to_"
							+input_file_basenames.last()
							+"_merged_tmp";
						current_output_files_[param_index].push_back(f);
					}
					else
					{
						foreach (const QString& str, input_file_basenames)
						{
							QString f = ts->getOutDir()
								+QDir::separator()
								+getOutputDir().toQString()
								+QDir::separator()
								+out_params[param_index].param_name.toQString()
								+QDir::separator()
								+str;
							if (!f.endsWith("_tmp"))
							{
								f += "_tmp";
							}
							current_output_files_[param_index].push_back(f);
						}
					}

					break; // we are done, don't push_back further file names for this parameter
				}
			}
		}
	}

	void TOPPASToolVertex::forwardTOPPOutput()
	{
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (!p)
		{
			return;
		}
		
		QString out = p->readAllStandardOutput();
		emit toppOutputReady(out);
	}
	
	void TOPPASToolVertex::setProgressColor(const QColor& c)
	{
		progress_color_ = c;
	}
	
	QColor TOPPASToolVertex::getProgressColor()
	{
		return progress_color_;
	}
	
	void TOPPASToolVertex::toolStartedSlot()
	{
		progress_color_ = Qt::yellow;
		update(boundingRect());
	}
	
	void TOPPASToolVertex::toolFinishedSlot()
	{
		progress_color_ = Qt::green;
		update(boundingRect());
	}
	
	void TOPPASToolVertex::toolFailedSlot()
	{
		progress_color_ = Qt::red;
		update(boundingRect());
	}

	void TOPPASToolVertex::toolCrashedSlot()
	{
		progress_color_ = Qt::red;
		update(boundingRect());
	}

	void TOPPASToolVertex::inEdgeHasChanged()
	{
		// something has changed --> tmp files might be invalid --> reset
		reset(true,true);
		
		TOPPASVertex::inEdgeHasChanged();
	}
	
	void TOPPASToolVertex::openInTOPPView()
	{
		QVector<IOInfo> out_infos;
		getOutputParameters(out_infos);
		
		if (out_infos.size() == all_written_output_files_.size())
		{
			foreach (const QStringList& files, all_written_output_files_)
			{
				if (files.size() > 0)
				{
					QProcess* p = new QProcess();
					p->setProcessChannelMode(QProcess::ForwardedChannels);
					p->start("TOPPView", files);
				}
			}
		}
	}
	
	String TOPPASToolVertex::getOutputDir()
	{
		String dir = String("TOPPAS_tmp")+String(QDir::separator())+get3CharsNumber_(topo_nr_)+"_"+getName();
		if (getType() != "")
		{
			dir += "_"+getType();
		}
		
		return dir;
	}
	
	void TOPPASToolVertex::createDirs()
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QDir current_dir(ts->getOutDir());
		
		foreach (const QStringList& files, current_output_files_)
		{
			if (!files.isEmpty())
			{
				QString dir = File::path(files.first()).toQString();
				if (!File::exists(dir))
				{
					if (!current_dir.mkpath(dir))
					{
						std::cerr << "Could not create path " << String(dir) << std::endl;
					}
				}
			}
		}
	}
	
	void TOPPASToolVertex::setTopoNr(UInt nr)
	{
		if (topo_nr_ != nr)
		{
			// topological number changes --> output dir changes --> reset
			reset(true,true);
			topo_nr_ = nr;
			emit somethingHasChanged();
		}
	}
	
	void TOPPASToolVertex::reset(bool reset_all_files, bool mergers_finished)
	{
		TOPPASVertex::reset(reset_all_files, mergers_finished);
		finished_ = false;
		current_output_files_.clear();
		progress_color_ = Qt::gray;
		update(boundingRect());
		
		if (reset_all_files)
		{
			all_written_output_files_.clear();
			QString remove_dir = qobject_cast<TOPPASScene*>(scene())->getOutDir() + QDir::separator() + getOutputDir().toQString();
			if (File::exists(remove_dir))
			{
				removeDirRecursively_(remove_dir);
			}
		}
	}
	
}

