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
#include <QtCore/QRegExp>
#include <QtGui/QImage>

#include <QDesktopServices>
#include <QUrl>
#include <QMessageBox>
#include <QCoreApplication>

namespace OpenMS
{
	UInt TOPPASToolVertex::uid_ = 1;
	
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
		if (initParam_()) {}
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
		if (initParam_()) {}
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
	
	bool TOPPASToolVertex::initParam_(const QString& old_ini_file)
	{
		Param tmp_param;
		QString ini_file = QDir::tempPath() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
		if (type_ != "")
		{
			ini_file += type_.toQString() + "_";
		}
		ini_file += File::getUniqueName().toQString() + "_tmp.ini";
    ini_file = QDir::toNativeSeparators(ini_file);

		String call = name_ + " -write_ini " + ini_file;
		if (type_ != "")
		{
			call += " -type " + type_;
		}
		if (old_ini_file != "")
		{
			if (!File::exists(old_ini_file))
			{
				QMessageBox::critical(0,"Error",(String("Could not open '")+old_ini_file+"'!").c_str());
				return false;
			}
			call += " -ini " + String(old_ini_file);
		}
		
		if (system(call.c_str()) != 0)
		{
			QMessageBox::critical(0,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
			return false;
		}
		if(!File::exists(ini_file))
		{
			QMessageBox::critical(0,"Error",(String("Could not open '")+ini_file+"'!").c_str());
			return false;
		}
		
		tmp_param.load(String(ini_file).c_str());
		param_=tmp_param.copy(name_+":1:",true);
		//param_.remove("log");
		//param_.remove("no_progress");
		//param_.remove("debug");
		//// handled by TOPPAS anyway:
		//param_.remove("type");
		
		writeParam_(param_,ini_file);
		bool changed = false;
		if (old_ini_file != "")
		{
			//check if ini file has changed (quick & dirty by file size)
			QFile q_ini(ini_file);
			QFile q_old_ini(old_ini_file);
			changed = q_ini.size() != q_old_ini.size();
			QFile::remove(old_ini_file);
		}
		QFile::remove(ini_file);
		
		return changed;
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
		
		QVector<String> hidden_entries;
		// remove type (should not be edited)
		if (edit_param.exists("type"))
		{
			hidden_entries.push_back("type");
		}
		// remove entries that are handled by edges already, user should not see them
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
				hidden_entries.push_back(name);
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
				hidden_entries.push_back(name);
			}
		}
		
    // remove entries explained by edges
    foreach (const String& name, hidden_entries)
		{
			edit_param.remove(name);
		}

		TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, name_, type_, hidden_entries);
		if (dialog.exec())
		{
      // take new values
      param_.update(edit_param);
			
			
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
				io_info.param_name = it.getName();
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
					std::cerr << "TOPPAS: Unexpected parameter value!" << std::endl;
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
		painter->drawEllipse(46,-52, 14, 14);
		
		//topo sort number
		qreal x_pos = -64.0;
		qreal y_pos = -41.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
		
		if (progress_color_ != Qt::gray)
		{
			QString text;
			if (in_parameter_has_list_type_)
			{
				text = QString::number(iteration_nr_ == 1 ? input_list_length_ : 0);
				text += QString(" / ") + QString::number(input_list_length_); 
			}
			else
			{
				text = QString::number(iteration_nr_)+" / "+QString::number(num_iterations_);
			}
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
			painter->drawText((int)(62.0-text_boundings.width()), 48, text);
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
	
	void TOPPASToolVertex::runToolIfInputReady()
	{
		__DEBUG_BEGIN_METHOD__
		
		//check if everything ready
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv && !tv->isFinished())
			{
				// some tool that we depend on has not finished execution yet --> do not start yet
				debugOut_("Not run (parent not finished)");
				
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		// all inputs are ready --> GO!
		updateCurrentOutputFileNames();
		
		createDirs();
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QString ini_file = File::getTempDirectory().toQString()
						           + QDir::separator()
						           + getOutputDir().toQString()
						           + QDir::separator()
						           + name_.toQString();
		if (type_ != "")
		{
			ini_file += "_" + type_.toQString();
		}
		ini_file += ".ini";
		// do not write the ini yet - we might need to alter it
		
		QStringList shared_args;
		shared_args	<< "-no_progress";
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
		
    // we might need to modify input/output file parameters before storing to INI
    Param param_tmp = param_;

		for (int i = 0; i < num_iterations_; ++i)
		{
			debugOut_(String("Enqueueing process nr ")+i);
			QStringList args = shared_args;
			
			// add all input file parameters
			QVector<IOInfo> in_params;
			getInputParameters(in_params);
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				int param_index = (*it)->getTargetInParam();
				if (param_index < 0)
				{
					std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
					break;
				}

        String param_name = in_params[param_index].param_name;

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        if (param_name.hasPrefix("ETool:")) store_to_ini = true;
        if (!store_to_ini) args << "-" + param_name.toQString();
				
        QStringList file_list;

				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (tv)
				{
					int out_param_index = (*it)->getSourceOutParam();
					if (out_param_index < 0)
					{
						std::cerr << "TOPPAS: Output parameter index out of bounds!" << std::endl;
						break;
					}
					file_list = tv->current_output_files_[out_param_index];
				}
				else if (mv)
				{
					file_list = mv->getCurrentOutputList();
				}
				else if (iflv)
				{
          file_list = iflv->getInputFilenames();
				}
        QStringList p_files = getFileArgument_(file_list, i, in_parameter_has_list_type_);
        if (store_to_ini)
        {
          if (param_tmp.getValue(param_name).valueType()==DataValue::STRING_LIST) param_tmp.setValue(param_name, StringList(p_files));
          else 
          { 
            if (p_files.size()>1) throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            param_tmp.setValue(param_name, String(p_files[0]));
          }
        }
        else args << p_files;

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
            String param_name = out_params[param_index].param_name;

            bool store_to_ini = false;
            // check for GenericWrapper input/output files and put them in INI file:
            if (param_name.hasPrefix("ETool:")) store_to_ini = true;
            if (!store_to_ini) args << "-" + param_name.toQString();

						const QStringList& output_files = current_output_files_[param_index];

            QStringList p_files = getFileArgument_(output_files, i, in_parameter_has_list_type_);
            if (store_to_ini)
            {
              if (param_tmp.getValue(param_name).valueType()==DataValue::STRING_LIST) param_tmp.setValue(param_name, StringList(p_files));
              else 
              { 
                if (p_files.size()>1) throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
                param_tmp.setValue(param_name, String(p_files[0]));
              }
            }
            else args << p_files;
						
						break; // (regardless of the number of out edges, every argument must appear only once)
					}
				}
			}
      
      // each iteration might have different params (input/output items which are registered in subsections (GenericWrapper stuff))
      QString ini_file_iteration = QDir::toNativeSeparators( ini_file + QString::number(num_iterations_) );
      writeParam_(param_tmp, ini_file_iteration);
      args << "-ini" << ini_file_iteration;

			//create process
			QProcess* p = new QProcess();
			p->setProcessChannelMode(QProcess::MergedChannels);
			connect(p,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(executionFinished(int,QProcess::ExitStatus)));
			connect(p,SIGNAL(readyReadStandardOutput()),this,SLOT(forwardTOPPOutput()));
			connect(ts,SIGNAL(terminateCurrentPipeline()),p,SLOT(kill()));
			
			//enqueue process
			std::cout << "Enqueue: " << name_ << " \"" << String(args.join("\" \"")) << "\"" << std::endl;
			ts->enqueueProcess(p, name_.toQString(), args);
		}
		
		__DEBUG_END_METHOD__
	}
	
  QStringList TOPPASToolVertex::getFileArgument_(const QStringList& source_files, const int index, const bool as_list) const
  {
    if (as_list)
    {
      return source_files;
    }
    else
    {
      if (index >= source_files.size())
			{
        LOG_ERROR << "TOPPAS: Input list too short!" << std::endl;
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, source_files.size());
			}
			return QStringList(source_files[index]);
    }

  }

	void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
	{
		__DEBUG_BEGIN_METHOD__
		
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
			__DEBUG_END_METHOD__
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
			__DEBUG_END_METHOD__
			return;
		}
		
		++iteration_nr_;
		update(boundingRect());
		debugOut_(String("Increased iteration_nr_ to ")+iteration_nr_+" / "+num_iterations_);
		
		// notify the scene that this process has finished (so the next pending one can run)
		ts->runNextProcess();
		
		if (iteration_nr_ == num_iterations_) // all iterations performed --> proceed in pipeline
		{
			debugOut_("All iterations finished!");
			
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
				TOPPASVertex* tv = (*it)->getTargetVertex();
				debugOut_(String("Starting child ")+tv->getTopoNr());
				
				TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
				if (ttv)
				{
					ttv->runToolIfInputReady();
					continue;
				}
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(tv);
				if (mv)
				{
					mv->forwardPipelineExecution();
					continue;
				}
				TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
				if (oflv)
				{
					oflv->finish();
					continue;
				}
			}
			
			debugOut_("All children started!");
		}
		
		//clean up
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (p)
		{
			delete p;
		}
		
		__DEBUG_END_METHOD__
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
		bool found_in_parameter = false;
		getInputParameters(in_params);
		QStringList input_file_basenames;
		
		bool force = false;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int param_index = (*it)->getTargetInParam();
			if (param_index < 0)
			{
				std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
				break;
			}
			
			if (in_params[param_index].param_name == "in" || force)
			{
				found_in_parameter = true;
				in_parameter_has_list_type_ = (in_params[param_index].type == IOInfo::IOT_LIST);
				
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
          const QStringList& input_files = iflv->getInputFilenames();
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
						std::cerr << "TOPPAS: Output parameter index out of bounds!" << std::endl;
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
			
			//if last iteration and still no "in" parameter found, repeat last iteration and treat the edge as input parameter (dirty - TODO)
			if (it == inEdgesEnd()-1 && !found_in_parameter)
			{
				--it;
				force = true;
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
						QString f = File::getTempDirectory().toQString()
							          + QDir::separator()
							          + getOutputDir().toQString()
							          + QDir::separator()
                        + out_params[param_index].param_name.remove(':').toQString().left(50) // max 50 chars per subdir
							          + QDir::separator();
						if (f.length()>150) LOG_WARN << "Warning: the temporary path '" << String(f) << "' used in TOPPAS has many characters.\n"
                                         << "         TOPPAS might not be able to write files properly.\n";

						f += QString(input_file_basenames.first()
                         + "_to_"
                         + input_file_basenames.last()
                         + "_merged");
            f = f.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
            f += "_tmp" + QString::number(uid_++);
            f = QDir::toNativeSeparators(f);
						current_output_files_[param_index].push_back(f);
					}
					else
          // single input file
					{
						foreach (const QString& input_file, input_file_basenames)
						{
							QString f = File::getTempDirectory().toQString()
								          + QDir::separator()
								          + getOutputDir().toQString()
								          + QDir::separator()
								          + out_params[param_index].param_name.remove(':').toQString().left(50) // max 50 chars per subdir
								          + QDir::separator();
							if (f.length()>150) LOG_WARN << "Warning: the temporary path '" << String(f) << "' used in TOPPAS has many characters.\n"
                                           << "         TOPPAS might not be able to write files properly.\n";
              f += input_file;
							QRegExp rx("_tmp\\d+$"); // remove "_tmp<number>" if its a suffix
							int tmp_index = rx.indexIn(f);
              if (tmp_index != -1)
							{
								f = f.left(tmp_index);
							}
              f = f.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
              f += "_tmp" + QString::number(uid_++);
              f = QDir::toNativeSeparators(f);
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
		reset(true);
		
		TOPPASVertex::inEdgeHasChanged();
	}
	

	void TOPPASToolVertex::openContainingFolder()
	{
		QVector<IOInfo> out_infos;
		getOutputParameters(out_infos);
		
		if (out_infos.size() == all_written_output_files_.size())
		{
			foreach (const QStringList& files, all_written_output_files_)
			{
				if (files.size() > 0)
				{
          QFileInfo fi(files[0]);
          QString path = QFileInfo(fi.canonicalFilePath()).path();
          if (!QDir(path).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + path, QUrl::TolerantMode))))
          {
            QMessageBox::warning(0, "Open Folder Error", String("The folder " + path + " could not be opened!").toQString());
          }
        }
			}
		}
	}
	  
	
	String TOPPASToolVertex::getOutputDir()
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		String workflow_dir = File::removeExtension(File::basename(ts->getSaveFileName()));
		if (workflow_dir == "")
		{
			workflow_dir = "Untitled_workflow";
		}
		String dir = String("TOPPAS_tmp")+
			String(QDir::separator())+
			workflow_dir+
			String(QDir::separator())+
			get3CharsNumber_(topo_nr_)+"_"+getName();
		if (getType() != "")
		{
			dir += "_"+getType();
		}
		
		return dir;
	}
	
	void TOPPASToolVertex::createDirs()
	{
		QDir current_dir(File::getTempDirectory().toQString());
		
		if (!current_dir.mkpath(getOutputDir().toQString()))
		{
			std::cerr << "TOPPAS: Could not create path " << getOutputDir() << std::endl;
		}
		
		foreach (const QStringList& files, current_output_files_)
		{
			if (!files.isEmpty())
			{
				QString dir = File::path(files.first()).toQString();
				if (!File::exists(dir))
				{
					if (!current_dir.mkpath(dir))
					{
						std::cerr << "TOPPAS: Could not create path " << String(dir) << std::endl;
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
			reset(true);
			topo_nr_ = nr;
			emit somethingHasChanged();
		}
	}
	
	void TOPPASToolVertex::reset(bool reset_all_files)
	{
		__DEBUG_BEGIN_METHOD__
		
		finished_ = false;
		current_output_files_.clear();
		progress_color_ = Qt::gray;
		
		if (reset_all_files)
		{
			all_written_output_files_.clear();
			QString remove_dir = File::getTempDirectory().toQString() + QDir::separator() + getOutputDir().toQString();
			if (File::exists(remove_dir))
			{
				removeDirRecursively_(remove_dir);
			}
			// reset UID for tmp files
			uid_ = 1;
		}
		
		TOPPASVertex::reset(reset_all_files);
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASToolVertex::checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run)
	{
		__DEBUG_BEGIN_METHOD__
		
		//all parents checked?
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			if (!((*it)->getSourceVertex()->isScListLengthChecked()))
			{
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		// do all input lists have equal length?
		int parent_per_round = (*inEdgesBegin())->getSourceVertex()->getScFilesPerRound();
		int parent_total = (*inEdgesBegin())->getSourceVertex()->getScFilesTotal();		
		
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			if ((*it)->getSourceVertex()->getScFilesPerRound() != parent_per_round)
			{
				unequal_per_round.push_back(QString::number(topo_nr_));
				break;
			}
		}
		
		// n:1 tool? --> files per round = 1
		QVector<IOInfo> input_infos;
		getInputParameters(input_infos);
		QVector<IOInfo> output_infos;
		getOutputParameters(output_infos);
		bool in_param_list_type = false;
		bool out_param_file_type = false;
		foreach (const IOInfo& io, input_infos)
		{
			if (io.param_name == "in" && io.type == IOInfo::IOT_LIST)
			{
				in_param_list_type = true;
			}
		}
		foreach (const IOInfo& io, output_infos)
		{
			if (io.param_name == "out" && io.type == IOInfo::IOT_FILE)
			{
				out_param_file_type = true;
			}
		}
		
		if (in_param_list_type && out_param_file_type)
		{
			sc_files_per_round_ = 1;
			sc_files_total_ = parent_total / parent_per_round;
		}
		else
		{
			sc_files_per_round_ = parent_per_round;
			sc_files_total_ = parent_total;
		}
		
		sc_list_length_checked_ = true;
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			tv->checkListLengths(unequal_per_round, unequal_over_entire_run);
		}
		
		__DEBUG_END_METHOD__
	}

	bool TOPPASToolVertex::refreshParameters()
	{
		QString old_ini_file = QDir::tempPath() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
		if (type_ != "")
		{
			old_ini_file += type_.toQString() + "_";
		}
		old_ini_file += File::getUniqueName().toQString() + "_tmp_OLD.ini";
		writeParam_(param_,old_ini_file);
		
		bool changed = initParam_(old_ini_file);
		
		return changed;
	}
	
	void TOPPASToolVertex::writeParam_(const Param& param, const QString& ini_file)
	{
		Param save_param;
		save_param.setValue(name_+":1:toppas_dummy", DataValue("blub"));
		save_param.insert(name_+":1:", param);
		save_param.remove(name_+":1:toppas_dummy");
		save_param.setSectionDescription(name_+":1", "Instance '1' section for '"+name_+"'");
		save_param.store(ini_file);
	}
}

