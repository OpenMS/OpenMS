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
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

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
			progress_color_(Qt::gray)
		{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		if (initParam_()) {}
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type)
		: TOPPASVertex(),
			name_(name),
			type_(type),
			param_(),
			progress_color_(Qt::gray)
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
			param_(rhs.param_),
			progress_color_(rhs.progress_color_)
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
		finished_ = rhs.finished_;
		progress_color_ = rhs.progress_color_;
		
		return *this;
	}
	
	bool TOPPASToolVertex::initParam_(const QString& old_ini_file)
	{
		Param tmp_param;
    // this is the only exception for writing directly to the tmpDir, instead of a subdir of tmpDir, as scene()->getTempDir() might not be available yet
		QString ini_file = File::getTempDirectory().toQString() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
		if (type_ != "")
		{
			ini_file += type_.toQString() + "_";
		}
		ini_file += File::getUniqueName().toQString() + "_tmp.ini";
    ini_file = QDir::toNativeSeparators(ini_file);

		String call = String("\"") + File::getExecutablePath() + name_ + "\"" + " -write_ini " + ini_file;
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
			QMessageBox::critical(0, "Error", (String("Could not execute '") + call + "'!\n\nMake sure the TOPP tools are present in '" + File::getExecutablePath() + "', that you have permission to write to the temporary file path, and that there is space left in the temporary file path.").c_str());
			return false;
		}
		if (!File::exists(ini_file))
		{
			QMessageBox::critical(0, "Error", (String("Could not open '")+ini_file+"'!").c_str());
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

		//topo sort number
		qreal x_pos = -64.0;
		qreal y_pos = -41.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
		
		if (progress_color_ != Qt::gray)
		{
			QString text = QString::number(round_counter_) + " / " + QString::number(round_total_);
			
      QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
			painter->drawText((int)(62.0-text_boundings.width()), 48, text);
		}

		// recycling status
    if (this->allow_output_recycling_)
    {
      painter->setPen(Qt::green);
      painter->drawChord(-7,-52, 14, 14, 0*16, 180*16);
    }
		
		// progress light
		painter->setPen(Qt::black);
		painter->setBrush(progress_color_);
		painter->drawEllipse(46,-52, 14, 14);
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
	
	const String& TOPPASToolVertex::getName() const
	{
		return name_;
	}
	
	const String& TOPPASToolVertex::getType() const
	{
		return type_;
	}
	
	void TOPPASToolVertex::run()
	{
		__DEBUG_BEGIN_METHOD__
		
		//check if everything ready (there might be more than one upstream node - ALL need to be ready)
		if (!isUpstreamReady())	return;
		
    if (finished_)
    {
      std::cerr << "This should not happen. Calling an already finished node '" << this->name_ << "' (#" << this->getTopoNr() << ")!\n";
      throw Exception::IllegalSelfOperation(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());

    QString ini_file = ts->getTempDir()
						           + QDir::separator()
						           + getOutputDir().toQString()
						           + QDir::separator()
						           + name_.toQString();
		if (type_ != "") ini_file += "_" + type_.toQString();
		// do not write the ini yet - we might need to alter it
		
    RoundPackages pkg;
    String error_msg("");
    bool success = buildRoundPackages(pkg, error_msg);
    if (!success)
    {
      std::cerr << "Could not retrieve input files from upstream nodes...\n";
      emit toolFailed(error_msg.toQString());
      return;
    }
		
		// all inputs are ready --> GO!
		if (!updateCurrentOutputFileNames(pkg, error_msg)) // based on input, we prepare output names
    {
      emit toolFailed(error_msg.toQString());
      return;
    }
		
		createDirs();

    emit toolStarted();
		
    /// update round status
    round_total_ = (int) pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0;        // once round_counter_ reaches round_total_, we are done

		QStringList shared_args;
		shared_args	<< "-no_progress";
		if (type_ != "") shared_args << "-type" << type_.toQString();

		// get *all* input|output file parameters (regardless if edge exists)
		QVector<IOInfo> in_params, out_params;
		getInputParameters(in_params);
    getOutputParameters(out_params);

    bool ini_round_dependent = false; // indicates if we need a new ini file for each round (usually GenericWrapper issue)

    for (int round = 0; round < round_total_; ++round)
		{
			debugOut_(String("Enqueueing process nr ") + round + "/" + round_total_);
			QStringList args = shared_args;

      // we might need to modify input/output file parameters before storing to INI
      Param param_tmp = param_;

			/// INCOMING EDGES
      for (RoundPackageConstIt ite = pkg[round].begin();
                               ite!= pkg[round].end();
                               ++ite)
			{
        TOPPASEdge incoming_edge = *(ite->second.edge);

        int param_index = incoming_edge.getTargetInParam();
				if (param_index < 0 || param_index >= in_params.size())
				{
					std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
					return;
				}

        String param_name = in_params[param_index].param_name;

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        if (param_name.hasPrefix("ETool:"))
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }
        if (!store_to_ini) args << "-" + param_name.toQString();
				
        QStringList file_list = ite->second.filenames;

        if (store_to_ini)
        {
          if (param_tmp.getValue(param_name).valueType()==DataValue::STRING_LIST) param_tmp.setValue(param_name, StringList(file_list));
          else 
          { 
            if (file_list.size()>1) throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            param_tmp.setValue(param_name, String(file_list[0]));
          }
        }
        else args << file_list;

			}

			// OUTGOING EDGES
      // ;output names are already prepared by 'updateCurrentOutputFileNames()'
      typedef RoundPackage::iterator EdgeIndexIt;
      for (EdgeIndexIt it_edge  = output_files_[round].begin();
                       it_edge != output_files_[round].end();
                       ++ it_edge)
			{
        int param_index = it_edge->first;
        String param_name = out_params[param_index].param_name;

        bool store_to_ini = false;
        // check for GenericWrapper input/output files and put them in INI file:
        if (param_name.hasPrefix("ETool:"))
        {
          store_to_ini = true;
          ini_round_dependent = true;
        }
        if (!store_to_ini) args << "-" + param_name.toQString();

        const QStringList& output_files = output_files_[round][param_index].filenames;

        if (store_to_ini)
        {
          if (param_tmp.getValue(param_name).valueType()==DataValue::STRING_LIST) param_tmp.setValue(param_name, StringList(output_files));
          else 
          { 
            if (output_files.size()>1) throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Multiple files were given to a param which supports only single files! ('" + param_name + "')");
            param_tmp.setValue(param_name, String(output_files[0]));
          }
        }
        else args << output_files;
			}
      
      // each iteration might have different params (input/output items which are registered in subsections (GenericWrapper stuff))
      QString ini_file_iteration;
      if (ini_round_dependent)
      {
        ini_file_iteration = QDir::toNativeSeparators( ini_file + QString::number(round) + ".ini" );
      }
      else
      {
        ini_file_iteration = QDir::toNativeSeparators( ini_file + ".ini" );
      }
      writeParam_(param_tmp, ini_file_iteration);
      args << "-ini" << ini_file_iteration;

			//create process
			QProcess* p;
      if (!ts->isDryRun())
      {
        p = new QProcess();
      }
      else
      {
        p = new FakeProcess();
      }

			p->setProcessChannelMode(QProcess::MergedChannels);
			connect(p, SIGNAL(readyReadStandardOutput()), this, SLOT(forwardTOPPOutput()));
			connect(ts, SIGNAL(terminateCurrentPipeline()), p, SLOT(kill()));
      // let this node know that round is done
			connect(p, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(executionFinished(int,QProcess::ExitStatus)));
      
			//enqueue process
			std::cout << "Enqueue: " << File::getExecutablePath() + name_ << " \"" << String(args.join("\" \"")) << "\"" << std::endl;
			ts->enqueueProcess(p, (File::getExecutablePath() + name_).toQString(), args);
		}

    // run pending processes
		ts->runNextProcess();

		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
	{
		__DEBUG_BEGIN_METHOD__
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());

    //** ERROR handling
		if (es != QProcess::NormalExit)
		{
			emit toolCrashed();
		}
		else if (ec != 0)
		{
			emit toolFailed();
		}
    else
    //** no error ... proceed
    {  
      ++round_counter_;
      std::cout << (String("Increased iteration_nr_ to ") + round_counter_ + " / " + round_total_ ) << " for " << this->name_ << std::endl;
		
      if (round_counter_ == round_total_) // all iterations performed --> proceed in pipeline
		  {
			  debugOut_("All iterations finished!");
			
        if (finished_)
        {
          std::cout << "SOMETHING is very fishy. The vertex is already set to finished, yet there was still a thread spawning...\n";
        }
        if (!ts->isDryRun()) renameOutput_(); // rename generated files by content
			  finished_ = true;
			  emit toolFinished();
			
			  // call all childs, proceed in pipeline
			  for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			  {
				  TOPPASVertex* tv = (*it)->getTargetVertex();
				  debugOut_(String("Starting child ") + tv->getTopoNr());
					tv->run();
			  }
			
			  debugOut_("All children started!");
		  }
    }
		
		//clean up
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (p)
		{
			delete p;
		}

    ts->processFinished();

		__DEBUG_END_METHOD__
	}
	
  bool TOPPASToolVertex::renameOutput_()
  {
    // get all output names
    QStringList files = this->getFileNames();

    std::set<String> unique;
    std::map<String, String> name_old_to_new;

    // create mapping from old to new filenames, while ensuring that they are unique
    foreach (QString file, files)
    {
      QFileInfo fi(file);
      String new_suffix = FileHandler::typeToName( FileHandler::getTypeByContent(file) );
      String new_prefix = String(fi.path() + "/" + fi.baseName()) + ".";
      String new_name = new_prefix + new_suffix;
      if (unique.count(new_name))
      { // make a new name
        Int counter(0);
        while (unique.count(new_prefix + counter + "." + new_suffix)) ++counter;
        new_name = new_prefix + counter + "." + new_suffix;
      }
      
      // filename is unique - use it
      unique.insert(new_name);
      name_old_to_new[file] = new_name;
    }

    for (Size i=0; i<output_files_.size(); ++i)
    {
      for (RoundPackageIt it = output_files_[i].begin();
                          it!= output_files_[i].end();
                          ++it)
      {
        for (int fi=0; fi < it->second.filenames.size(); ++fi)
        {
          // rename file and update record
          QFile file(it->second.filenames[fi]);
          if (File::exists(name_old_to_new[it->second.filenames[fi]]))
          {
            bool success = File::remove(name_old_to_new[it->second.filenames[fi]]);
            if (!success)
            {
              std::cerr << "Could not remove " << name_old_to_new[it->second.filenames[fi]] << "\n";
              return false;
            }
          }
          bool success = file.rename(name_old_to_new[it->second.filenames[fi]].toQString());
          if (!success)
          {
            std::cerr << "Could not rename " << String(it->second.filenames[fi]) << " to " << name_old_to_new[it->second.filenames[fi]] << "\n";
            return false;
          }
          it->second.filenames[fi] = name_old_to_new[it->second.filenames[fi]].toQString();
        }
      }
    }
    return true;
  }

	const Param& TOPPASToolVertex::getParam()
	{
		return param_;
	}
	
	void TOPPASToolVertex::setParam(const Param& param)
	{
		param_ = param;
	}
	
	bool TOPPASToolVertex::updateCurrentOutputFileNames(const RoundPackages& pkg, String& error_msg)
	{
    if (pkg.size() < 1)
    {
      error_msg = "Less than one round received from upstream tools. Something is fishy!\n";
      std::cerr << error_msg;
      return false;
    }

    // look for the input with the most files in round 0 (as this is the maximal number of output files we can produce)
    // we assume the number of files is equal in all rounds...
    // however, we upstream nodes which use 'recycling' of input, as the names will always be the same
    int max_size_index = -1;
    int max_size = -1;
    for (RoundPackageConstIt it  = pkg[0].begin();
                             it != pkg[0].end();
                             ++it)
    {
      if (it->second.filenames.size() > max_size && !it->second.edge->getSourceVertex()->isRecyclingEnabled())
      {
        max_size_index = it->first;
        max_size       = it->second.filenames.size();
      }
    }

    if (max_size_index == -1)
    {
      error_msg = "Did not find upstream nodes with unrecycled names. Something is fishy!\n";
      std::cerr << error_msg;
      return false;
    }

    // use input names from the selected upstream vertex (hoping that this is the maximal number of files we are going to produce)
    std::vector<QStringList> per_round_basenames;
    for (int i = 0; i < pkg.size(); ++i)
    {
      per_round_basenames.push_back( pkg[i].find(max_size_index)->second.filenames );
    }

		// now, update the output file names (each output parameter gets its name from the same input parameter, i.e. the longest one):
		QVector<IOInfo> out_params;
		getOutputParameters(out_params);
		
    // clear output file list
		output_files_.clear();
    output_files_.resize(pkg.size()); // #rounds
		
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		for (int i = 0; i < out_params.size(); ++i)
		{
			// search for an out edge for this parameter (not required to exist)
      bool found(false);
      int param_index;
      TOPPASEdge* param_edge;
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				param_index = (*it)->getSourceOutParam();
				if (i == param_index) // corresponding out edge found
				{
          param_edge = *it;
          found=true;
          break;
        }
      }
      if (!found) continue;

      // store edge for this param for all rounds
      for (Size r=0; r<per_round_basenames.size(); ++r)
      {
        VertexRoundPackage vrp;
        vrp.edge = param_edge;
        output_files_[r][param_index] = vrp; // index by index of source-out param
      }

			QString f = ts->getTempDir()
							    + QDir::separator()
							    + getOutputDir().toQString()
							    + QDir::separator()
                  + out_params[param_index].param_name.remove(':').toQString().left(50) // max 50 chars per subdir
							    + QDir::separator();
			if (f.length()>150) LOG_WARN << "Warning: the temporary path '" << String(f) << "' used in TOPPAS has many characters.\n"
                                   << "         TOPPAS might not be able to write files properly.\n";

      for (Size r=0; r<per_round_basenames.size(); ++r)
      {
        QString fn = f;
          
			  // check if tool consumes list and outputs single file (such as IDMerger or FileMerger)
			  if (per_round_basenames[r].size()>1 && out_params[param_index].type == IOInfo::IOT_FILE)
			  {
				  fn += QString( QFileInfo(per_round_basenames[r].first()).fileName()
                        + "_to_"
                        + QFileInfo(per_round_basenames[r].last()).fileName()
                        + "_merged");
          fn = fn.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
          fn += "_tmp" + QString::number(uid_++);
          fn = QDir::toNativeSeparators(fn);
          output_files_[r][param_index].filenames.push_back(fn);
        }
        else // each input file will have a corresponding output file
        {
          foreach (const QString& input_file, per_round_basenames[r])
				  {
            fn += QFileInfo(input_file).fileName(); // discard directory
					  QRegExp rx("_tmp\\d+$"); // remove "_tmp<number>" if its a suffix
					  int tmp_index = rx.indexIn(fn);
            if (tmp_index != -1) fn = fn.left(tmp_index);
            fn = fn.left(220); // allow max of 220 chars per path+filename (~NTFS limit)
            fn += "_tmp" + QString::number(uid_++);
            fn = QDir::toNativeSeparators(fn);
					  output_files_[r][param_index].filenames.push_back(fn);
          }
        }
			}
		}
    return true;
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
    QString path = getFullOutputDirectory().toQString();
    if (!QDir(path).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + path, QUrl::TolerantMode))))
    {
      QMessageBox::warning(0, "Open Folder Error", "The folder " + path + " could not be opened!");
    }
	}
	  
	String TOPPASToolVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    return QDir::toNativeSeparators( ts->getTempDir() + QDir::separator() + getOutputDir().toQString() );
  }
	
	String TOPPASToolVertex::getOutputDir() const
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
		QDir dir;
    bool ok = dir.mkpath(getFullOutputDirectory().toQString());

		if (!ok)
		{
			std::cerr << "TOPPAS: Could not create path " << getFullOutputDirectory() << std::endl;
		}
		
    // subsdirectories named after the output parameter name
    QStringList files = this->getFileNames();
		foreach (const QString& file, files)
		{
			QString sdir = File::path(file).toQString();
			if (!File::exists(sdir))
			{
				if (!dir.mkpath(sdir))
				{
					std::cerr << "TOPPAS: Could not create path " << String(sdir) << std::endl;
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
		output_files_.clear();
		progress_color_ = Qt::gray;
		
		if (reset_all_files)
		{
			QString remove_dir = getFullOutputDirectory().toQString();
			if (File::exists(remove_dir))
			{
				File::removeDirRecursively(remove_dir);
			}
			// reset UID for tmp files
			uid_ = 1;
		}
		
		TOPPASVertex::reset(reset_all_files);
		
		__DEBUG_END_METHOD__
	}
	

	bool TOPPASToolVertex::refreshParameters()
	{
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QString old_ini_file = ts->getTempDir() + QDir::separator() + "TOPPAS_" + name_.toQString() + "_";
		if (type_ != "")
		{
			old_ini_file += type_.toQString() + "_";
		}
		old_ini_file += File::getUniqueName().toQString() + "_tmp_OLD.ini";
		writeParam_(param_, old_ini_file);
		
		bool changed = initParam_(old_ini_file);
		QFile::remove(old_ini_file);

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

