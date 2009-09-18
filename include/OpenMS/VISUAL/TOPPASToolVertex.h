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

#ifndef OPENMS_VISUAL_TOPPASTOOLVERTEX_H
#define OPENMS_VISUAL_TOPPASTOOLVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtCore/QVector>

namespace OpenMS
{
	/**
		@brief A vertex representing a TOPP tool
		
		Besides TOPPASScene, this class contains most of the remaining functionality of
		TOPPAS regarding the execution of pipelines. Once a pipeline run is started
		from TOPPASScene, the execution is propagated from tool to tool and the
		TOPP tools are actually called from here.
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASToolVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			/// Stores the information for input/output files/lists
			struct IOInfo
			{
				///Standard constructor
				IOInfo()
					:	type(IOT_FILE),
						param_name(),
						valid_types(),
						listified(false)
				{
				}
				
				///Copy constructor
				IOInfo(const IOInfo& rhs)
					:	type(rhs.type),
						param_name(rhs.param_name),
						valid_types(rhs.valid_types),
						listified(rhs.listified)
				{
				}
				
				///The type
				enum IOType
				{
					IOT_FILE,
					IOT_LIST
				};
				
				///Comparison operator
				bool operator< (const IOInfo& rhs) const
				{
					if (type != rhs.type)
					{
						return type == IOT_FILE;
					}
					else
					{
						return param_name.compare(rhs.param_name) < 0;
					}
				}
				
				///Assignment operator
				IOInfo& operator= (const IOInfo& rhs)
				{
					type = rhs.type;
					param_name = rhs.param_name;
					valid_types = rhs.valid_types;
					listified = rhs.listified;
					
					return *this;
				}
				
				///The type of the parameter
				IOType type;
				///The name of the parameter
				String param_name;
				///The valid file types for this parameter
				StringList valid_types;
				///Is the parameter actually a single file parameter but is used in list iteration?
				bool listified;
			};
			
			/// Default constructor
			TOPPASToolVertex();
			/// Constructor
			TOPPASToolVertex(const String& name, const String& type = "", const String& tmp_path = "");
			/// Copy constructor
			TOPPASToolVertex(const TOPPASToolVertex& rhs);
			/// Destructor
			virtual ~TOPPASToolVertex();
			/// Assignment operator
			TOPPASToolVertex& operator= (const TOPPASToolVertex& rhs);
			
			/// Returns the name of the tool
			const String& getName();
			/// Returns the type of the tool
			const String& getType();
			/// Fills @p input_infos with the required input file/list parameters together with their valid types.
			void getInputParameters(QVector<IOInfo>& input_infos);
			/// Fills @p output_infos with the required output file/list parameters together with their valid types.
			void getOutputParameters(QVector<IOInfo>& output_infos);
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			/// Returns whether this node has already been processed during the current pipeline execution
			bool isFinished();
			/// Set whether this node has already been processed during the current pipeline execution
			void setFinished(bool b);
			/// Sets the Param object of this tool
			void setParam(const Param& param);
			/// Returns the Param object of this tool
			const Param& getParam();
			/// Starts the pipeline execution recursively
			void runRecursively();
			/// Checks if all parent nodes have finished the tool execution and, if so, runs the tool
			void runToolIfInputReady();
			/// Returns a vector containing the lists of output files for all output parameters
			const QVector<QStringList>& getOutputFileNames();
			/// Updates the vector containing the lists of output files for all output parameters
			void updateOutputFileNames();
			/// Sets whether the currently running pipeline has already been started at this vertex
			void setStartedHere(bool b);
			/// Sets the progress color
			void setProgressColor(const QColor& c);
			/// Returns the progress color
			QColor getProgressColor();
			/// Lets the user edit the parameters of the tool
			void editParam();
			/// Returns the number of iterations this tool has to perform
			int numIterations();
			/// Returns whether the list iteration mode is enabled
			bool listModeActive();
			/// (Un)sets the list iteration mode
			void setListModeActive(bool b);
			/// Returns the directory where this tool stores its output files
			String getOutputDir();
			/// Creates all necessary directories (called by the scene before the pipeline is run)
			void createDirs(const QString& out_dir);
			/// Sets the topological sort number and removes invalidated tmp files
			virtual void setTopoNr(UInt nr);
			
		public slots:
		
			/// Called when the execution of this tool has finished
			void executionFinished(int ec, QProcess::ExitStatus es);
			/// Called when the running TOPP tool produces output
			void forwardTOPPOutput();
			/// Called when the tool is started
			void toolStartedSlot();
			/// Called when the tool has finished
			void toolFinishedSlot();
			/// Called when the tool has crashed
			void toolCrashedSlot();
			/// Called when the tool has failed
			void toolFailedSlot();
			/// Called by an incoming edge when it has changed
			virtual void inEdgeHasChanged();
		
		signals:
		
			/// Emitted when the tool is started
			void toolStarted();
			/// Emitted when the tool is finished
			void toolFinished();
			/// Emitted when the tool crashes
			void toolCrashed();
			/// Emitted when the tool execution fails
			void toolFailed();
			/// Emitted from forwardTOPPOutput() to forward the signal outside
			void toppOutputReady(const QString& out);
		
		protected:
		
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
			
			/// Initializes the parameters with standard values (from -write_ini)
			void initParam_();
			/// Fills @p io_infos with the required input/output file/list parameters. If @p input_params is true, input params are returned, otherwise output params.
			void getParameters_(QVector<IOInfo>& io_infos, bool input_params);
			
			/// The name of the tool
			String name_;
			/// The type of the tool, or "" if it does not have a type
			String type_;
			/// The temporary path
			String tmp_path_;
			/// The parameters of the tool
			Param param_;
			/// Stores whether this node has already been processed during the current pipeline execution
			bool finished_;
			/// Stores whether the currently running pipeline has already been started at this vertex
			bool started_here_;
			/// Stores the file names of the different output parameters
			QVector<QStringList> output_file_names_;
			/// Color representing the progress (red = waiting, yellow = processing, green = finished, else: gray)
			QColor progress_color_;
			/// Stores whether we are currently in list mode
			bool list_mode_;
			/// The symbol for the list mode
			static QImage symbol_image_;
			/// The number of the current iteration
			int iteration_nr_;
			/// The length of (all) input lists
			int input_list_length_;
	};
}

#endif
