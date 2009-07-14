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
				enum IOType
				{
					IOT_FILE,
					IOT_LIST
				};
				
				IOType type;
				String param_name;
				StringList valid_types;
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
			/// Runs the tool
			void compute();
			/// Returns whether this node has already been processed during the current pipeline execution
			bool isFinished();
			/// Set whether this node has already been processed during the current pipeline execution
			void setFinished(bool b);
			/// Sets the Param object of this tool
			void setParam(const Param& param);
			/// Returns the Param object of this tool
			const Param& getParam();
			
		protected:
		
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
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
			/// Stores the file names of the different output parameters
			QVector<QStringList> output_file_names_;
			
	};
}

#endif
