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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{

	TransformationXMLFile::TransformationXMLFile()
		: XMLHandler("", "1.0"),
			XMLFile("/SCHEMAS/TrafoXML_1_0.xsd", "1.0"),
			trafo_(0),
			param_(),
			pairs_()
	{
	}

  void TransformationXMLFile::load(const String& filename, TransformationDescription& transformation)
  {
  	//Filename for error messages in XMLHandler
  	file_ = filename;

  	transformation.clear();
		trafo_ = &transformation;
  	param_.clear();
		pairs_.clear();
		
		parse_(filename,this);
	}
  					 
  void TransformationXMLFile::store(String filename, const TransformationDescription& transformation)
  {
		if ( transformation.getName() == "" )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "will not write a transformation with empty name");
		}
		
  	//open stream
		std::ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		os.precision(writtenDigits<DoubleReal>());
		
		//write header
		os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
		//add XSLT file if it can be found
		os << "<TrafoXML version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/TrafoXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << "\n";
				
		// open tag
		os << "\t<Transformation name=\"" << transformation.getName() << "\">\n";

		// write parameters
		for ( Param::ParamIterator it = transformation.getParameters().begin(); it != transformation.getParameters().end(); ++it )
		{
			if(it->value.valueType()!=DataValue::EMPTY_VALUE)
			{
				switch( it->value.valueType() )
				{
				case DataValue::INT_VALUE:
					os << "\t\t<Param  type=\"int\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
					break;
				case DataValue::DOUBLE_VALUE:
					os << "\t\t<Param  type=\"float\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
					break;
				case DataValue::STRING_VALUE:
				case DataValue::STRING_LIST:
				case DataValue::INT_LIST:
				case DataValue::DOUBLE_LIST:
					os << "\t\t<Param  type=\"string\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
					break;
				default: // no other value types are supported!
					fatalError(STORE, String("Unsupported parameter type of parameter '") + it->name + "' with value '" + it->value.toString() +"'");
					break;
				};
			}
		}
		
		//write pairs
		Size pairs_size = transformation.getPairs().size();
		if (pairs_size!=0)
		{
			os << "\t\t<Pairs count=\"" << pairs_size << "\">\n";
			for (Size i=0; i<pairs_size; ++i)
			{
				os << "\t\t\t<Pair from=\"" << transformation.getPairs()[i].first << "\" to=\"" << transformation.getPairs()[i].second << "\"/>\n";
			}
			os << "\t\t</Pairs>\n";
		}
		
		// close tag
		os << "\t</Transformation>\n";

		//write footer
		os << "</TrafoXML>" << "\n";
		
		//close stream
		os.close();
 }

  
	void TransformationXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{		

		String element = sm_.convert(qname);
	
		if (element == "TrafoXML")
		{
			//check file version against schema version
			DoubleReal file_version = attributeAsDouble_(attributes,"version");
			if (file_version>version_.toDouble())
			{
				warning(LOAD, String("The XML file (") + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
			}
		}
		else if ( element == "Transformation" )
		{
			trafo_->setName(attributeAsString_(attributes,"name"));
		}
		else if ( element == "Param" )
		{
			String type = attributeAsString_(attributes,"type");
			if ( type == "int" )
			{
				param_.setValue(attributeAsString_(attributes,"name"),attributeAsInt_(attributes,"value"));
			}
			else if ( type == "float" )
			{
				param_.setValue(attributeAsString_(attributes,"name"),attributeAsDouble_(attributes,"value"));
			}
			else if ( type == "string" )
			{
				param_.setValue(attributeAsString_(attributes,"name"),String(attributeAsString_(attributes,"value")));
			}
			else
			{
			  error(LOAD, String("Unsupported parameter type: '")+type+"'");
			}

		}
		else if ( element == "Pairs" )
		{
			pairs_.reserve(attributeAsInt_(attributes,"count"));
		}
		else if ( element == "Pair" )
		{
			pairs_.push_back(make_pair((Real)attributeAsDouble_(attributes,"from"),(Real)attributeAsDouble_(attributes,"to")));
		}
		else
		{
		  warning(LOAD,String("Unknown element: '")+element+"'");
		}
	}
	
	void TransformationXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String element = sm_.convert(qname);
	
		if ( element == "Transformation" )
		{
			trafo_->setParameters(param_);
			trafo_->setPairs(pairs_);
		}
	}

} // namespace OpenMS
