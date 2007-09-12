// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

namespace OpenMS
{	
	namespace Internal
	{
		//--------------------------------------------------------------------------------
	  void FeatureXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << "</" << sm_.convert(qname) << ">\n";
				if (String(sm_.convert(qname)) != enum2str_(TAGMAP,DESCRIPTION))
				{
					return;
				}
			}
	
			int tag = leaveTag(qname);
	
			// Do something depending on the tag
			switch(tag) {
				case DESCRIPTION:
					// delegate control to ExperimentalSettings handler
					{
						// initialize parser
						xercesc::XMLPlatformUtils::Initialize();
						xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
						parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
						parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
						
						MzDataExpSettHandler handler( *((ExperimentalSettings*)map_),file_);
						handler.resetErrors();
						parser->setContentHandler(&handler);
						parser->setErrorHandler(&handler);
						
						String tmp(exp_sett_.str().c_str());
						
						xercesc::MemBufInputSource source((const XMLByte*)(tmp.c_str()), tmp.size(), "dummy");
		      	parser->parse(source);
		      	delete(parser);
					}
					break;
				case FEATURE:
					if ((!options_.hasRTRange() || options_.getRTRange().encloses(feature_->getPosition()[0]))
							&&	(!options_.hasMZRange() || options_.getMZRange().encloses(feature_->getPosition()[1]))
							&&	(!options_.hasIntensityRange() || options_.getIntensityRange().encloses(feature_->getIntensity())))
					{
						map_->push_back(*feature_);
					}
					delete feature_;
					break;
				case FEATMODEL:
					model_desc_->setParam(*param_);
					feature_->setModelDescription(*model_desc_);
					delete param_;
					delete model_desc_;
					break;
				case HULLPOINT:
					current_chull_->addPoint(*hull_position_);
					delete hull_position_;
					break;
				case CONVEXHULL:
					feature_->getConvexHulls().push_back(*current_chull_);
					delete current_chull_;
					break;
			}
	  }
	
	  void FeatureXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << '<' << sm_.convert(qname);
				UInt n=attributes.getLength();
				for (UInt i=0; i<n; ++i)
				{
					exp_sett_ << ' ' << sm_.convert(attributes.getQName(i)) << "=\""	<< sm_.convert(attributes.getValue(i)) << '\"';
				}
				exp_sett_ << '>';
				return;
			}
	
			int tag = enterTag(qname, attributes);
	
			String tmp_str;
			// Do something depending on the tag
			switch(tag) {
				case DESCRIPTION: 
					exp_sett_ << '<' << sm_.convert(qname) << '>'; 
					break;
				case FEATURELIST:
					if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
					break;
	   		case FEATURE: 	 
	   			feature_        = new Feature();
	   			break;
				case QUALITY:
					tmp_str = getAttributeAsString_(DIM, true, qname);
					current_qcoord_ = asUInt_(tmp_str); 
					break;
				case POSITION:
					tmp_str = getAttributeAsString_(DIM, true, qname);
					current_pcoord_ = asUInt_(tmp_str); 
					break;
				case CONVEXHULL: 
					current_chull_  = new ConvexHullType(); 
					break;
				case HULLPOINT:  
					hull_position_  = new DPosition<2>(); 
					break;
				case HPOSITION:  
					tmp_str = getAttributeAsString_(DIM, true, qname);
					current_hcoord_ = asUInt_(tmp_str); 
					break;
				case FEATMODEL:
					model_desc_ = new ModelDescription<2>();
			  	param_ = new Param();
					tmp_str = getAttributeAsString_(NAME, true, qname);
			  	if (tmp_str != "")
			  	{
			  		model_desc_->setName(tmp_str);
			  	}
			  	break;
			  case PARAM:
			  {
			  	String name = getAttributeAsString_(NAME, true, qname);
					String value = getAttributeAsString_(VALUE, true, qname);
			  	if (name != "" && value != "")
			  		param_->setValue(name, value);
			  	break;
			  }
			}
		}
	
	  void FeatureXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
	  {
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << sm_.convert(chars);
				return;
			}
	
			// find the tag that the parser is in right now
	 		for (UInt i=0; i<is_parser_in_tag_.size(); i++)
	 		{
				if (is_parser_in_tag_[i])
				{
					switch(i) 
					{
						case FEATINTENSITY: 
							feature_->setIntensity(asDouble_(sm_.convert(chars))); 
							break;
						case POSITION:
							feature_->getPosition()[current_pcoord_] = asDouble_(sm_.convert(chars));
							break;
						case QUALITY:       
							feature_->setQuality(current_qcoord_, asDouble_(sm_.convert(chars)));
								break;
						case OVERALLQUALITY:  
							feature_->setOverallQuality(asDouble_(sm_.convert(chars))); break;
						case CHARGE:          
							feature_->setCharge(asInt_(sm_.convert(chars)));
							break;
						case HPOSITION:       
							(*hull_position_)[current_hcoord_] = asDouble_(sm_.convert(chars)); 
							break;
						case META:						
							feature_->setMetaValue(3,String(sm_.convert(chars))); 
							break;
					}
				}
			}
	  }
	
	
	 	void FeatureXMLHandler::writeTo(std::ostream& os)
		{
			UniqueIdGenerator id_generator = UniqueIdGenerator::instance();
	
			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
			   << "<featureMap xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
	
			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cmap_),"");
			handler.writeTo(os);
	
			os << "\t<featureList count=\"" << cmap_->size() << "\">\n";
	
			// write features with their corresponding attributes
			for (UInt s=0; s<cmap_->size(); s++)
			{
				const Feature& feat = (*cmap_)[s];
	
				os << "\t\t<feature id=\"" << id_generator.getUID() << "\">" << std::endl;
	
				for (UInt i=0; i<2;i++)
				{
					os <<	"\t\t\t<position dim=\"" << i << "\">" << feat.getPosition()[i] << "</position>" << 	std::endl;
				}
	
				os << "\t\t\t<intensity>" << feat.getIntensity() << "</intensity>" << std::endl;
	
				for (UInt i=0; i<2;i++)
				os << "\t\t\t<quality dim=\"" << i << "\">" << feat.getQuality(i) << "</quality>" << std:: endl;
	
				if(feat.getMetaValue(3)!=DataValue::EMPTY)
					os << "\t\t\t<meta>" << feat.getMetaValue(3) << "</meta>" << std:: endl;
	
				os << "\t\t\t<overallquality>" << feat.getOverallQuality() << "</overallquality>" << std:: endl;
				os << "\t\t\t<charge>" << feat.getCharge() << "</charge>" << std:: endl;
	
				// write model description
				ModelDescription<2> desc = feat.getModelDescription();
				os << "\t\t\t<model name=\"" << desc.getName() << "\">" << std:: endl;
				Param modelp = desc.getParam();
				Param::ParamIterator piter = modelp.begin();
				while (piter != modelp.end())
				{
					os << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\">";
					os << "</param>" << std::endl;
					piter++;
				}
				os << "\t\t\t</model>" << std::endl;
	
				// write convex hull
				ConvexHullVector hulls = feat.getConvexHulls();
				ConvexHullVector::iterator citer = hulls.begin();
	
				UInt hulls_count = hulls.size();
	
				for (UInt i=0;i<hulls_count; i++)
				{
					os << "\t\t\t<convexhull nr=\"" << i << "\">" << std:: endl;
	
					ConvexHullType current_hull = hulls[i];
					UInt hull_size       = current_hull.getPoints().size();
	
					for (UInt j=0;j<hull_size;j++)
					{
						os << "\t\t\t\t<hullpoint>" << std::endl;
	
						DPosition<2> pos = current_hull.getPoints()[j];
						UInt pos_size = pos.size();
						for (UInt k=0; k<pos_size; k++)
						{
							os << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << std::endl;
						}
	
						os << "\t\t\t\t</hullpoint>" << std::endl;
					} // end for (..hull_size..)
	
					os << "\t\t\t</convexhull>" << std::endl;
				} // end  for ( ... hull_count..)
	
				os << "\t\t</feature>\n";
	
			} // end for ( features )
	
			os << "\t</featureList>\n</featureMap>\n";
			os <<
				"<!-- Local Variables: -->\n"
				"<!-- mode: nxml -->\n"
				"<!-- tab-width: 2 -->\n"
				"<!-- End: -->\n";
		}
		
	} //namespace Internal

} // namespace OpenMS






