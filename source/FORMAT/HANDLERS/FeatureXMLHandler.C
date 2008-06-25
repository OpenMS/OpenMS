// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

using namespace std;

namespace OpenMS
{	
	namespace Internal
	{
	  void FeatureXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			static const XMLCh* s_dim = xercesc::XMLString::transcode("dim");
			static const XMLCh* s_name = xercesc::XMLString::transcode("name");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_type = xercesc::XMLString::transcode("type");
			static const XMLCh* s_int = xercesc::XMLString::transcode("int");
			static const XMLCh* s_float = xercesc::XMLString::transcode("float");
			static const XMLCh* s_string = xercesc::XMLString::transcode("string");	
			static const XMLCh* s_featuremap = xercesc::XMLString::transcode("featureMap");

			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//cout << "Start: " << tag << endl;
			
			// collect Experimental Settings
			if (in_description_)	
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
			
			if (tag=="feature")
			{
				feature_ = Feature();
			}
			else if (tag=="description")
			{
				//cout << "RT range : " << options_.getRTRange() << endl;
				//cout << "MZ range : " << options_.getMZRange() << endl;
				//cout << "Int range: " << options_.getIntensityRange() << endl;
				
				exp_sett_ << "<description>"; 
				in_description_ = true;
			}
			else if (tag=="featureList")
			{
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			else if (tag=="quality" || tag=="hposition" || tag=="position")
			{
				dim_ = attributeAsInt_(attributes,s_dim); 
			}
			else if (tag=="convexhull")
			{
				current_chull_ = ConvexHull2D(); 
			}
			else if (tag=="hullpoint")
			{
				hull_position_ = DPosition<2>::zero; 
			}
			else if (tag=="model")
			{
				model_desc_ = new ModelDescription<2>();
		  	param_.clear();
		  	model_desc_->setName(attributeAsString_(attributes,s_name));
			}
			else if (tag=="param")
			{
				String name = attributeAsString_(attributes,s_name);
				String value = attributeAsString_(attributes,s_value);
				if (name != "" && value != "") param_.setValue(name, value);
			}
			else if (tag == "userParam")
			{				
				const XMLCh* value = attributes.getValue(s_value);
				const XMLCh* type = attributes.getValue(s_type);
				String name = sm_.convert(attributes.getValue(s_name));
				
				if(*type==*s_int)
				{
					feature_.setMetaValue(name, xercesc::XMLString::parseInt(value));
				}
				else if (*type==*s_float)
				{
					feature_.setMetaValue(name, atof(sm_.convert(value)) );
				}
				else if (*type==*s_string)
				{
					feature_.setMetaValue(name, (String)sm_.convert(value));
				}
				else
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invalid userParam type '") + sm_.convert(type) + "'" );
				}
			}
			else if (equal_(qname,s_featuremap))
			{
				//check file version against schema version
				String file_version="1.0";
				optionalAttributeAsString_(file_version,attributes,"version");
				if (file_version.toDouble()>version_.toDouble())
				{
					warning(String("The XML file (") + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
			}	
		}
		
	  void FeatureXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static const XMLCh* s_description = xercesc::XMLString::transcode("description");
			static const XMLCh* s_feature = xercesc::XMLString::transcode("feature");
			static const XMLCh* s_model = xercesc::XMLString::transcode("model");
			static const XMLCh* s_hullpoint = xercesc::XMLString::transcode("hullpoint");
			static const XMLCh* s_convexhull = xercesc::XMLString::transcode("convexhull");
		
			//cout << "End: '" << sm_.convert(qname) <<"' - '"<< sm_.convert(s_description) << "'" << endl;
			
			open_tags_.pop_back();
			
			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << "</" << sm_.convert(qname) << ">\n";
				if (!equal_(qname,s_description)) return;
			}
			
			if (equal_(qname,s_description))
			{
				in_description_ = false;
				// call MzDataExpSett parser
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
			else if (equal_(qname,s_feature))
			{
				if ((!options_.hasRTRange() || options_.getRTRange().encloses(feature_.getRT()))
					&&	(!options_.hasMZRange() || options_.getMZRange().encloses(feature_.getMZ()))
					&&	(!options_.hasIntensityRange() || options_.getIntensityRange().encloses(feature_.getIntensity())))
				{
					map_->push_back(feature_);
				}
			}
			else if (equal_(qname,s_model))
			{
				model_desc_->setParam(param_);
				feature_.setModelDescription(*model_desc_);
				delete model_desc_;
			}
			else if (equal_(qname,s_hullpoint))
			{
				current_chull_.addPoint(hull_position_);
			}
			else if (equal_(qname,s_convexhull))
			{
				feature_.getConvexHulls().push_back(current_chull_);
			}
	  }
	
	  void FeatureXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
	  {
	  	// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << sm_.convert(chars);
				return;
			}
			
			String& current_tag = open_tags_.back();
			if (current_tag == "intensity")
			{
				feature_.setIntensity(asDouble_(sm_.convert(chars))); 
			}
			else if (current_tag == "position")
			{
				feature_.getPosition()[dim_] = asDouble_(sm_.convert(chars));
			}
			else if (current_tag == "quality")
			{
				feature_.setQuality(dim_, asDouble_(sm_.convert(chars)));
			}
			else if (current_tag == "overallquality")
			{
				feature_.setOverallQuality(asDouble_(sm_.convert(chars)));
			}
			else if (current_tag == "charge")
			{
				feature_.setCharge(asInt_(chars));
			}
			else if (current_tag == "hposition")
			{
				hull_position_[dim_] = asDouble_(sm_.convert(chars)); 
			}
	  }
	
	
	 	void FeatureXMLHandler::writeTo(ostream& os)
		{
			UInt identifier = 0;

			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
			   << "<featureMap version=\"" << version_ << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_2.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
	
			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cmap_),"");
			handler.writeTo(os);
	
			os << "\t<featureList count=\"" << cmap_->size() << "\">\n";
			
			// write features with their corresponding attributes
			for (UInt s=0; s<cmap_->size(); s++)
			{
				const Feature& feat = (*cmap_)[s];
				os << "\t\t<feature id=\"f_" << (identifier++) << "\">" << endl;
				for (UInt i=0; i<2;i++)
				{
					os <<	"\t\t\t<position dim=\"" << i << "\">" << feat.getPosition()[i] << "</position>" << 	endl;
				}
				os << "\t\t\t<intensity>" << feat.getIntensity() << "</intensity>" << endl;
				for (UInt i=0; i<2;i++)
				{
					os << "\t\t\t<quality dim=\"" << i << "\">" << feat.getQuality(i) << "</quality>" <<  endl;
				}
				os << "\t\t\t<overallquality>" << feat.getOverallQuality() << "</overallquality>" <<  endl;
				os << "\t\t\t<charge>" << feat.getCharge() << "</charge>" <<  endl;
	
				// write model description
				ModelDescription<2> desc = feat.getModelDescription();
				os << "\t\t\t<model name=\"" << desc.getName() << "\">" <<  endl;
				Param modelp = desc.getParam();
				Param::ParamIterator piter = modelp.begin();
				while (piter != modelp.end())
				{
					os << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\">";
					os << "</param>" << endl;
					piter++;
				}
				os << "\t\t\t</model>" << endl;
	
				// write convex hull
				vector<ConvexHull2D> hulls = feat.getConvexHulls();
				vector<ConvexHull2D>::iterator citer = hulls.begin();
	
				UInt hulls_count = hulls.size();
	
				for (UInt i=0;i<hulls_count; i++)
				{
					os << "\t\t\t<convexhull nr=\"" << i << "\">" <<  endl;
	
					ConvexHull2D current_hull = hulls[i];
					UInt hull_size  = current_hull.getPoints().size();
	
					for (UInt j=0;j<hull_size;j++)
					{
						os << "\t\t\t\t<hullpoint>" << endl;
	
						DPosition<2> pos = current_hull.getPoints()[j];
						UInt pos_size = pos.size();
						for (UInt k=0; k<pos_size; k++)
						{
							os << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << endl;
						}
	
						os << "\t\t\t\t</hullpoint>" << endl;
					} // end for (..hull_size..)
	
					os << "\t\t\t</convexhull>" << endl;
				} // end  for ( ... hull_count..)
				
				writeUserParam_("userParam", os, feat, 3);
				
				os << "\t\t</feature>\n";
	
			} // end for ( features )
	
			os << "\t</featureList>\n</featureMap>\n";
		}

		void FeatureXMLHandler::setOptions(const PeakFileOptions& options)
		{ 
			options_ = options; 
		}

		
	} //namespace Internal
} // namespace OpenMS
