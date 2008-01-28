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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FeaturePairsHandler.h>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    void FeaturePairsHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
			static const XMLCh* s_dim = xercesc::XMLString::transcode("dim");
			static const XMLCh* s_name = xercesc::XMLString::transcode("name");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_type = xercesc::XMLString::transcode("type");
			static const XMLCh* s_int = xercesc::XMLString::transcode("int");
			static const XMLCh* s_float = xercesc::XMLString::transcode("float");
			static const XMLCh* s_string = xercesc::XMLString::transcode("string");	
							
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			if (tag=="pair")
			{
				pair_ = ElementPair<Feature>();
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
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invlid userParam type '") + sm_.convert(type) + "'" );
				}
			}
    }

    void FeaturePairsHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
			static const XMLCh* s_first = xercesc::XMLString::transcode("first");
			static const XMLCh* s_second = xercesc::XMLString::transcode("second");
			static const XMLCh* s_pair = xercesc::XMLString::transcode("pair");
			static const XMLCh* s_model = xercesc::XMLString::transcode("model");
			static const XMLCh* s_hullpoint = xercesc::XMLString::transcode("hullpoint");
			static const XMLCh* s_convexhull = xercesc::XMLString::transcode("convexhull");
			
			//cout << "End: '" << sm_.convert(qname) <<"' - '"<< sm_.convert(s_description) << "'" << endl;
			
			open_tags_.pop_back();
			
			if (equal_(qname,s_first))
			{
				pair_.setFirst(feature_);
				feature_ = Feature();
			}
			else if (equal_(qname,s_second))
			{
				pair_.setSecond(feature_);
				feature_ = Feature();
			}
			else if (equal_(qname,s_pair))
			{
				pairs_->push_back(pair_);
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

    void FeaturePairsHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
    {
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
			else if (current_tag == "pairquality")
			{
				pair_.setQuality(asDouble_(sm_.convert(chars)));
			}
		}

    void FeaturePairsHandler::writeTo(ostream& os)
    {

      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
      os << "<featurePairs xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeaturePairsXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl;

      // write features with their attributes
      for (UInt s=0; s<cpairs_->size(); s++)
      {
        const ElementPair< Feature >& pair = (*cpairs_)[s];

        os << "\t<pair nr=\"" << s << "\">" << endl;
        os << "\t\t<pairquality>" << pair.getQuality() << "</pairquality>" << endl;

        os << "\t\t<first>" << endl;
        Feature first = pair.getFirst();
        writeFeature_(os,first);
        os << "\t\t</first>" << endl;

        os << "\t\t<second>" << endl;
        Feature seco  = pair.getSecond();
        writeFeature_(os,seco);
        os << "\t\t</second>" << endl;

        os << "\t</pair>" << endl;

      } // end for ( features )

      os << "</featurePairs>" << endl;
      os <<
      "<!-- Local Variables: -->\n"
      "<!-- mode: nxml -->\n"
      "<!-- tab-width: 2 -->\n"
      "<!-- End: -->\n"
      ;
    }

    void FeaturePairsHandler::writeFeature_(ostream& os, Feature dfeat)
    {
      os << "\t\t<feature id=\"" << id_generator_.getUID() << "\">" << endl;

      Feature::PositionType pos = dfeat.getPosition();
      UInt dpos_size = pos.size();

      for (UInt i=0; i<dpos_size;i++)
      {
        os << "\t\t\t<position dim=\"" << i << "\">" << pos[i] << "</position>" <<  endl;
      }

      os << "\t\t\t<intensity>" << dfeat.getIntensity() << "</intensity>" << endl;

      for (UInt i=0; i<dpos_size;i++)
      {
        os << "\t\t\t<quality dim=\"" << i << "\">" << dfeat.getQuality(i) << "</quality>" <<  endl;
      }

      os << "\t\t\t<overallquality>" << dfeat.getOverallQuality() << "</overallquality>" <<  endl;
      os << "\t\t\t<charge>" << dfeat.getCharge() << "</charge>" <<  endl;

      // write model description
      ModelDescription<2> desc = dfeat.getModelDescription();
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
      vector<ConvexHull2D> hulls = dfeat.getConvexHulls();
      vector<ConvexHull2D>::iterator citer = hulls.begin();

      UInt hulls_count = hulls.size();

      for (UInt i=0;i<hulls_count; i++)
      {
        os << "\t\t\t<convexhull nr=\"" << i << "\">" << endl;

        ConvexHull2D current_hull = hulls[i];
        UInt hull_size = current_hull.getPoints().size();

        for (UInt j=0;j<hull_size;j++)
        {
          os << "\t\t\t\t<hullpoint>" << endl;

          Feature::PositionType pos = current_hull.getPoints()[j];
          UInt pos_size = pos.size();
          for (UInt k=0; k<pos_size; k++)
          {
            os << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << endl;
          }

          os << "\t\t\t\t</hullpoint>" << endl;
        } // end for (..hull_size..)

        os << "\t\t\t</convexhull>" << endl;
      } // end  for ( ... hull_count..)
			writeUserParam_("userParam", os, dfeat, 3);
      os << "\t\t</feature>\n";
    }
    
  } // namespace Internal
} // namespace OpenMS
