// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>

namespace OpenMS
{
  namespace Internal
  {
  	
    void ConsensusXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
    	//std::cout << "END: " << sm_.convert(qname) << std::endl;
    	static XMLCh* s_consensuselement = xercesc::XMLString::transcode("consensusElement");
      
      if (equal_(qname,s_consensuselement))
      {
				consensus_map_->push_back(act_cons_element_);
      }
    }

    void ConsensusXMLHandler::characters(const XMLCh* const /*chars*/, unsigned int /*length*/)
  	{
  	}

    void ConsensusXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      //std::cout << "BEGIN: " << sm_.convert(qname) << std::endl;
      static XMLCh* s_consensuselement = xercesc::XMLString::transcode("consensusElement");
      static XMLCh* s_name = xercesc::XMLString::transcode("name");
      static XMLCh* s_map = xercesc::XMLString::transcode("map");
      static XMLCh* s_element = xercesc::XMLString::transcode("element");
      static XMLCh* s_range = xercesc::XMLString::transcode("range");
      static XMLCh* s_centroid = xercesc::XMLString::transcode("centroid");
      static XMLCh* s_rt = xercesc::XMLString::transcode("rt");
      static XMLCh* s_mz = xercesc::XMLString::transcode("mz");
      static XMLCh* s_it = xercesc::XMLString::transcode("it");
      static XMLCh* s_id = xercesc::XMLString::transcode("id");
      static XMLCh* s_rtmin = xercesc::XMLString::transcode("rtMin");
      static XMLCh* s_rtmax = xercesc::XMLString::transcode("rtMax");
      static XMLCh* s_mzmin = xercesc::XMLString::transcode("mzMin");
      static XMLCh* s_mzmax = xercesc::XMLString::transcode("mzMax");
      static XMLCh* s_itmin = xercesc::XMLString::transcode("itMin");
      static XMLCh* s_itmax = xercesc::XMLString::transcode("itMax");
      static XMLCh* s_consensusxml = xercesc::XMLString::transcode("consensusXML");
      	
      String tmp_str;
      if (equal_(qname,s_map))
    	{
      	consensus_map_->setFileName(attributeAsInt_(attributes,s_id),attributeAsString_(attributes,s_name));
    	}
      else if (equal_(qname,s_consensuselement))
    	{
        act_cons_element_ = ConsensusFeature();
        consensus_element_range_ = true;
    	}
      else if (equal_(qname,s_centroid))
    	{
          tmp_str = attributeAsString_(attributes,s_rt);
          if (tmp_str != "")
          {
            pos_[Peak2D::RT] = asDouble_(tmp_str);
          }

          tmp_str = attributeAsString_(attributes,s_mz);
          if (tmp_str != "")
          {
            pos_[Peak2D::MZ] = asDouble_(tmp_str);
          }

          tmp_str = attributeAsString_(attributes,s_it);
          if (tmp_str != "")
          {
            it_ = asDouble_(tmp_str);
          }

    	}
      else if (equal_(qname,s_range))
    	{
        if (consensus_element_range_)
        {
          tmp_str = attributeAsString_(attributes, s_rtmin);
          if (tmp_str != "")
          {
            pos_range_.setMinX(asDouble_(tmp_str));

            tmp_str = attributeAsString_(attributes, s_rtmax);
            if (tmp_str != "")
            {
              pos_range_.setMaxX(asDouble_(tmp_str));

              tmp_str = attributeAsString_(attributes, s_mzmin);
              if (tmp_str != "")
              {
                pos_range_.setMinY(asDouble_(tmp_str));

                tmp_str = attributeAsString_(attributes, s_mzmax);
                if (tmp_str != "")
                {
                  pos_range_.setMaxY(asDouble_(tmp_str));

                  tmp_str = attributeAsString_(attributes, s_itmin);
                  if (tmp_str != "")
                  {
                    it_range_.setMin(asDouble_(tmp_str));

                    tmp_str = attributeAsString_(attributes, s_itmax);
                    if (tmp_str != "")
                    {
                      it_range_.setMax(asDouble_(tmp_str));

                      consensus_element_range_ = false;
                    }
                  }
                }
              }
            }
          }
        }
    	}
      else if (equal_(qname,s_element))
    	{
        FeatureHandle act_index_tuple;
        tmp_str = attributeAsString_(attributes, s_map);
        if (tmp_str != "")
        {
          UInt map_index = asUInt_(tmp_str);
          tmp_str = attributeAsString_(attributes, s_id);

          if (tmp_str != "")
          {
            UInt element_index = asUInt_(tmp_str);

            act_index_tuple.setMapIndex(map_index);
            act_index_tuple.setElementIndex(element_index);

            tmp_str = attributeAsString_(attributes, s_rt);
            DPosition<2> pos;
            pos[0] = asDouble_(tmp_str);
            tmp_str = attributeAsString_(attributes, s_mz);
            pos[1] = asDouble_(tmp_str);

            act_index_tuple.setPosition(pos);
            act_index_tuple.setIntensity(attributeAsDouble_(attributes,s_it));
            act_cons_element_.insert(act_index_tuple,false);
          }
        }
        act_cons_element_.getPosition() = pos_;
        act_cons_element_.getPositionRange() = pos_range_;
        act_cons_element_.setIntensity(it_);
        act_cons_element_.getIntensityRange() = it_range_;
    	}
      else if (equal_(qname,s_consensusxml))
    	{
	   		//check file version against schema version
				String file_version="1.0";
				optionalAttributeAsString_(file_version,attributes,"version");
				if (file_version.toDouble()>version_.toDouble())
				{
					warning("The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
			}
    }

    void ConsensusXMLHandler::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      << "<consensusXML version=\"" << version_ << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/ConsensusXML_1_2.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

      const Map<UInt,String>& name_vector = consensus_map_->getFileNames();
      os << "\t<mapList count=\"" << name_vector.size() << "\">\n";
      for (Map<UInt,String>::const_iterator it=name_vector.begin(); it!=name_vector.end(); ++it)
      {
        os << "\t\t<map id=\"" << it->first << "\" name =\"" << it->second << "\"/>\n";
      }
      os << "\t</mapList>\n";

      os << "\t<consensusElementList>\n";
      for (UInt i = 0; i < consensus_map_->size(); ++i)
      {
      	const ConsensusFeature& elem = consensus_map_->operator[](i);
        os << "\t\t<consensusElement id=\""<< i << "\">\n";
        os << "\t\t\t<centroid rt=\"" << elem.getRT()
        << "\" mz=\"" << elem.getMZ()
        << "\" it=\"" << elem.getIntensity() <<"\"/>\n";
        os << "\t\t\t<range rtMin=\"" << elem.getPositionRange().min()[0]
        << "\" rtMax=\"" << elem.getPositionRange().max()[0]
        << "\" mzMin=\"" << elem.getPositionRange().min()[1]
        << "\" mzMax=\"" << elem.getPositionRange().max()[1]
        << "\" itMin=\"" << elem.getIntensityRange().min()
        << "\" itMax=\"" << elem.getIntensityRange().max() <<"\"/>\n";

        os << "\t\t\t<groupedElementList>\n";
        for (ConsensusFeature::HandleSetType::const_iterator it = elem.begin(); it != elem.end(); ++it)
        {
          os  << "\t\t\t\t<element id=\"" << it->getElementIndex()
          << "\" map=\"" << it->getMapIndex()
          << "\" rt=\"" << it->getPosition()[0]
          << "\" mz=\"" << it->getPosition()[1]
          << "\" it=\"" << it->getIntensity() << "\"/>\n";
        }
        os << "\t\t\t</groupedElementList>\n";
        os << "\t\t</consensusElement>\n";
      }
      os << "\t</consensusElementList>\n";
      os << "</consensusXML>"<< std::endl;
    }

  } // namespace Internal
} // namespace OpenMS






