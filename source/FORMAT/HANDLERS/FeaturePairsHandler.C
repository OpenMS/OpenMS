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


namespace OpenMS
{
  namespace Internal
  {

    void FeaturePairsHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      int tag = enterTag(qname, attributes);

      String tmp_str;
      switch(tag)
      {
      case FEATURE:  feature_        = new Feature(); break;
      case PAIR:     pair_           = new ElementPair <Feature> (); break;
      case QUALITY:
        tmp_str = getAttributeAsString_(DIM, true, qname);
        current_qcoord_ = asUInt_(tmp_str);
        break;
      case POSITION:
        tmp_str = getAttributeAsString_(DIM, true, qname);
        current_pcoord_ = asUInt_(tmp_str);
        break;
      case CONVEXHULL: current_chull_  = new ConvexHull2D(); break;
      case HULLPOINT:  hull_position_  = new Feature::PositionType(); break;
      case HPOSITION:
        tmp_str = getAttributeAsString_(DIM, true, qname);
        current_hcoord_ = asUInt_(tmp_str);
        break;
      case FEATMODEL:
        model_desc_ = new ModelDescription<2>();
        param_ = new Param();
        tmp_str = getAttributeAsString_(NAME, true, qname);
        if (tmp_str != "")
          model_desc_->setName(tmp_str);
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

    void FeaturePairsHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      int tag = leaveTag(qname);
      switch(tag)
      {
      case FIRST:
        pair_->setFirst(*feature_);
        delete feature_;
        break;
      case SECOND:
        pair_->setSecond(*feature_);
        delete feature_;
        break;
      case PAIR:
        pairs_->push_back(*pair_);
        delete pair_;
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

    void FeaturePairsHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
    {
      for (UInt i=0; i<is_parser_in_tag_.size(); i++)
      {
        if (is_parser_in_tag_[i])
        {
          switch(i)
          {
          case FEATINTENSITY:   feature_->setIntensity(asDouble_(sm_.convert(chars))); break;
          case POSITION:        feature_->getPosition()[current_pcoord_] = asDouble_(sm_.convert(chars)); break;
          case QUALITY:         feature_->setQuality(current_qcoord_,asDouble_(sm_.convert(chars))); break;
          case OVERALLQUALITY:  feature_->setOverallQuality(asDouble_(sm_.convert(chars))); break;
          case CHARGE:          feature_->setCharge(asInt_(sm_.convert(chars))); break;
          case HPOSITION:       (*hull_position_)[current_hcoord_] = asDouble_(sm_.convert(chars)); break;
          case PAIRQUALITY:     pair_->setQuality(asDouble_(sm_.convert(chars)));
          }
        }
      }
    }

    void FeaturePairsHandler::writeTo(std::ostream& os)
    {

      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
      os << "<featurePairs xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeaturePairsXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << std::endl;

      // write features with their attributes
      for (UInt s=0; s<cpairs_->size(); s++)
      {
        const ElementPair< Feature >& pair = (*cpairs_)[s];

        os << "<pair nr=\"" << s << "\">" << std::endl;
        os << "\t<pairquality>" << pair.getQuality() << "</pairquality>" << std::endl;

        os << "\t<first>" << std::endl;
        Feature first = pair.getFirst();
        writeFeature_(os,first);
        os << "\t</first>" << std::endl;

        os << "\t<second>" << std::endl;
        Feature seco  = pair.getSecond();
        writeFeature_(os,seco);
        os << "\t</second>" << std::endl;

        os << "</pair>" << std::endl;

      } // end for ( features )

      os << "</featurePairs>" << std::endl;
      os <<
      "<!-- Local Variables: -->\n"
      "<!-- mode: nxml -->\n"
      "<!-- tab-width: 2 -->\n"
      "<!-- End: -->\n"
      ;
    }

    void FeaturePairsHandler::writeFeature_(std::ostream& os, Feature dfeat)
    {
      os << "\t<feature id=\"" << id_generator_.getUID() << "\">" << std::endl;

      Feature::PositionType pos = dfeat.getPosition();
      UInt dpos_size = pos.size();

      for (UInt i=0; i<dpos_size;i++)
      {
        os << "\t\t<position dim=\"" << i << "\">" << pos[i] << "</position>" <<  std::endl;
      }

      os << "\t\t<intensity>" << dfeat.getIntensity() << "</intensity>" << std::endl;

      for (UInt i=0; i<dpos_size;i++)
      {
        os << "\t\t<quality dim=\"" << i << "\">" << dfeat.getQuality(i) << "</quality>" << std:: endl;
      }

      os << "\t\t<overallquality>" << dfeat.getOverallQuality() << "</overallquality>" << std:: endl;
      os << "\t\t<charge>" << dfeat.getCharge() << "</charge>" << std:: endl;

      // write model description
      ModelDescription<2> desc = dfeat.getModelDescription();
      os << "\t\t<model name=\"" << desc.getName() << "\">" << std:: endl;
      Param modelp = desc.getParam();
      Param::ParamIterator piter = modelp.begin();
      while (piter != modelp.end())
      {
        os << "\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\">";
        os << "</param>" << std::endl;
        piter++;
      }
      os << "\t\t</model>" << std::endl;

      // write convex hull
      Feature::ConvexHullVector hulls = dfeat.getConvexHulls();
      Feature::ConvexHullVector::iterator citer = hulls.begin();

      UInt hulls_count = hulls.size();

      for (UInt i=0;i<hulls_count; i++)
      {
        os << "\t\t<convexhull nr=\"" << i << "\">" << std:: endl;

        ConvexHull2D current_hull = hulls[i];
        UInt hull_size = current_hull.getPoints().size();

        for (UInt j=0;j<hull_size;j++)
        {
          os << "\t\t\t<hullpoint>" << std::endl;

          Feature::PositionType pos = current_hull.getPoints()[j];
          UInt pos_size = pos.size();
          for (UInt k=0; k<pos_size; k++)
          {
            os << "\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << std::endl;
          }

          os << "\t\t\t</hullpoint>" << std::endl;
        } // end for (..hull_size..)

        os << "\t\t</convexhull>" << std::endl;
      } // end  for ( ... hull_count..)

      os << "\t</feature>\n";
    }
  } // namespace Internal

} // namespace OpenMS






