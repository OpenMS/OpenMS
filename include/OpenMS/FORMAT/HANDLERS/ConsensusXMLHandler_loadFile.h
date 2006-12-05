// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>

namespace OpenMS
{


  template < typename  AlignmentT >
  void ConsensusXMLHandler<AlignmentT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    int tag = enterTag(qname, attributes);

    String tmp_str;
    // Do something depending on the tag
    switch(tag)
    {
        case MAPLIST:
        tmp_str = getAttributeAsString(COUNT);
        if (tmp_str != "")
        {
          UnsignedInt count = asUnsignedInt_(tmp_str);
          consensus_map_->getMapVector().resize(count);
          consensus_map_->getFilenames().resize(count);
        }
        break;
        case MAPTYPE:
        tmp_str = getAttributeAsString(NAME);
        if (tmp_str != "")
        {
          if (getAttributeAsString(ID) == "feature_map")
          {
            feature_map_flag_ = true;
            consensus_map_flag_ = false;
          }
          else
            if (getAttributeAsString(ID) == "consensus_map")
            {
              consensus_map_flag_ = true;
              feature_map_flag_ = false;
            }

        }
        break;
        case MAP:
        tmp_str = getAttributeAsString(ID);
        if (tmp_str != "")
        {
          UnsignedInt id = asUnsignedInt_(tmp_str);
          tmp_str = getAttributeAsString(NAME);
          if (tmp_str != "")
          {
            String act_filename = tmp_str;
            // load FeatureMapXML
            if (feature_map_flag_)
            {
              DFeatureMapFile feature_file;
              loadFile_(feature_file,act_filename,id);
            }
            // load MzData
            else
              if (consensus_map_flag_)
              {
                ConsensusXMLFile cons_file;
                loadFile_(cons_file,act_filename,id);
              }
              else
              {
                MzDataFile mzdata_file;
                loadFile_(mzdata_file,act_filename,id);
              }

            consensus_map_->getFilenames()[id]=act_filename;
          }
        }
        break;
        case CONSENSUSELEMENT:
        act_cons_element_ = ConsensusElementType();
        consensus_element_range_ = true;
        break;
        case CENTROID:
        tmp_str = getAttributeAsString(RT_ATT);
        if (tmp_str != "")
        {
          pos_[RT] = asDouble_(tmp_str);
        }

        tmp_str = getAttributeAsString(MZ_ATT);
        if (tmp_str != "")
        {
          pos_[MZ] = asDouble_(tmp_str);
        }

        tmp_str = getAttributeAsString(IT);
        if (tmp_str != "")
        {
          it_ = asDouble_(tmp_str);
        }
        break;
        case RANGE:
        if (consensus_element_range_)
        {
          tmp_str = getAttributeAsString(RTMIN);
          if (tmp_str != "")
          {
            pos_range_.setMinX(asDouble_(tmp_str));

            tmp_str = getAttributeAsString(RTMAX);
            if (tmp_str != "")
            {
              pos_range_.setMaxX(asDouble_(tmp_str));

              tmp_str = getAttributeAsString(MZMIN);
              if (tmp_str != "")
              {
                pos_range_.setMinY(asDouble_(tmp_str));

                tmp_str = getAttributeAsString(MZMAX);
                if (tmp_str != "")
                {
                  pos_range_.setMaxY(asDouble_(tmp_str));

                  tmp_str = getAttributeAsString(ITMIN);
                  if (tmp_str != "")
                  {
                    it_range_.setMin(asDouble_(tmp_str));

                    tmp_str = getAttributeAsString(ITMAX);
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
        break;
        case ELEMENT:
        IndexTuple< ElementContainerType > act_index_tuple;
        tmp_str = getAttributeAsString(MAP_ATT);
        if (tmp_str != "")
        {
          UnsignedInt map_index = asUnsignedInt_(tmp_str);
          tmp_str = getAttributeAsString(ID);
          if (tmp_str != "")
          {
            UnsignedInt element_index = asUnsignedInt_(tmp_str);

            act_index_tuple.getMapIndex() = map_index;
            act_index_tuple.getElementIndex() = element_index;
            act_index_tuple.setElement(((consensus_map_->getMapVector())[map_index])[element_index]);
            act_cons_element_.insert(act_index_tuple);
            act_cons_element_.getPosition() = pos_;
            act_cons_element_.getPositionRange() = pos_range_;
            act_cons_element_.getIntensity() = it_;
            act_cons_element_.getIntensityRange() = it_range_;
          }
        }
        break;
    }
  }

  /// Load the peaks
  template <>
  template <>
  void ConsensusXMLHandler< StarAlignment< ConsensusFeature<> > >::loadFile_< DFeatureMapFile >(DFeatureMapFile& feature_file, const String& file_name, UnsignedInt id)
  {
    feature_file.load(file_name,(consensus_map_->getMapVector())[id]);
  }

  // load MzData
  template <>
  template <>
  void ConsensusXMLHandler< StarAlignment< ConsensusPeak<> > >::loadFile_<MzDataFile>(MzDataFile& mzdata_file, const String& file_name, UnsignedInt id)
  {
    MSExperiment< Peak > ms_exp;
    mzdata_file.load(file_name,ms_exp);
    ms_exp.get2DData((consensus_map_->getMapVector())[id]);
  }

  // load consensusXML
  template <>
  template <>
  void ConsensusXMLHandler< StarAlignment< ConsensusFeature< ConsensusMap< ConsensusFeature<> > > > >::loadFile_<ConsensusXMLFile>(ConsensusXMLFile& cons_file, const String& file_name, UnsignedInt id)
  {
    cons_file.load(file_name,(consensus_map_->getMapVector())[id]);
  }
} // namespace OpenMS

#endif // OPENMS_FOMAT_MZXMLFILE_H
