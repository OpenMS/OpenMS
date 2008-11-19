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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

#ifdef ITRAQ_NAIVECORRECTION
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#endif

namespace OpenMS
{
	ItraqQuantifier::ItraqQuantifier()
		:DefaultParamHandler("ItraqQuantifier"),
		 itraq_type_(FOURPLEX),
		 isotope_corrections_()
	{
		initIsotopeCorrections_();
		setDefaultParams_();
	}
	
	ItraqQuantifier::ItraqQuantifier(Int itraq_type)
		:DefaultParamHandler("ItraqQuantifier"),
		 itraq_type_(itraq_type),
		 isotope_corrections_()
	{
		initIsotopeCorrections_();
		setDefaultParams_();
	}

	ItraqQuantifier::ItraqQuantifier(Int itraq_type, const Param& param)
		:DefaultParamHandler("ItraqQuantifier"),
		 itraq_type_(itraq_type),
		 isotope_corrections_()
	{
		initIsotopeCorrections_();
		setDefaultParams_();
		setParameters(param);
		updateMembers_();
	}

	ItraqQuantifier::ItraqQuantifier(const ItraqQuantifier& cp)
	: DefaultParamHandler(cp),
		ItraqConstants(cp),
		itraq_type_(cp.itraq_type_), 
    channel_map_(cp.channel_map_),
		isotope_corrections_(cp.isotope_corrections_)
	{
	}

	ItraqQuantifier& ItraqQuantifier::operator = (const ItraqQuantifier& rhs)
	{
		if (this == &rhs) return *this;
		
		DefaultParamHandler::operator = (rhs);
		ItraqConstants::operator = (rhs);
		itraq_type_ = rhs.itraq_type_; 
    channel_map_ = rhs.channel_map_;
		isotope_corrections_ = rhs.isotope_corrections_;
		
		return *this;
	}

	/**
	 *	@brief using the raw iTRAQ intensities we apply isotope correction, normalization (using median) and protein inference
	 *	
	 *	@param consensus_map_in Raw iTRAQ intensities from previous step
	 *	@param consensus_map_out Postprocessed iTRAQ ratios for peptides
	 *
	 *	@throws Exception::FailedAPICall is least-squares fit fails
	 *	@throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
	 */
	void ItraqQuantifier::run(const ConsensusMap& consensus_map_in, 
					 ConsensusMap& consensus_map_out
					 )
	{
		run(consensus_map_in, 
				std::vector< PeptideIdentification >(),
				std::vector< ProteinIdentification >(),
				consensus_map_out);
		return;
	}
	/**
	 *	@brief using the raw iTRAQ intensities we apply isotope correction, normalization (using median) and protein inference
	 *	
	 *	@param consensus_map_in Raw iTRAQ intensities from previous step
	 *	@param peptide_ids List of peptides identified by a search engine on the same MS² dataset
	 *	@param protein_ids List of proteins inferred from peptides
	 *	@param consensus_map_out Postprocessed iTRAQ ratios for Proteins (if provided) or Peptides otherwise
	 *
	 *	@throws Exception::FailedAPICall is least-squares fit fails
	 *	@throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
	 */
	void ItraqQuantifier::run(const ConsensusMap& consensus_map_in, 
					 const std::vector< PeptideIdentification > &peptide_ids,
					 const std::vector< ProteinIdentification > &protein_ids, 
					 ConsensusMap& consensus_map_out
					 )
	{

		reconstruct_channel_info_(consensus_map_in);

		consensus_map_out = consensus_map_in;

		// first do isotope correction
		if (String(param_.getValue("isotope_correction")) == "true")
		{
			// translate isotope_corrections_ to a channel_frequency matrix
			Matrix<double> channel_frequency(CHANNEL_COUNT[itraq_type_], CHANNEL_COUNT[itraq_type_]);
			for (Int i=0; i < CHANNEL_COUNT[itraq_type_]; ++i)
			{
				for (Int j=0; j < CHANNEL_COUNT[itraq_type_]; ++j)
				{
					// diagonal (should be close to 1 = 100%)
					if (i==j)
					{
						double val = 1.0;

						// subtract all isotope deviations of row i
						for (Int col_idx=0; col_idx < 4; ++col_idx) 
						{
							val += -isotope_corrections_[itraq_type_].getValue(i,col_idx) / 100;
						}
						channel_frequency.setValue(i,j,val);
					}
					else
					{ // from mass i to mass j (directly copy the deviation)
						if (i-j<=2 && i-j>0)
						{
							channel_frequency.setValue(j,i, isotope_corrections_[itraq_type_].getValue(i,j-i+2) / 100);
						}
						else if (j-i<=2 && j-i>0)
						{
							channel_frequency.setValue(j,i, isotope_corrections_[itraq_type_].getValue(i,j-i+1) / 100);
						}
					}
				}
			}

			#ifdef ITRAQ_DEBUG
			std::cout << "channel_frequency matrix: \n" << channel_frequency << "\n" << std::endl;
			#endif

			#ifdef ITRAQ_NAIVECORRECTION
			// this solves the system naively
			int gsl_status = 0;
			gsl_matrix* gsl_m = channel_frequency.toGslMatrix();
			gsl_permutation* gsl_p = gsl_permutation_alloc (channel_frequency.rows());
			int* gsl_sign = new int(0);
			gsl_vector* gsl_b = gsl_vector_alloc (CHANNEL_COUNT[itraq_type_]);
			gsl_vector* gsl_x = gsl_vector_alloc (CHANNEL_COUNT[itraq_type_]);
			// lets see if the matrix is invertible
			gsl_status = gsl_linalg_LU_decomp (gsl_m, gsl_p, gsl_sign);
			if (gsl_status!=0)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; the Matrix is not invertible!");
			}
			#endif
			Matrix<double> m_b(CHANNEL_COUNT[itraq_type_], 1);
			Matrix<double> m_x(CHANNEL_COUNT[itraq_type_], 1);
			
			// correct all consensus elements
			for (size_t i=0; i< consensus_map_out.size(); ++i)
			{
				#ifdef ITRAQ_DEBUG
				std::cout << "\nMAP element  #### " << i << " #### \n" << std::endl;
				#endif

				consensus_map_out[i].clear(); // delete only the consensus handles
				// fill b vector
				for (ConsensusFeature::HandleSetType::const_iterator it_elements = consensus_map_in[i].getFeatures().begin();
						 it_elements != consensus_map_in[i].getFeatures().end();
						 ++it_elements)
				{

					
					//find channel_id of current element
					Int index = Int(consensus_map_in.getFileDescriptions() [it_elements->getMapIndex()].getMetaValue("channel_id"));
					#ifdef ITRAQ_DEBUG
					std::cout << "	map_index " << it_elements->getMapIndex() << "-> id " << index << " with intensity " << it_elements->getIntensity() <<"\n" << std::endl;
					#endif
					
					#ifdef ITRAQ_NAIVECORRECTION
					// this is deprecated, but serves as quality measurement
					gsl_vector_set (gsl_b, index, it_elements->getIntensity());
					#endif
					m_b(index,0) = it_elements->getIntensity();
				}

				// solve
				#ifdef ITRAQ_NAIVECORRECTION
				gsl_status = gsl_linalg_LU_solve (gsl_m, gsl_p, gsl_b, gsl_x);
				if (gsl_status!=0)
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; the Matrix is not invertible!");
				}
				#endif
				Int status = NonNegativeLeastSquaresSolver::solve(channel_frequency,m_b,m_x);
				if (status!=NonNegativeLeastSquaresSolver::SOLVED)
				{
					throw Exception::FailedAPICall(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Failed to find least-square fit!");
				}

				// write back the values to the map
				for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_in[i].begin();
						 it_elements != consensus_map_in[i].end();
						 ++it_elements)
				{
					FeatureHandle handle = *it_elements;
					//find channel_id of current element
					Int index = Int(consensus_map_out.getFileDescriptions() [it_elements->getMapIndex()].getMetaValue("channel_id"));

					#ifdef ITRAQ_NAIVECORRECTION
					//this has become useless (even for comparison of methods)
					//handle.setIntensity ( Peak2D::IntensityType(gsl_vector_get (gsl_x, index)) );
					#endif
					handle.setIntensity ( Peak2D::IntensityType( m_x(index, 0)) );
					
					consensus_map_out[i].insert(handle);
					
					#ifdef ITRAQ_DEBUG
					std::cout <<  it_elements->getIntensity() << " -> " << handle.getIntensity () << std::endl;
					#endif
				}

			}
			
			#ifdef ITRAQ_NAIVECORRECTION
			// clean up
			gsl_matrix_free(gsl_m);
			gsl_permutation_free(gsl_p);
			delete gsl_sign;
			gsl_vector_free(gsl_b);
			gsl_vector_free(gsl_x);
			#endif
		} // ! isotope_correction


		// ** NORMALIZATION ** //

		// normalize median of channel-to-reference ratio to 1
		Int reference_channel = Int(param_.getValue("channel_reference"));
		#ifdef ITRAQ_DEBUG
		std::cout << "reference_channel is: " << reference_channel  << std::endl;
		#endif

		Map <Int, Int> map_to_vectorindex;
		Int ref_mapid = 0;
		Int index = 0;
		for (ConsensusMap::FileDescriptions::const_iterator file_it = consensus_map_out.getFileDescriptions().begin(); 
			 file_it!=consensus_map_out.getFileDescriptions().end(); 
			 ++file_it)
		{
			if ((Int) file_it->second.getMetaValue ("channel_name") == reference_channel) 
			{
				ref_mapid = file_it->first;
				#ifdef ITRAQ_DEBUG
				std::cout << "reference_map_id is: " << ref_mapid <<  std::endl;
				#endif	
			}
			map_to_vectorindex[file_it->first] = index;
			index++;
		}

		if (channel_map_.has(reference_channel))
		{
			std::vector< std::vector<double> > peptide_ratios;
			// this is a control (the normalization factors should be about the same)
			std::vector< std::vector<double> > peptide_intensities;

			// build mapping of map_index to ratio_array_index
			peptide_ratios.resize(channel_map_.size());
			peptide_intensities.resize(channel_map_.size());

			//build up ratios for each peptide of non-reference channels
			ConsensusFeature::HandleSetType::iterator ref_it;
			Peak2D::IntensityType ref_intensity;
			for (size_t i=0; i< consensus_map_out.size(); ++i)
			{
				// find reference index (this is inefficient to do every time, 
				// but the most robust against anyone who tries to change the internals of ConsensusFeature):
				ref_it=consensus_map_out[i].end();
				for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
						 it_elements != consensus_map_out[i].end();
						 ++it_elements)
				{
					if ((Int) consensus_map_out.getFileDescriptions() [it_elements->getMapIndex()].getMetaValue("channel_name") == reference_channel)
					{
						ref_it = it_elements;
						break;
					}
				}

				// reference channel not found in this ConsensusFeature
				if (ref_it==consensus_map_out[i].end())
				{
					std::cerr << "ItraqQuantifier::run() WARNING: ConsensusFeature " << i << " does not have a reference channel! Skipping" << std::endl;
					continue;
				}

				ref_intensity = ref_it->getIntensity();

				// now collect the ratios and intensities
				for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
						 it_elements != consensus_map_out[i].end();
						 ++it_elements)
				{
					peptide_ratios[map_to_vectorindex[it_elements->getMapIndex()]].push_back(it_elements->getIntensity() / ref_intensity);
					// control
					peptide_intensities[map_to_vectorindex[it_elements->getMapIndex()]].push_back(it_elements->getIntensity());
				}					
			} // ! collect ratios

			double max_deviation_from_control = 0;
			// find MEDIAN of ratios for each channel (store as 0th element in sorted vector)
			for (Map <Int, Int>::const_iterator it_map = map_to_vectorindex.begin(); it_map != map_to_vectorindex.end(); ++it_map)
			{
				// sort vector (partial_sort might improve performance here)
				std::sort(peptide_ratios[it_map->second].begin(), peptide_ratios[it_map->second].end());
				// save median as first element 
				peptide_ratios[it_map->second][0] = peptide_ratios[it_map->second][peptide_ratios[it_map->second].size()/2];

				// sort control (intensities)
				std::sort(peptide_intensities[it_map->second].begin(), peptide_intensities[it_map->second].end());
				// find MEDIAN of control-method (intentities) for each channel
				peptide_intensities[it_map->second][0] = peptide_intensities[it_map->second][peptide_intensities[it_map->second].size()/2] /
																								 peptide_intensities[ref_mapid][peptide_intensities[ref_mapid].size()/2];
				#ifdef ITRAQ_DEBUG
					std::cout << "iTRAQ-normalize:  map-id " << (it_map->first) << " has factor " << (peptide_ratios[it_map->second][0]) << std::endl;
					std::cout << "iTRAQ-normalize:  map-id " << (it_map->first) << " has factor " << (peptide_intensities[it_map->second][0]) << " (control)" << std::endl;
				#endif
				double dev = (peptide_ratios[it_map->second][0]-peptide_intensities[it_map->second][0])/peptide_ratios[it_map->second][0];
				if	(fabs(max_deviation_from_control) < fabs(dev))
				{
					max_deviation_from_control = dev;
				}
			}

			std::cout << "iTRAQ-normalization: max ratio deviation of alternative method is " << (max_deviation_from_control*100) << "%\n";

			#ifdef ITRAQ_DEBUG
			std::cout << "debug OUTPUT\n";
			for (UInt i = 1; i<peptide_ratios[0].size(); ++i)
			{
				if (i == peptide_intensities[0].size()/2)
				{
					std::cout << "++++++++++ median: \n";
				}
				for (UInt j = 0; j<peptide_ratios.size(); ++j)
				{
					std::cout << peptide_ratios[j][i] << " ";
				}
				std::cout << " -- int -- ";
				for (UInt j = 0; j<peptide_intensities.size(); ++j)
				{
					std::cout << peptide_intensities[j][i] << " ";
				}
				if (i == peptide_intensities[0].size()/2)
				{
					std::cout << "\n----------- median: ";
				}
				std::cout << "\n";
			}
			#endif

			// adjust intensity ratios 
			for (size_t i=0; i< consensus_map_out.size(); ++i)
			{
				// find reference index (this is inefficient to do every time, 
				// but the most robust against anyone who tries to change the internals of ConsensusFeature):
				ref_it=consensus_map_out[i].end();
				for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
						 it_elements != consensus_map_out[i].end();
						 ++it_elements)
				{
					if ((Int) consensus_map_out.getFileDescriptions() [it_elements->getMapIndex()].getMetaValue("channel_name") == reference_channel)
					{
						ref_it = it_elements;
						break;
					}
				}

				// reference channel not found in this ConsensusFeature
				if (ref_it==consensus_map_out[i].end())
				{
					continue;
				}

				ref_intensity = ref_it->getIntensity();

				// now adjust the ratios
				ConsensusFeature cf = consensus_map_out[i];
				cf.clear(); // delete its handles
				for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map_out[i].begin();
						 it_elements != consensus_map_out[i].end();
						 ++it_elements)
				{
					FeatureHandle hd = *it_elements;
					if (it_elements == ref_it) 
					{
						hd.setIntensity(1);
					}
					else
					{ // divide current intensity by normalization factor (which was stored at position 0)
						hd.setIntensity(hd.getIntensity() / peptide_ratios[map_to_vectorindex[it_elements->getMapIndex()]][0]);
					}
					cf.insert(hd);
				}					
				// replace consensusFeature with updated intensity
				consensus_map_out[i] = cf;
			} // ! adjust ratios

		} // ! ref_channel valid
		else
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier::run() Parameter 'channel_reference' does not name a valid channel!");
		}



		// ** PEPTIDE PROTEIN MAPPING ** //

		// find unique peptides of a protein and use a robust peptide ratio as protein ratio
		if ((!peptide_ids.empty() && !protein_ids.empty()))
		{
			
			// annotate consensusMap with identifications
			IDMapper mapper;
			mapper.setRTDelta(0.005);
			mapper.setMZDelta(0.0005);
			mapper.annotate(consensus_map_out, peptide_ids, protein_ids, false);
			
			// put quantitative info on Proteins
			ProteinInference inferrer;
			inferrer.infer(consensus_map_out, ref_mapid);

		}

		consensus_map_out.setExperimentType("itraq");
		
		return;
	}

	void ItraqQuantifier::setDefaultParams_()
	{
		defaults_.setValue("isotope_correction", "true", "enable isotope correction (highly recommended)", StringList::create("advanced")); 
		defaults_.setValidStrings("isotope_correction", StringList::create("true,false"));

		defaults_.setValue("isotope_correction_values", "", "override default values (see Documentation); use the following format: <channel>:<v1>/<v2>/<v3>/<v4> ; e.g. '114:0/0.3/4/0 , 116:0.1/0.3/3/0.2' ", StringList::create("advanced"));

		if (itraq_type_ == ItraqConstants::FOURPLEX)
		{
			defaults_.setValue("channel_reference", 114, "number of the reference channel"); 
			defaults_.setMinInt("channel_reference",114);
			defaults_.setMaxInt("channel_reference",117);
		}
		else
		{
			defaults_.setValue("channel_reference", 113, "number of the reference channel"); 
			defaults_.setMinInt("channel_reference",113);
			defaults_.setMaxInt("channel_reference",121);			
		}			

		defaultsToParam_();
	}


	void ItraqQuantifier::updateMembers_()
	{

		// update isotope_corrections_ Matrix with custom values
		StringList channels = StringList::create(param_.getValue("isotope_correction_values"));
		if (channels.size()>0)
		{
			// split the channels key:name pairs apart
			for (StringList::const_iterator it=channels.begin();it!=channels.end();++it)
			{
				StringList result;
				it->split(':',result);
				if (result.size()!=2)
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected one ':', got this: '" + (*it) + "'");
				}
				result[0] = result[0].trim(); // hold channel name
				result[1] = result[1].trim(); // holds 4 values

				Int channel = result[0].toInt();
				Int line = (itraq_type_ == FOURPLEX ? channel-114 : channel-113);
				if ((itraq_type_ == FOURPLEX && (line<0 || line>3))
						|| 
						(itraq_type_ == EIGHTPLEX && (line<0 || line>8) || channel==120))
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,String("ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; channel-name is not valid for ") + String(itraq_type_==FOURPLEX ? "4plex": "8plex") + String(": '") + result[0] + String("'"));
				}
				// if set to 121 we still want to change line 7 of the matrix
				if (line==8) line=7;

				StringList corrections;
				result[1].split('/',corrections);
				if (corrections.size()!=4)
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected four correction values separated by '&', got this: '" + result[1] + "'");
				}

				// overwrite line in Matrix with custom values
				isotope_corrections_[itraq_type_].setValue(line,0, corrections[0].toDouble());
				isotope_corrections_[itraq_type_].setValue(line,1, corrections[1].toDouble());
				isotope_corrections_[itraq_type_].setValue(line,2, corrections[2].toDouble());
				isotope_corrections_[itraq_type_].setValue(line,3, corrections[3].toDouble());

				#ifdef ITRAQ_DEBUG
				std::cout << "Channel " << channel << " has values " << corrections << std::endl;
				#endif


			}
		}
	} // ! update_members_

	/// initialize
	void ItraqQuantifier::initIsotopeCorrections_() 
	{
		isotope_corrections_.resize(2);
		isotope_corrections_[0].setMatrix<4,4>(ISOTOPECORRECTIONS_FOURPLEX);
		isotope_corrections_[1].setMatrix<8,4>(ISOTOPECORRECTIONS_EIGHTPLEX);
	}

	/// extract channel information (active channels, names, etc) from ConsensusMap
	void ItraqQuantifier::reconstruct_channel_info_(const ConsensusMap& consensus_map)
	{
		channel_map_.clear();

		for (ConsensusMap::FileDescriptions::const_iterator file_it = consensus_map.getFileDescriptions().begin(); 
				 file_it!=consensus_map.getFileDescriptions().end(); 
				 ++file_it)
		{
			if (file_it->second.metaValueExists("channel_name"))
			{
				ChannelInfo info;
				// fill info
				info.name = file_it->second.getMetaValue ("channel_name");
				info.id = file_it->second.getMetaValue ("channel_id");
				info.description = file_it->second.getMetaValue ("channel_description");
				info.center = file_it->second.getMetaValue ("channel_center");
				info.active = (String(file_it->second.getMetaValue ("channel_active")) == "true" ? true : false);
				channel_map_[info.name] = info;
				#ifdef ITRAQ_DEBUG
				std::cout << " setting info.name " << (info.name) << " and id " << (info.id) << std::endl;
				#endif
			}
			else
			{
				throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "ItraqQuantifier::reconstruct_channel_info_ The ConsensusMap provided is missing MetaInfo from ItraqChannelExtractor!");
			}
		}
	}
		
	// currently from http://www.matrixscience.com/help/quant_config_help.html
	const double ItraqQuantifier::ISOTOPECORRECTIONS_FOURPLEX[4][4] = {
		{0.0, 1.0, 5.9, 0.2},		//114
		{0.0, 2.0, 5.6, 0.1},
		{0.0, 3.0, 4.5, 0.1},
		{0.1, 4.0, 3.5, 0.0}		//117
	};
	
	//taken from Applied Biosystems Website
	// http://faqs.appliedbiosystems.com/cgi-bin/appliedbio.cfg/php/enduser/std_adp.php?p_faqid=3671
	const double ItraqQuantifier::ISOTOPECORRECTIONS_EIGHTPLEX[8][4] = {
		{0.00, 2.50, 6.89, 0.22},		//113
		{0.00, 0.94, 5.90, 0.16},
		{0.00, 1.88, 4.90, 0.10},
		{0.00, 2.82, 3.90, 0.07},
		{0.06, 3.77, 2.88, 0.00},
		{0.09, 4.71, 1.88, 0.00},
		{0.14, 5.66, 0.87, 0.00},
		{0.27, 7.44, 0.18, 0.00}		//121
	};
 
	// number of channels for iTRAQ types. (make sure it corresponds to enum ITRAQ_TYPES)
	const Int ItraqQuantifier::CHANNEL_COUNT[2] = {4, 8};
}
 
