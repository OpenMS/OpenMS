/// -*- Mode: C++; tab-width: 2; -*-
/// vi: set ts=2:
//
/// --------------------------------------------------------------------------
///                   OpenMS Mass Spectrometry Framework
/// --------------------------------------------------------------------------
///  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
///  This library is free software; you can redistribute it and/or
///  modify it under the terms of the GNU Lesser General Public
///  License as published by the Free Software Foundation; either
///  version 2.1 of the License, or (at your option) any later version.
//
///  This library is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
///  Lesser General Public License for more details.
//
///  You should have received a copy of the GNU Lesser General Public
///  License along with this library; if not, write to the Free Software
///  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
/// --------------------------------------------------------------------------
/// $Id: Outfile.C,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Outfile.h>

namespace OpenMS 
{
	Outfile::Outfile()
	: queries_(),
		peptide_hits_(),
		protein_hits_(),
		precursor_retention_times_(),
		precursor_mz_values_(),
		ok_(false)
	{}
	
	Outfile::Outfile(const Outfile& source)
		: queries_(source.queries_),
			protein_ids_(source.protein_ids_),
    	peptide_hits_(source.peptide_hits_),
    	protein_hits_(source.protein_hits_),
    	precursor_retention_times_(source.precursor_retention_times_),
    	precursor_mz_values_(source.precursor_mz_values_),
    	ok_(source.ok_),
			curr_query_(source.curr_query_),
			curr_peptide_hit_(source.curr_peptide_hit_),
			curr_protein_hit_(source.curr_protein_hit_)
  {
  }

  Outfile::~Outfile()
  {
		queries_.clear();
		protein_ids_.clear();
		peptide_hits_.clear();
		protein_hits_.clear();
		precursor_retention_times_.clear();
		precursor_mz_values_.clear();
  	peptide_hits_.clear();
  	protein_hits_.clear();
  }

  bool Outfile::ok() const
  {
  	return ok_;
  }
  
  Outfile& Outfile::operator>>(Identification& query)
  {
  	query.clear();
		if ( curr_query_ == queries_.end() ) return *this;
		
		query = Identification(*curr_query_);
		++curr_query_;
    
    return *this;
  }
  

  Outfile& Outfile::operator>>(PeptideHit& peptide_hit)
  {
    peptide_hit.clear();
		
		if (curr_peptide_hit_ == peptide_hits_.end())
		{
      return *this;			
		}
		
		/// copy values from current hit
		peptide_hit = PeptideHit(*curr_peptide_hit_);
    ++curr_peptide_hit_;
    
    return *this;
  }

  Outfile& Outfile::operator>>(ProteinHit& protein_hit)
  {
    protein_hit.clear();
		
		if (curr_protein_hit_ == protein_hits_.end())
		{
      return *this;			
		}
		/// copy values from current hit
		protein_hit = ProteinHit(*curr_protein_hit_);
    ++curr_protein_hit_;
    
    return *this;
  }

	/// Assignment operator
  Outfile& Outfile::operator=(const Outfile& source)
  {
  	if (this == &source)
  	{
  		return *this;
  	}
  	precursor_retention_times_ = source.precursor_retention_times_;  			
  	precursor_mz_values_ = source.precursor_mz_values_;  			
  	queries_ = source.queries_;
  	peptide_hits_ = source.peptide_hits_;
  	protein_hits_ = source.protein_hits_;

		return *this;
  }
  
  /// returns the retention time of the  search
  const std::vector<float>& Outfile::getPrecursorRetentionTimes() const
  {
  	return precursor_retention_times_;
  }

  /// sets the retention time of the  search
  void Outfile::setPrecursorRetentionTimes(const std::vector<float>& precursor_retention_times)
	{
		precursor_retention_times_ = precursor_retention_times;
	}
	
  /// returns the m/z of the precursor peak of the  search
  const std::vector<float>& Outfile::getPrecursorMZValues() const
  {
  	return precursor_mz_values_;
  }

  /// sets the m/z of the precursor peak of the  search
  void Outfile::setPrecursorMZValues(const std::vector<float>& precursor_mz_values)
  {
  	precursor_mz_values_ = precursor_mz_values;
  }
  
  const std::vector< Identification >& Outfile::getIdentifications() const
	{
		return queries_; 
	}
	
  void Outfile::setIdentifications(const std::vector< Identification >& queries)
  {
  	queries_ = queries;
  }
  
  const ProteinIdentification& Outfile::getProteinIdentification() const
	{
		return protein_ids_; 
	}
	
  void Outfile::setProteinIdentification(const ProteinIdentification& protein_ids)
  {
  	protein_ids_ = protein_ids;
  }
	
	/// get the accession and accession type from a line
	void Outfile::get_ac_and_ac_type(String line, const std::string& filename, std::string& accession, std::string& accession_type) throw (Exception::ParseError)
	{
		std::pair<std::string, std::string> p;
		/// if it's a FASTA line
		if ( line.hasPrefix(">") )
		{
			line.erase(0,1);
		}
		line.trim();
		
		/// if it's a swissprot accession
		if ( line.hasPrefix("tr") || line.hasPrefix("sp") )
		{
			accession = line.substr(3, line.find('|', 3)-3);
			accession_type = "SwissProt";
		}
		else if ( line.hasPrefix("gi") )
		{
			unsigned int snd = line.find('|', 3);
			unsigned int third = line.find('|', ++snd);
			++third;
			accession = line.substr(third, line.find('|', third)-third);
			accession_type = line.substr(snd, third-1-snd);
			if ( accession_type == "gb" )
			{
				accession_type = "GenBank";
			}
			else if ( accession_type == "emb" )
			{
				accession_type = "EMBL";
			}
			else if ( accession_type == "dbj" )
			{
				accession_type = "DDBJ";
			}
		}
		else if ( line.hasPrefix("ref") )
		{
			accession = line.substr(4, line.find('|', 4) - 4);
			accession_type = "NCBI";
		}
		/// if it's a swissprot line
		else if ( line.hasPrefix("AC") )
		{
			line.erase(0,2);
			accession = line.trim();
			accession_type = "SwissProt";
		}
		else
		{
			accession = line.trim();
			accession_type = "unknown";
			/*std::string error = "unkown accession identifier!";
			error.append(line);
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error.c_str(), filename);*/
		}
		/* GenBank                           gi|gi-number|gb|accession|locus
 EMBL Data Library                 gi|gi-number|emb|accession|locus
 DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
 NCBI Reference Sequence           gi|gi-number|ref|accession|locus
 SWISS-PROT                        gi|gi-number|sp|accession|name
 General database identifier       gnl|database|identifier
 Local Sequence identifier         lcl|identifier*/
	}
	
	/// given a vector of peptide hits, either insert the new peptide hit or update its ProteinHits, returns whether an update took place
	bool Outfile::updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits)
	{
		std::vector< PeptideHit >::iterator pep_hit_i;
		for ( pep_hit_i = peptide_hits.begin(); pep_hit_i != peptide_hits.end(); ++pep_hit_i)
		{
			if ( (pep_hit_i->getSequence() == peptide_hit.getSequence()) && (pep_hit_i->getScore() == peptide_hit.getScore()) ) break;
		}
		
		/// a peptide hit may only be inserted if it's score type matches the one of the existing hits
		if ( (peptide_hits.empty()) || (peptide_hits[0].getScoreType() == peptide_hit.getScoreType()) )
		{
			/// if a new peptide is found, insert it
			if ( pep_hit_i == peptide_hits.end() )
			{
				peptide_hits.push_back(peptide_hit);
				return false;
			}
			/// if the peptide has already been inserted, insert additional protein hits
			else
			{
				std::vector< std::pair< String, String > >::iterator prot_hit_i1, prot_hit_i2;
				/// remove all protein hits from the peptide that are already in the list
				for ( prot_hit_i1 = peptide_hit.getProteinIndices().begin(); prot_hit_i1 != peptide_hit.getProteinIndices().end(); )
				{
					prot_hit_i2 = std::find(pep_hit_i->getProteinIndices().begin(), pep_hit_i->getProteinIndices().end(), *prot_hit_i1);
					if ( prot_hit_i2 != pep_hit_i->getProteinIndices().end() ) peptide_hit.getProteinIndices().erase(prot_hit_i1);
					else ++prot_hit_i1;
				}
				/// add the additional protein hits
				for ( prot_hit_i2 = peptide_hit.getProteinIndices().begin(); prot_hit_i2 != peptide_hit.getProteinIndices().end(); ++prot_hit_i2 )
				{
					pep_hit_i->addProteinIndex(*prot_hit_i2);
				}
				return true;
			}
		}
		return false;
	}

} //namespace OpenMS
