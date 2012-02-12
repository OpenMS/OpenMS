// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticCompressed.h>
#include <stack>
#include <fstream>
#include <cmath>
#include <cstring>
#include <algorithm>

#include <OpenMS/CHEMISTRY/ModifierRep.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

namespace OpenMS
{

/**
@brief comperator for two substings represented as pair of ints

holds a reference of the string and compairs two substrings. It will be used for sorting the indices.
*/
struct SubstringLess : public binary_function<pair<SignedSize,SignedSize> , pair<SignedSize,SignedSize>, bool>
{
	/**
	@brief constructor
	@param str const reference to the string
	*/
	SubstringLess(const String & str) : str_(str) {}
	/**
	@brief copy constructor
	*/
	SubstringLess(const SubstringLess & rhs ) : str_(rhs.str_) {}

	/**
	@brief implementation of the '<' operator for two given substrings
	@param left first substring
	@param right second substring
	@return true if first substring '<' second substring	
	*/
	bool operator()( pair<SignedSize,SignedSize> left, pair<SignedSize,SignedSize> right) const
	{
		return ((str_.substr(left.first,left.second))<(str_.substr(right.first,right.second)));
	}

	protected:
	String const & str_; ///< string
};

/**
@brief comparator for two doubles with a tolerance value
*/
struct FloatsWithTolLess : public binary_function<DoubleReal , DoubleReal, bool>
{
	/**
	@brief constructor
	@param t const reference to the tolerance
	*/
	FloatsWithTolLess(const DoubleReal & t) : tol_(t) {}
	/**
	@brief copy constructor
	*/
	FloatsWithTolLess(const FloatsWithTolLess & rhs ) : tol_(rhs.tol_) {}
	
	/**
	@brief implementation of the '<' operator for two doubles with the tolerance value
	@param f1 first DoubleReal
	@param f2 second DoubleReal
	@return true if first DoubleReal '<' second DoubleReal-tolerance
	*/
	bool operator()( DoubleReal f1, DoubleReal f2) const
	
	{
		return (f1<(f2-tol_));
		//return (fabs(f1 - f2) < tol_);
	}

	protected:
	DoubleReal const & tol_; ///< tolerance value
};


// gets the index of the next separator character
SignedSize SuffixArrayTrypticCompressed::getNextSep_(const SignedSize p) const
{
	for (Size i = (p+1); i < s_.length();++i)
	{
		if (s_[i]=='$')
		{
			return (i);
		}
	}
	return (-1);
}


// getting lowest common prefix of two entrys of suffix array
SignedSize SuffixArrayTrypticCompressed::getLCP_(const pair<SignedSize,SignedSize> & last_point, const pair<SignedSize,SignedSize> & current_point)
{
	SignedSize lastBegin = last_point.first;
	SignedSize currentBegin = current_point.first;
	for (SignedSize i = 0;i<(last_point.second);++i)
	{
		if (i>current_point.second)
		{
			return i;
		}
		if (s_[lastBegin+i]!=s_[currentBegin+i])
		{
			return i;
		}
	}
	return last_point.second;
}


// constructor 
SuffixArrayTrypticCompressed::SuffixArrayTrypticCompressed(const String & st, const String & sa_file_name, const WeightWrapper::WEIGHTMODE weight_mode)
 :WeightWrapper(weight_mode), 
	s_(st),
	tol_(0.5),
	number_of_modifications_(0)
{
	use_tags_=false;
	if (st[0]!='$')
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "String has to start with empyt string ($)","");
	}
	if (st[st.length()-1]!='$')
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "String has to end with separator ($)","");
	}
	//creating array with aminoacid masses
	ResidueDB* rdb = ResidueDB::getInstance();
		
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
	
	for (Size i = 0; i<256;++i)
	{
		masse_[i]=0;
	}

	for (Size i = 0; i<strlen(aa);++i)
	{
		const Residue* r = rdb->getResidue(aa[i]);
		masse_[(int)aa[i]]= getWeight(*r, Residue::Internal);
	}
	
	if (sa_file_name!="")
	{
		open(sa_file_name);
	} 
	else
	{
		//creating unsorted suffix array with every tryptic suffix
		Size next_pos = getNextSep_(0);
		bool is_at_start=true;
		for (Size i = 0 ; i < (s_.length());++i)
		{
			if (next_pos<i) next_pos = getNextSep_(i);
			const char start_char = s_[i];
			// here we have to pay attension that we are not running out of the string, so if we do we set the next character to a character that will allow a digestion before (e.i. for trypsin everything but P)
			const char next_char = (i<s_.length())?s_[i+1]:'T';
			if (start_char=='$')
			{
				is_at_start=true;
			} else 
			{
				// if we have been at a start postion in last step or if we reached a digesting site we add the index the the suffix array
				if (is_at_start||isDigestingEnd(start_char,next_char)) 
				{
					SignedSize start_pos = (is_at_start)?i:(i+1);
					pair<SignedSize,SignedSize> p (start_pos,next_pos-start_pos);
					if (p.second!=0)
					{
						indices_.push_back(p);
					}
					// if we are at start and now at a digesting site we must assure not to forget the index of the distested part
					if (is_at_start&&isDigestingEnd(start_char,next_char))
					{
						
						SignedSize start_pos = (i+1);
						pair<SignedSize,SignedSize> p (start_pos,next_pos-start_pos);
						if (p.second!=0)
						{
							indices_.push_back(p);
						}
					}
					is_at_start=false;
				}
			}
		}
		//sorting array using function using SubstringLess as comparator
		sort(indices_.begin(),indices_.end(),SubstringLess(s_));
		
		//creating lcp
		for (Size i = 1;i<indices_.size();++i)
		{
			lcp_.push_back(getLCP_(indices_.at(i-1),indices_.at(i)));
		}
		lcp_.push_back(0);

		//creating skip entry
		for (Size i = 0; i < lcp_.size();++i)
		{
			SignedSize lcp_value = lcp_.at(i);
			SignedSize j = 0;
			while ((lcp_value>0) && (lcp_.at(j+i)>=lcp_value)&&((Size) j < (lcp_.size()-i))) 
			{
				++j;
			}
			skip_.push_back(j);
		}
		//cout<<"sa_size: "<<indices_.size()<<endl;
	}
	
}

//Copy constructor
SuffixArrayTrypticCompressed::SuffixArrayTrypticCompressed(const SuffixArrayTrypticCompressed & sa)
 :SuffixArray(sa),
	WeightWrapper(sa),
	s_(sa.s_),
	tol_(sa.tol_),
	indices_(sa.indices_),
	lcp_(sa.lcp_),
	skip_(sa.skip_),
	number_of_modifications_(sa.number_of_modifications_) 
{
	use_tags_=(sa.use_tags_);
}

bool SuffixArrayTrypticCompressed::isDigestingEnd(const char aa1, const char aa2) const
	{
		return ((aa1 == 'K' || aa1 == 'R') && aa2 != 'P');
	}

bool SuffixArrayTrypticCompressed::save(const String & file_name)
{
	ofstream file_INDICES;
	ofstream file_LCP;
	ofstream file_SKIP;
	file_INDICES.open ((file_name+".sa2").c_str());
	file_LCP.open ((file_name+".lcp2").c_str());
	file_SKIP.open ((file_name+".skip2").c_str());
	if (!file_INDICES.is_open())
	{
		throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".sa2"));
	}
	if (!file_LCP.is_open())
	{
		throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".lcp2"));
	}
	if (!file_SKIP.is_open())
	{
		throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".skip2"));
	}
	for (Size i = 0; i < indices_.size();++i)
	{
		file_INDICES << indices_.at(i).first<<"\n"<<indices_.at(i).second;
		file_LCP << lcp_.at(i);
		file_SKIP << skip_.at(i);
		
		if (i<(indices_.size()-1))
		{
			file_INDICES<<"\n";
			file_LCP<<"\n";
			file_SKIP<<"\n";
		}
	}
	file_INDICES.close();
	file_LCP.close();
	file_SKIP.close();
	return true;
}

bool SuffixArrayTrypticCompressed::open(const String & file_name)
{
	indices_.clear();
	lcp_.clear();
	skip_.clear();
	ifstream file_INDICES;
	ifstream file_LCP;
	ifstream file_SKIP;
	file_INDICES.open ((file_name+".sa2").c_str());
	file_LCP.open ((file_name+".lcp2").c_str());
	file_SKIP.open ((file_name+".skip2").c_str());
	if (!file_INDICES.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".sa2"));
	}
	if (!file_LCP.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".lcp2"));
	}
	if (!file_SKIP.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, (file_name+".skip2"));
	}
	SignedSize int_INDICES_FIRST;
	while (file_INDICES>>int_INDICES_FIRST)
	{
		SignedSize int_INDICES_SECOND;
		SignedSize int_SKIP;
		SignedSize int_LCP;
		
		file_INDICES >> int_INDICES_SECOND;
		file_LCP >> int_LCP;
		file_SKIP >> int_SKIP;
		pair<SignedSize,SignedSize> p(int_INDICES_FIRST, int_INDICES_SECOND);
		indices_.push_back(p);
		lcp_.push_back(int_LCP);
		skip_.push_back(int_SKIP);
		
	}
	return true;
}

//destructor
SuffixArrayTrypticCompressed::~SuffixArrayTrypticCompressed(){}

//prints the suffix array and lcp and skip table
String SuffixArrayTrypticCompressed::toString()
{
	ostringstream ss;
	vector<pair<SignedSize,SignedSize> >::iterator it = indices_.begin();
	for (;it!=indices_.end();++it)
	{
		SignedSize it_start = (*it).first;
		SignedSize it_l = (*it).second;
		ss<<s_.substr(it_start,it_l);
		ss<<endl;
	}
	ss<<"lcp: ";
	for (Size i = 0;i < lcp_.size();++i)
	{
		ss<<lcp_.at(i);
	}
	ss<<endl<<"skip: ";
	for (Size i = 0;i < skip_.size();++i)
	{
		ss<<skip_.at(i);
	}
	return ss.str();
}

SignedSize SuffixArrayTrypticCompressed::findFirst_ (const vector<DoubleReal> & spec, DoubleReal & m,SignedSize start, SignedSize  end) {
	
	if (end-start<=1) return (spec.at(start)<m-tol_)?end:start;
	SignedSize middle = ((end-start)/2)+start;
	
	if (spec.at(middle)<m-tol_){
		return findFirst_(spec,m,middle,end);
	} 
	if (spec.at(middle)>m+tol_){
		return findFirst_(spec,m,start,middle);
	}
	while (middle>=0&&spec.at(middle)>=m-tol_){
		middle--;
	}
	return (middle+1);
}

SignedSize SuffixArrayTrypticCompressed::findFirst_ (const vector<DoubleReal> & spec, DoubleReal & m) {
	return findFirst_ (spec,m,0,spec.size()-1);
}

// finds all occurences of a given spectrum
void SuffixArrayTrypticCompressed::findSpec(vector<vector<pair<pair<SignedSize,SignedSize>,DoubleReal > > >& candidates, const vector<DoubleReal> & spec )
{
	//time_t t0 (time(NULL));
	if (spec.empty())
	{
		return;
	}

	ModifierRep modifier;
	modifier.setNumberOfModifications(number_of_modifications_);

	Size number_of_posible_mods = modifier.getMaxModificationMasses();

	//check if spectrum is sorted
	for (Size i = 1; i < spec.size();++i)
	{
		if (spec[i - 1] > spec[i])
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Spectrum has to be sorted ascendingly","");
		}
	}
	
	//preparing result table
	for (Size i = 0; i < spec.size();++i)
	{
		vector<pair<pair<SignedSize,SignedSize>,DoubleReal > > v ;
		candidates.push_back(v);
	}
	DoubleReal mmax = spec.back();
	
	//history contains three values: a position within the indices vector, a position for how far we 'walked' into candidate, and the mass at this position, so we initialize with the length of the indices vector, -1 and 0
	stack<pair<pair<SignedSize,map<DoubleReal,SignedSize> >, pair<SignedSize,DoubleReal> > > history ;
	
	
	SignedSize tag_pos = 0;
	
	history.push(pair<pair <SignedSize, map<DoubleReal,SignedSize> >, pair<SignedSize, DoubleReal> >(pair<SignedSize, map<DoubleReal, SignedSize> > (indices_.size() + 1, map<DoubleReal, SignedSize>()), pair<SignedSize, DoubleReal>(-1, getWeight(EmpiricalFormula("H2O")) )));
	
	SignedSize steps = 0;
	SignedSize nres = 0;
	map<DoubleReal, SignedSize> mod_map_start;
	for (SignedSize i = 0; i < (SignedSize)indices_.size(); ++i)
	{
		SignedSize str_len = indices_[i].second;
		// we are looking for the next history entry representing a position we have not been yet
		while (history.top().first.first < i)
		{
			history.pop();
		}
		// mass at this position
		DoubleReal m = history.top().second.second;
		
		map<DoubleReal, SignedSize> modification_map(history.top().first.second);
		
		//if (history.size()==1) 
		//{
		//	modification_map.clear();
		//}
		
		if (tag_pos > history.top().second.first + 1 - 3) 
		{
			tag_pos = -1;
		}
		
		// j indicates how far we have walked into the string
		for (SignedSize j = (history.top().second.first + 1); j < str_len; ++j)
		{
			
			// if we walked out of string something went wrong
			if ((Size)(indices_[i].first + j) >= s_.length()) 
			{
				Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,(indices_[i].first)+j,s_.length());
			}
			
			if (j < 2) tag_pos = -1;
			if (use_tags_ && tag_pos==-1 && j>=2)
			{
				if (binary_search(tags_.begin(),tags_.end(), s_.substr((indices_[i].first)+j-2,3))) 
				{
					tag_pos = j - 2;
				}
			}
			const char c = s_[(indices_[i].first) + j];
			
			if (modification_map.size() < number_of_posible_mods) 
			{
				modifier.refreshModificationList(modification_map, c);
			}
			char cn = ((Size)(indices_[i].first + j + 1) == s_.length() - 1) ? 'R' : s_[(indices_[i].first) + j + 1];
			m += masse_[(int)c];
			
			// there is one special case if we are at a node where the last character before this node is a digesting start and the first outgoing char prevents digestion but not one of the left childs. then we have to pay attention on not skipping this edge that could possible be a peptide candidate
			bool have_to_go_in = false;
			
			if (j + 1 <= lcp_[i] && !isDigestingEnd(c, cn) && isDigestingEnd(c, '$'))
			{
				SignedSize i_copy = i;
				while ((Size)i_copy < indices_.size() && isDigestingEnd(c, s_[(indices_[i_copy].first) + j + 1])) 
				{
					++i_copy;
				}
				if ((Size)i_copy < indices_.size() && i_copy <= i + skip_[i]) 
				{
					have_to_go_in = true;
				}
			}
			
			//cerr << c << " " << cn << " " << isDigestingEnd(c, cn) << " " << (j == (str_len - 1)) << " " << str_len << " " << have_to_go_in << endl;
		
			/*
			try
			{
				cerr << "1.>" << s_.substr(indices_[i].first, j + 1) << " " << getWeight(AASequence(s_.substr(indices_[i].first, j + 1))) << endl;
			}
			catch (...)
			{
			}
			*/

			
			// if we are at a digesting site or at the end of the peptide
			if (isDigestingEnd(c, cn) || j == (str_len - 1) || have_to_go_in)
			{
				if (!use_tags_ || (tag_pos >= 0 && tag_pos <= j - 2)) {
					
					vector<DoubleReal> found_masses;
					if (binary_search(spec.begin(),spec.end(), m, FloatsWithTolLess(tol_)))
					{
						found_masses.push_back(0);
					}
					// if the mass is in spectrum we will add the entry to all matching masses
					map<DoubleReal,SignedSize>::iterator it;
					for (it = modification_map.begin(); it!= modification_map.end();++it)
					{
						if (binary_search(spec.begin(),spec.end(), m + (DoubleReal)it->first, FloatsWithTolLess(tol_)))
						{
							found_masses.push_back(it->first);
						}
					}
				
					for (Size o = 0; o < found_masses.size(); o++) 
					{
						DoubleReal mass_with_mods = (found_masses[o] + m);
						Size first_occ = findFirst_(spec, mass_with_mods);
						Size first_occ_copy = first_occ;
						if (!have_to_go_in)
						{
							++steps;
							++nres;
							pair<pair<SignedSize, SignedSize>, DoubleReal> pnew(pair<SignedSize, SignedSize>(indices_[i].first, j + 1), found_masses[o]);
/*
							try
							{
							cerr << "2.>" << s_.substr(indices_[i].first, j + 1) << " " << getWeight(AASequence(s_.substr(indices_[i].first, j + 1))) << endl;
							}
							catch(...)
							{
							}
	*/						
							while (first_occ_copy < spec.size() && spec[first_occ_copy] <= mass_with_mods + tol_)
							{
											/*
								try
								{
									cerr << "3.>" << s_.substr(indices_[i].first, j + 1) << endl;
								}
								catch(...)
								{
								}*/
								candidates[first_occ_copy].push_back(pnew);
								++first_occ_copy;
							}
						} 
						// if lcp value is bigger than we have walked into the string we add the next entry (indicated by skip vector)
						// isDigestingEnd assures the no wrong hits are added (when we are at the end of a entry and it was no digesting site we must not add the next sequences)
						if ((j+1) <= lcp_[i] && (isDigestingEnd(c, cn) || have_to_go_in))
						{
							for (SignedSize z = 1; z <= skip_[i]; ++z)
							{
								
								char cn_new = ((Size)(indices_[i+z].first + j + 1) == s_.length() - 1) ? 'R' : s_[(indices_[i+z].first) + j + 1];
								if (isDigestingEnd(c, cn_new))
								{
									++nres;
									pair<pair<SignedSize,SignedSize>,DoubleReal> pnew(pair<SignedSize, SignedSize>(indices_[i + z].first, j + 1), found_masses[o]);
									Size first_occ_copy = first_occ;
									while (first_occ_copy<spec.size()&&spec[first_occ_copy] <= mass_with_mods+tol_)
									{
										candidates[first_occ_copy].push_back(pnew);
										++first_occ_copy;
									}
								}
							}
						}
					}
				}
			}
				
			// if we are reaching a lcp postion we add this entry to history
			if (j == (lcp_[i] - 1) && lcp_[i] > 0)
			{
				history.push(pair<pair<SignedSize, map<DoubleReal,SignedSize> >, pair<SignedSize, DoubleReal> >(pair<SignedSize, map<DoubleReal, SignedSize> >(i + skip_[i], map<DoubleReal,SignedSize>(modification_map)), pair<SignedSize, DoubleReal> (j, m)));
			}
			// if mass is to big we can skip the sub tree
			if (m > mmax + tol_) 
			{
				if ((j + 1) <= lcp_[i]) 
				{
					i += skip_[i];
				}
				break;
			}			
		}
	}
	//time_t t1(time(NULL));
	//cout<<"time for search: " << t1 - t0 << endl;
	//cout<<"hits: "<< nres << endl;
	return;
}

void SuffixArrayTrypticCompressed::setTolerance (DoubleReal t)
{
	if (t < 0)
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Tolerance value must not be negative",(String)t);
	}
	tol_ = t;
}

DoubleReal SuffixArrayTrypticCompressed::getTolerance () const 
{
	return (tol_);
}

void SuffixArrayTrypticCompressed::setTags (const vector<String> & tags)
{
	tags_ = tags;
	for (Size i = 0; i < tags.size();++i)
	{
		if (tags.at(i).size()!=3)
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Tag must have size 3", (String)tags.at(i));
		}
	}
	sort(tags_.begin(),tags_.end());
	use_tags_ = true;
}

const vector<String> & SuffixArrayTrypticCompressed::getTags ()
{
	return (tags_);
}

void SuffixArrayTrypticCompressed::setUseTags (bool use_tags)
{
	use_tags_ = use_tags;
	if (tags_.empty()) use_tags_ = false;
}

bool SuffixArrayTrypticCompressed::getUseTags ()
{
	return (use_tags_);
}


void SuffixArrayTrypticCompressed::setNumberOfModifications(Size number_of_mods)
{
	number_of_modifications_ = number_of_mods;
}

Size SuffixArrayTrypticCompressed::getNumberOfModifications ()
{
	return (number_of_modifications_);
}

void SuffixArrayTrypticCompressed::printStatistic ()
{
	progress_=0;
	cout << "Number of suffices: " << indices_.size()<<endl;
	vector<pair<SignedSize,SignedSize> > out_number;
	vector<pair<SignedSize,SignedSize> > edge_length;
	vector<SignedSize> leafe_depth;
	parseTree_(0,indices_.size()-1,1,0,1,out_number,edge_length,leafe_depth);
	for (Size i = 0; i < leafe_depth.size();i++){
		cout<<leafe_depth.at(i)<<",";
	}
	cout<<endl;
	for (Size i = 0; i < out_number.size();i++){
		cout<<"("<<out_number.at(i).first<<","<<out_number.at(i).second<<") ; ";
	}
	cout<<endl;
	for (Size i = 0; i < edge_length.size();i++){
		cout<<"("<<edge_length.at(i).first<<","<<edge_length.at(i).second<<") ; ";
	}
	cout<<endl;
}

void SuffixArrayTrypticCompressed::parseTree_ (SignedSize start_index, SignedSize stop_index, SignedSize depth, SignedSize walked_in, SignedSize edge_len, vector<pair<SignedSize,SignedSize> > & out_number, vector<pair<SignedSize,SignedSize> > & edge_length, vector<SignedSize> & leafe_depth)
{
	//to start walked_in set to 0, depth=1, edge_len = 1
	
	if ((SignedSize)((DoubleReal)leafe_depth.size()/(DoubleReal)indices_.size()*100)>progress_){
		cout<<(DoubleReal)leafe_depth.size()/(DoubleReal)indices_.size()*100<<"%"<<endl;
		progress_++;
	
	}
	

	if (start_index>stop_index) {
		return;
	}
	if (start_index==stop_index) {
		leafe_depth.push_back(depth);
		edge_length.push_back(pair<SignedSize,SignedSize>(depth,indices_.at(start_index).second-walked_in+1));
		return;
	}
	char last_char ='*';
	char actual_char = '*';
	SignedSize start_index_copy = start_index;
	SignedSize number_of_outgoings = 0;
	for (SignedSize i = start_index; i<=stop_index;i++){
		if (indices_.at(i).second<=walked_in){
			leafe_depth.push_back(depth-1);
			++start_index_copy;
		} else {
			if (last_char=='*') {
				last_char=s_[indices_.at(i).first+walked_in];
			}
			actual_char=s_[indices_.at(i).first+walked_in];
			if (actual_char != last_char){
				++number_of_outgoings;
				if (!hasMoreOutgoings_(start_index_copy,i-1,walked_in+1)){
					parseTree_ (start_index_copy,i-1,depth,walked_in+1,edge_len+1,out_number,edge_length,leafe_depth);
				} else {
					parseTree_ (start_index_copy,i-1,depth+1,walked_in+1,1,out_number,edge_length,leafe_depth);
					edge_length.push_back(pair<SignedSize,SignedSize>(depth,edge_len));
				}
				start_index_copy = i;
				last_char=actual_char;
			} else {

			}
		}
	}
	++number_of_outgoings;
	if (!hasMoreOutgoings_(start_index_copy,stop_index,walked_in+1)){
		parseTree_ (start_index_copy,stop_index,depth,walked_in+1,edge_len+1,out_number,edge_length,leafe_depth);
	} else {
		parseTree_ (start_index_copy,stop_index,depth+1,walked_in+1,1,out_number,edge_length,leafe_depth);
		if (number_of_outgoings>1) edge_length.push_back(pair<SignedSize,SignedSize>(depth,edge_len));
	}
	if (number_of_outgoings>1) out_number.push_back(pair<SignedSize,SignedSize>(depth-1,number_of_outgoings));
}

bool SuffixArrayTrypticCompressed::hasMoreOutgoings_ (SignedSize start_index, SignedSize stop_index, SignedSize walked_in){
	SignedSize n_occ =0;
	char last_char ='*';
	char actual_char = '*';
	for (SignedSize i = start_index; i<=stop_index;i++){
		if (indices_.at(i).second<=walked_in){
			n_occ=1;
		} else  {
			if (last_char=='*') {
				last_char=s_[indices_.at(i).first+walked_in];
				++n_occ;
			}
			actual_char=s_[indices_.at(i).first+walked_in];
			if (actual_char != last_char){
				++n_occ;
			} else {

			}
			if (n_occ>1) return (true);
		}
	}
	return false;
}

} // namespace OpenMS

