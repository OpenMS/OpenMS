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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;

namespace OpenMS
{
	
	AScore::AScore()
	{
	}

	AScore::~AScore()
	{
	}

	PeptideHit AScore::compute(PeptideHit& hit, RichPeakSpectrum& real_spectrum, DoubleReal fmt)
	{
		if(real_spectrum.size() == 0)
		{
			return hit;
		}
		//TODO Masse berechnen um zu wissen, ob es sich um multiple handelt oder nicht!!
		TheoreticalSpectrumGenerator spectrum_generator;
		vector<RichPeakSpectrum> th_spectra; //typedef MSSpectrum<RichPeak1D> RichPeakSpectrum;

		for(RichPeakSpectrum::iterator it = real_spectrum.begin(); it < real_spectrum.end(); ++it)
		{
			//cout<<it->getMZ()<< " ";
		}	
		//cout<<endl;
		//cout<<"produce ts"<<endl;	
		///produce theoretical spectra
		{
			const AASequence sequence(hit.getSequence().toUnmodifiedString());
			if((sequence.getNumberOf("S") + sequence.getNumberOf("T") + sequence.getNumberOf("Y")) <= 1) //only compute score if there is more than one possibility
			{
				return hit;
			}
			th_spectra.resize((sequence.getNumberOf("S") + sequence.getNumberOf("T") + sequence.getNumberOf("Y")));
			const String& unmodified = sequence.toUnmodifiedString();
			Size next_th=0;
			//cout<<AASequence(String("T(Phospho)TAT")).getNumberOf("T")<<endl;
			//cout<<"S: "<<sequence.getNumberOf("S") <<" T: "<<sequence.getNumberOf("T") <<" Y: "<<sequence.getNumberOf("Y") <<" total: " <<(sequence.getNumberOf("S") + sequence.getNumberOf("T") + sequence.getNumberOf("Y"))<<endl;
		//	cout<<"#theoretical: "<<th_spectra.size()<<endl;
			for(Size i = 0; i < unmodified.size(); ++i)
			{
				if(unmodified[i] == 'Y' || unmodified[i] == 'T' || unmodified[i] == 'S')
				{
					//cout<<AASequence(String(unmodified.prefix(i+1) + "(Phospho)" + unmodified.suffix(unmodified.size() - i - 1))).toString()<<endl;
					//cout<<AASequence(String(unmodified.prefix(i+1) + "(Phospho)" + unmodified.suffix(unmodified.size() - i - 1))).isValid()<<endl;
				/*	Param param(spectrum_generator.getDefaults());
					param.setValue("a_intensity",0.0);
					param.setValue("x_intensity",0.0);
					param.setValue("c_intensity",0.0);
					param.setValue("z_intensity",0.0);
					param.setValue("relative_loss_intensity",0.0);
					param.setValue("precursor_intensity",0.0);
					param.setValue("precursor_H2O_intensity",0.0);
					param.setValue("precursor_NH3_intensity",0.0);
					spectrum_generator.setParameters(param);*/
					String str(unmodified.prefix(i+1) + "(Phospho)" + unmodified.suffix(unmodified.size() - i - 1));
					spectrum_generator.getSpectrum(th_spectra[next_th], AASequence(str),hit.getCharge() );
					th_spectra[next_th].setName(str);
						//cout<<str<<" , #peaks: "<<th_spectra[next_th].size()<<endl;
						for(RichPeakSpectrum::iterator it = th_spectra[next_th].begin(); it < th_spectra[next_th].end(); ++it)
						{
							//cout<<it->getMZ()<< " ";
						}	
						//cout<<endl;
					//cout<<next_th<<endl;
					//th_spectra[next_th].sortByPosition();
					++next_th;
				}
			}
		}
		///produce theoretical spectra - END
		
		if(!real_spectrum.isSorted())
		{
			real_spectrum.sortByPosition();
		}
		//cout<<"prepare peak depth"<<endl;
		vector< vector<Real> > peptide_site_scores(th_spectra.size());
		vector< RichPeakSpectrum > windows_with_all_peak_depths; 
		///prepare peak depth for all windows in the actual spectrum
		{
			//DPosition<1> highest_peak = real_spectrum.getMax();
		//	cout<<"real_size(): "<<real_spectrum.size()<<endl;
		//	cout<<"groesster?: "<<real_spectrum[real_spectrum.size()-1].getMZ()<<endl;
			//cout<<"highest_peak: "<<highest_peak[0]<<endl;
			Real biggest_window = real_spectrum[real_spectrum.size()-1].getMZ();
			//cout<<"biggest_window: "<<biggest_window<<endl;
			Size number_of_windows = 1;
			if(biggest_window > 100) number_of_windows = ceil(biggest_window/100);
		//	cout<<"number_of_windows: "<<number_of_windows<<endl;
			windows_with_all_peak_depths.resize(number_of_windows);
			RichPeakSpectrum::Iterator begin_window = real_spectrum.MZBegin(0);
			RichPeakSpectrum::Iterator end_window = real_spectrum.MZBegin(100);
			for(Size current_window = 0; current_window < number_of_windows;++current_window)
			{
				RichPeakSpectrum real_window;
				RichPeakSpectrum::Iterator runner = begin_window;
				while(runner <= end_window)
				{
					real_window.push_back(*runner);
					++runner;
				}
				
				real_window.sortByIntensity(true);
				for(Size i = 0; i < 10 ; ++i)
				{
					if(i < real_window.size()) windows_with_all_peak_depths[current_window].push_back(real_window[i]);
				}
				begin_window = ++end_window;
				end_window = real_spectrum.MZBegin((current_window+1) *100);	
			}	
		}
		//cout<<"windows....size: "<<windows_with_all_peak_depths.size()<<endl;
		///prepeare peak depth for all windows in the actual spectrum - END
	//	cout<<"compute matched ions"<<endl;
		UInt N;
		vector< vector<Real> >::iterator side_scores = peptide_site_scores.begin();
		for(vector<RichPeakSpectrum>::iterator it = th_spectra.begin(); it < th_spectra.end(); ++it, ++side_scores) //each theoretical spectrum
		{
			N = it->size();//real oder theo!!
			//cout<<it->getName()<<":  "<<endl;
			UInt i= 1;
			side_scores->resize(10);
			while(i <= 10)
			{		
				//cout<<"depth: "<<i<<endl;
				//Auslagern
				UInt n = 0;
				for(Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth ) //each 100 m/z window 
				{
					UInt nmi = numberOfMatchedIons_(*it, windows_with_all_peak_depths[depth],depth,fmt);
					n += numberOfMatchedIons_(*it, windows_with_all_peak_depths[depth],depth,fmt);
					//cout<<nmi<<" ";
				}
				//cout<<"n: "<<n <<endl;
				//auslagern_end
				Real p = (Real)i/100;
			//	cout<<"p: "<<p <<"/////////////////////////////////////////"<<endl;
		//		cout<<"dudl"<<endl;
				Real cum_socre = computeCumulativeScore(N,n,p);
				//cout<<"Cumulative Score: "<<cum_socre << " avec N: "<<N <<" n: " <<n << " p: "<<p<< " -10*log:" <<(-10* log10(cum_socre))<<endl;
				(*side_scores)[i-1]= (-10* log10(cum_socre));//computeCumulativeScore(N,n,p);
		//		cout<<"wtf"<<endl;
				++i;
			}
		}
	//	cout<<"twohighestpeaks"<<endl;
		Size first,second, peak_depth;
		computeTwoHighestPeptides(peptide_site_scores,first,second,peak_depth);
		//cout<<"first: "<<first<<" ,second: "<<second<<" ,peak_depth: "<<peak_depth<<endl; 
		///Calculate AScore
	//	cout<<"site_determing_ions"<<endl;
		vector<RichPeakSpectrum> site_determining_ions;
		compute_site_determining_ions(hit, first,second, site_determining_ions);
		N = site_determining_ions[0].size(); // all possiblities have the same number
		//cout<<"site_determining_ions"<<endl;
		for(RichPeakSpectrum::iterator it = site_determining_ions[0].begin(); it < site_determining_ions[0].end(); ++it)
		{
			//cout<<it->getMZ()<< " ";
		}
		//cout<<endl;
		for(RichPeakSpectrum::iterator it = site_determining_ions[1].begin(); it < site_determining_ions[1].end(); ++it)
		{
			//cout<<it->getMZ()<< " ";
		}
		//cout<<endl;
		Real p = (Real)peak_depth/100;
		Int n = 0;
	//	cout<<"P_first"<<endl;		
		//Auslagern
		for(Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth ) //each 100 m/z window 
		{
			n += numberOfMatchedIons_(site_determining_ions[0], windows_with_all_peak_depths[depth],depth,fmt);
		}
		//cout<<"N: "<<N << " n: "<< n <<" p: "<<p <<" peak_depth: "<<peak_depth<<endl;
		//auslagern_end
		Real P_first = computeCumulativeScore(N,n,p);
		n = 0;
	//	cout<<"P_second"<<endl;	
		//Auslagern
		for(Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth ) //each 100 m/z window 
		{
			n += numberOfMatchedIons_(site_determining_ions[1], windows_with_all_peak_depths[depth],depth,fmt);
		}
		//auslagern_end
		//cout<<"n: "<<n<<endl;
		//cout<<"N: "<<N << " n: "<< n <<" p: "<<p <<" peak_depth: "<<peak_depth<<endl;
		Real P_second = computeCumulativeScore(N,n,p);
		Real score_first = -10* log10(P_first);
		Real score_second = -10* log10(P_second);
		Real AScore_first = score_first - score_second;
		Real AScore_second =score_second - score_first;
		
		/*	cout<<"****************************************************"<<endl;
			cout<<"Peptide: "<<hit.getSequence().toString()<< " , charge: "<<hit.getCharge() <<" , Precursor: "<<real_spectrum.getPrecursors()[0].getMZ()<<" , MonoWeight: "<< hit.getSequence().getMonoWeight()<<", #peaks: "<<real_spectrum.size()<< " , mascot score: "<<hit.getScore()<<endl;
			cout<<"P_first: "<<P_first<<endl;
			cout<<"P_second: "<<P_second<<endl;
			cout<<"score_first: "<<score_first <<" , score_second: "<<score_second<<endl;
			cout<<"\nfirst: "<<th_spectra[first].getName()<<", Peptide Score: " << score_first<<" , AScore: "<<AScore_first<<endl;
			cout<<"\nsecond: "<<th_spectra[second].getName()<<", Peptide Score: " << score_second<<" , AScore: "<<AScore_second<<endl;*/
		if(AScore_first > AScore_second)
		{
			hit.setMetaValue("PeptideScore",score_first);
			hit.setMetaValue("AScore",AScore_first);
			hit.setMetaValue("Phospho-Sites",th_spectra[first].getName());
		}
		else
		{
			hit.setMetaValue("PeptideScore",score_second);
			hit.setMetaValue("AScore",AScore_second);
			hit.setMetaValue("Phospho-Sites",th_spectra[second].getName());		
		}
		return hit;
		
		/****/
	/*	Real highest_peak = real_spectrum.getMax().getIntensity();
		Int biggest_window = ceil(highest_peak);
		Int number_of_windows = 1;
		if(biggest_window > 100) number_of_windows = biggest_window/100;
		Int N = real_spectrum.size();
		UInt i = 1;
		MSSpectrum<>::Iterator begin_window = real_spectrum.MZBegin(0);
		MSSpectrum<>::Iterator end_window = real_spectrum.MZBegin(100);
		while(i <= 10)
		{
			Int n = 0;
			for(Int current_window = 1; current_window <= biggest_window;++current_window)
			{
				MSSpectrum<> window(begin_window, end_window);
				window.sortIntensity(true);
				Int matched_ions = 0;
				for(UInt match = 0 ; match < i && match < window.size();++match)
				{
					Size nearest_peak = th_spectrum.findNearest(window[match].getMZ());
					if(fabs(th_spectrum[nearest_peak].getMZ() - window[match].getMZ()) < fmt) ++matched_ions;
				}
				n += matched_ions;
				
				begin_window = ++end_window;
				end_window = real_spectrum.MZBegin(current_window*100);
			}
			Real p = i/100;
			Real P = ...
			begin_window = real_spectrum.MZBegin(0)		;
			end_window = real_spectrum.MZBegin(100);			
			++i;
		}*/
		
//!erstmal das Spektrum in 100 m/z Fenster aufteilen.

//!N = total number of fragment ions for the peptide
//!i = 1 // peak depth
//!while(i<= 10)				// peak depth testen von 1 bis 10 so wie in Figure 3b				
//!	n =  0  // number of ions matched to the spectrum
//!	for_each mass window
//!		overlay top i peaks with predicted b- and y- ions(by intensity)
//!		n = matched_ions + n
//!	p = i/100
//!	P_i(X)  =  sum_(k=n)^N  (N über k) p^k(1-p)^(N-k) // am besten auf Seite 1286 nochmal nachschauen :D
//!	Score_i = -10 * log(P_i(X))
//!	i++

//!Am Ende kann man Score_i wie auf Seite 1291unten zusammenfügen. // 1)

//!Dann nimmt man die zwei Möglichkeiten mit dem höchsten score, wie bei 1) ausgerechnet.

//!Eigentlicher AScore:

//!Man nimmt von den besten Möglichkeiten des Peptides die "peak depth", bei der Score unterschied am größten ist.
//!Und berechnet P(X) dafür nochmals. Allerdings mit dem Unterschied, dass man nur die "site-determining ions"  verwendet.
		
		
	}
	
	Real AScore::computeCumulativeScore(UInt N,UInt n, Real p)
	{
	//	cout<<"1 ";
		Real score = 0.0;
	//	cout<<"2 ";
		if(n > N)	return score;
	//	cout<<"3 ";
		for(UInt k = n; k <= N ; ++k)
		{
		//	cout<<"p: "<<p<<" ,k: "<<k<<" , "<<boost::math::binomial_coefficient<Real>(N, k)<<"*"<<pow((double)p,(double)k)<<"*"<<pow(double(1-p),double(N-k))<<"="<<boost::math::binomial_coefficient<Real>(N, k)*pow((double)p,(double)k)*pow(double(1-p),double(N-k))<<endl;
		//	cout<<"3.1 ";
		//	cout.flush();
			DoubleReal coeff = boost::math::binomial_coefficient<DoubleReal>((Real)N, k);
		//	cout<<"3.2 ";
		//	cout.flush();
			Real pow1 = pow((double)p,(double)k);
		//	cout<<"3.3 ";
		//	cout.flush();
			Real pow2 = pow(double(1-p),double(N-k));
		//	cout<<"3.4 ";
		//	cout.flush();
			score += coeff*pow1*pow2;
		//	cout<<"3.5 ";
		//	cout.flush();
		}
		//cout<<"4 ";
	//	cout<<"cumulativeScore: "<<score<<endl;
		return score;
	}
	
	void AScore::computeTwoHighestPeptides( vector< vector<Real> >& peptide_site_scores,Size& first,Size& second, Size& peak_depth)
	{
		first = 0;	
		second = 0;
		Real first_score = 0.0;
		Real second_score = 0.0;
		peak_depth = 1;
		for(Size it = 0;it < peptide_site_scores.size(); ++it)
		{
			Real current_score = (peptide_site_scores[it][0]*0.5
														+peptide_site_scores[it][1]*0.75
														+peptide_site_scores[it][2]//*1
														+peptide_site_scores[it][3]//*1
														+peptide_site_scores[it][4]//*1
														+peptide_site_scores[it][5]//*1
														+peptide_site_scores[it][6]*0.75
														+peptide_site_scores[it][7]*0.5
														+peptide_site_scores[it][8]*0.25
														+peptide_site_scores[it][9]*0.25)
														/10;
	//		cout<<"current_score:"<< current_score<<endl;
			if(current_score > first_score)
			{
				second_score = first_score;
				first_score = current_score;
				second = first;
				first = it;
			
			}
			else if(current_score > second_score)
			{
				second_score = current_score;
				second = it;
			}
		}
		if(second == first) ++second;
		Real current_peak_depth = 0.0;
		vector<Real>::iterator first_it = peptide_site_scores[first].begin();
		Size depth = 1;
		//cout<<"---"<<endl;
		for(vector<Real>::iterator second_it = peptide_site_scores[second].begin(); second_it < peptide_site_scores[second].end(); ++second_it, ++first_it)
		{
		//	cout<<"current: "<<current_peak_depth<<" ,first: "<<*first_it<<" ,second: "<<*second_it <<" ,fabs: "<<fabs(*first_it - *second_it)<<endl;
			if(fabs(*first_it - *second_it) > current_peak_depth)
			{
				current_peak_depth = fabs(*first_it - *second_it);
				peak_depth = depth;
				//cout<<current_peak_depth<<endl;
			}
			++depth;
		}
	}
	void AScore::compute_site_determining_ions(PeptideHit& hit, Size first,Size second, vector<RichPeakSpectrum>& site_determining_ions)
	{
		const String& unmodified = hit.getSequence().toUnmodifiedString();
		Size first_AA = 0;
		Size second_AA = 0;
		Size count = 0;
		for(Size i = 0; i < unmodified.size(); ++i)
		{
			if(unmodified[i] == 'Y' || unmodified[i] == 'T' || unmodified[i] == 'S')
			{

				if(count == first) first_AA = i;
				if(count == second) second_AA = i;
				++count;
			}
		}
		site_determining_ions.resize(2);
		TheoreticalSpectrumGenerator spectrum_generator;
	/*	Param param(spectrum_generator.getDefaults());
		param.setValue("a_intensity",0.0);
		param.setValue("x_intensity",0.0);
		param.setValue("c_intensity",0.0);
		param.setValue("z_intensity",0.0);
		param.setValue("relative_loss_intensity",0.0);
		param.setValue("precursor_intensity",0.0);
		param.setValue("precursor_H2O_intensity",0.0);
		param.setValue("precursor_NH3_intensity",0.0);
		spectrum_generator.setParameters(param);		*/
		//cout<<"Site determining ions von: "<<unmodified<<endl;
		//cout<<"first: "<<first<<" , " <<first_AA<<" , second: "<<second <<" , "<<second_AA<<endl;
		if(first_AA < second_AA)
		{
			AASequence pref(unmodified.prefix(first_AA));
			AASequence suf(unmodified.suffix(unmodified.size()-second_AA-1));
			AASequence pref_with_phospho_first(unmodified.prefix(first_AA+1) + "(Phospho)"+ unmodified.substr(first_AA+1,second_AA - first_AA) );
			AASequence pref_with_phospho_second(unmodified.prefix(first_AA+1)+ unmodified.substr(first_AA+1,second_AA - first_AA) + "(Phospho)" );
			AASequence suf_with_phospho_first(unmodified.substr(first_AA,1)+"(Phospho)" + unmodified.suffix(unmodified.size()-first_AA-1));
			AASequence suf_with_phospho_second(unmodified.substr(first_AA,second_AA - first_AA+1) +"(Phospho)"+ unmodified.suffix(unmodified.size()-second_AA-1));
			//cout<<"prefix: "<<pref.toString()<<endl;
			//cout<<"suffix: "<<suf.toString()<<endl;
			//cout<<"pref_with_phospho_first: "<<pref_with_phospho_first.toString()<<" pref_with_phospho_second: "<<pref_with_phospho_second.toString()<<endl;
			//cout<<"suff_with_phospho_first: "<<suf_with_phospho_first.toString()<<"suff_with_phospho_second: "<<suf_with_phospho_second.toString()<<endl;
			RichPeakSpectrum prefix, suffix,prefix_with_phospho_first, prefix_with_phospho_second, suffix_with_phospho_second, suffix_with_phospho_first;
			spectrum_generator.addPeaks(prefix, pref, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix, suf, Residue::YIon, hit.getCharge());
			spectrum_generator.addPeaks(prefix_with_phospho_first, pref_with_phospho_first, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(prefix_with_phospho_second, pref_with_phospho_second, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix_with_phospho_first, suf_with_phospho_first, Residue::YIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix_with_phospho_second, suf_with_phospho_second, Residue::YIon,hit.getCharge());
			if(prefix.size() > 0)
			{
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
				{
					if(it->getMZ() > prefix[prefix.size()-1].getMZ()) site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
				{
					if(it->getMZ() > prefix[prefix.size()-1].getMZ()) site_determining_ions[1].push_back(*it);
				}
			}
			else
			{
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
				{
					site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
				{
					site_determining_ions[1].push_back(*it);
				}				
			}
			if(suffix.size() > 0)
			{
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_first.begin(); it < suffix_with_phospho_first.end(); ++it)
				{
					if(it->getMZ() > suffix[suffix.size()-1].getMZ()) site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
				{
					if(it->getMZ() > suffix[suffix.size()-1].getMZ()) site_determining_ions[1].push_back(*it);
				}
			}
			else
			{
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_first.begin(); it < suffix_with_phospho_first.end(); ++it)
				{
					site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
				{
					site_determining_ions[1].push_back(*it);
				}
			}
			
			//cout<<"größe: "<<site_determining_ions[0].size() <<" "<<site_determining_ions[1].size()<<endl;
	//		spectrum_generator.getSpectrum(site_determining_ions[0], AASequence(String(unmodified.substr(first_AA,1)+ "(Phospho)" + unmodified.substr(first_AA+1,second_AA - first_AA))),hit.getCharge() );
	//		cout<<unmodified.substr(first_AA,second_AA - first_AA+1) + "(Phospho)"<<endl;
	//		spectrum_generator.getSpectrum(site_determining_ions[1], AASequence(String(unmodified.substr(first_AA,second_AA - first_AA+1) + "(Phospho)")),hit.getCharge() );
		}
		else
		{
			AASequence pref(unmodified.prefix(second_AA));
			AASequence suf(unmodified.suffix(unmodified.size()-first_AA-1));
			AASequence pref_with_phospho_second(unmodified.prefix(second_AA+1) + "(Phospho)"+ unmodified.substr(second_AA+1,first_AA - second_AA) );
			AASequence pref_with_phospho_first(unmodified.prefix(second_AA+1)+ unmodified.substr(second_AA+1,first_AA - first_AA) + "(Phospho)" );
			AASequence suf_with_phospho_second(unmodified.substr(second_AA,1)+"(Phospho)" + unmodified.suffix(unmodified.size()-second_AA-1));
			AASequence suf_with_phospho_first(unmodified.substr(second_AA,first_AA - second_AA+1) +"(Phospho)"+ unmodified.suffix(unmodified.size()-first_AA-1));
			//cout<<"prefix: "<<pref.toString()<<endl;
			//cout<<"suffix: "<<suf.toString()<<endl;
			//cout<<"pref_with_phospho_first: "<<pref_with_phospho_first.toString()<<" pref_with_phospho_second: "<<pref_with_phospho_second.toString()<<endl;
			//cout<<"suff_with_phospho_first: "<<suf_with_phospho_first.toString()<<"suff_with_phospho_second: "<<suf_with_phospho_second.toString()<<endl;
			RichPeakSpectrum prefix, suffix,prefix_with_phospho_first, prefix_with_phospho_second, suffix_with_phospho_second, suffix_with_phospho_first;
			spectrum_generator.addPeaks(prefix, pref, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix, suf, Residue::YIon, hit.getCharge());
			spectrum_generator.addPeaks(prefix_with_phospho_first, pref_with_phospho_first, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(prefix_with_phospho_second, pref_with_phospho_second, Residue::BIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix_with_phospho_first, suf_with_phospho_first, Residue::YIon,hit.getCharge());
			spectrum_generator.addPeaks(suffix_with_phospho_second, suf_with_phospho_second, Residue::YIon,hit.getCharge());
			if(prefix.size() > 0)
			{
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
				{
					if(it->getMZ() > prefix[prefix.size()-1].getMZ()) site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
				{
					if(it->getMZ() > prefix[prefix.size()-1].getMZ()) site_determining_ions[1].push_back(*it);
				}
			}
			else
			{
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
				{
					site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
				{
					site_determining_ions[1].push_back(*it);
				}				
			}
			if(suffix.size() > 0)
			{
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_first.begin(); it < suffix_with_phospho_first.end(); ++it)
				{
					if(it->getMZ() > suffix[suffix.size()-1].getMZ()) site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
				{
					if(it->getMZ() > suffix[suffix.size()-1].getMZ()) site_determining_ions[1].push_back(*it);
				}
			}
			else
			{
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_first.begin(); it < suffix_with_phospho_first.end(); ++it)
				{
					site_determining_ions[0].push_back(*it);
				}
				for(RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
				{
					site_determining_ions[1].push_back(*it);
				}
			}			
			
			//spectrum_generator.getSpectrum(site_determining_ions[0], AASequence(String(unmodified.substr(second_AA,1) + "(Phospho)" + unmodified.substr(second_AA+1,first_AA - second_AA))),hit.getCharge() );
			//spectrum_generator.getSpectrum(site_determining_ions[1], AASequence(String(unmodified.substr(second_AA,first_AA - second_AA+1) + "(Phospho)")),hit.getCharge() );
		}
	}
	Int AScore::numberOfMatchedIons_(const RichPeakSpectrum& th,const RichPeakSpectrum& windows ,Size depth, DoubleReal fmt)
	{
		Int n = 0;
		for(Size i = 0; i < windows.size() && i < depth; ++i)
		{
				Size nearest_peak = th.findNearest(windows[i].getMZ());
				if(nearest_peak < th.size() && fabs(th[nearest_peak].getMZ() - windows[i].getMZ()) < fmt) ++n;
		}
		return n;
	}
	
} // namespace OpenMS

