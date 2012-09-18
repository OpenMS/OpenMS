// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

using namespace std;

namespace OpenMS
{
	
	PILISScoring::PILISScoring()
		:	DefaultParamHandler("PILISScoring")
	{
		defaults_.setValue("use_local_scoring", 1, "If set to 1, a E-Value of an identification run of one spectrum is used additionally");
		defaults_.setValue("survival_function_bin_size", 20, "Bin size of the survival function", StringList::create("advanced"));
		defaults_.setValue("global_linear_fitting_threshold", 0.1, "Fitting threshold of the survival function of the global E-Value calculation", StringList::create("advanced"));
		defaults_.setValue("local_linear_fitting_threshold", 0.5, "Fitting threshold of the survival function of the local E-Value calculation", StringList::create("advanced"));
		defaults_.setValue("score_default_value", 10e10, "If no score can be assigned use this one", StringList::create("advanced"));
		defaultsToParam_();
	}

	PILISScoring::PILISScoring(const PILISScoring& rhs)
		: DefaultParamHandler(rhs)
	{
	}

	PILISScoring& PILISScoring::operator = (const PILISScoring& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator=(rhs);
		}
		return *this;
	}

	PILISScoring::~PILISScoring()
	{
	}

	void PILISScoring::getScore(PeptideIdentification& id)
	{
    if (id.getHits().empty())
    {
      return;
    }

    if (id.getHits().size() > 2)
    {
      vector<double> scores;
      vector<PeptideHit>::const_iterator it = id.getHits().begin();
      for (++it; it != id.getHits().end(); ++it)
      {
        scores.push_back(it->getScore());
      }

      double slope(0);
      double intercept(0);

      getFitParameter_(slope, intercept, scores, (double)param_.getValue("local_linear_fitting_threshold"));

      if (slope != 0 && intercept != 0)
      {
				id.setScoreType("PILIS-E-value");
				vector<PeptideHit> tmp_hits = id.getHits();
        for (vector<PeptideHit>::iterator it = tmp_hits.begin(); it != tmp_hits.end(); ++it)
        {
          double evalue = exp(intercept + slope * log(it->getScore()));
          it->setScore(evalue);
        }
				id.setHits(tmp_hits);
      }
    }
	}

	void PILISScoring::getScores(vector<PeptideIdentification>& ids)
	{
		// get all but the first scores
  	vector<double> global_scores;
  	for (vector<PeptideIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
  	{
			if (it->getHits().empty())
			{
				break;
			}
    	vector<PeptideHit>::const_iterator it1 = it->getHits().begin();
    	for (++it1; it1 != it->getHits().end(); ++it1)
    	{
      	global_scores.push_back(it1->getScore());
    	}
  	}

		// get the fit parameter for the global survival function
		double global_slope = 0;
		double global_intercept = 0;
		getFitParameter_(global_slope, global_intercept,  global_scores, (double)param_.getValue("global_linear_fitting_threshold"));

		// annotate the ProteinIdentification with both scores (global and single identification)
		for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
		{
			getScore_(*it, global_slope, global_intercept);
		}

		return;
	}

	void PILISScoring::getScore_(PeptideIdentification& id, double global_slope, double global_intercept)
	{
		if (id.getHits().empty())
		{
			return;
		}

		bool use_local_scoring(true);
		if ((Size)param_.getValue("use_local_scoring") == 0)
		{
			use_local_scoring = false;
		}

		// if possible and allowed using local scoring 
		if (id.getHits().size() > 2 && use_local_scoring)
		{
			vector<double> scores;
			vector<PeptideHit>::const_iterator it = id.getHits().begin();
			for (++it; it != id.getHits().end(); ++it)
			{
				scores.push_back(it->getScore());
			}
		
			double slope(0);
			double intercept(0);

			getFitParameter_(slope, intercept, scores, (double)param_.getValue("local_linear_fitting_threshold"));

			if (slope != 0 && intercept != 0)
			{
				id.setScoreType("PILIS-E-value");
				vector<PeptideHit> tmp_hits = id.getHits();
				for (vector<PeptideHit>::iterator it = tmp_hits.begin(); it != tmp_hits.end(); ++it)
				{
					double local_evalue = exp(intercept + slope * log(it->getScore()));
					double global_evalue = exp(global_intercept + global_slope * log(it->getScore()));
					it->setScore(local_evalue + global_evalue);
				}
				id.setHits(tmp_hits);
			}
			else
			{
				double score_default_value = (double)param_.getValue("score_default_value");
				id.setScoreType("PILIS-E-value");
				vector<PeptideHit> tmp_hits = id.getHits();
				for (vector<PeptideHit>::iterator it = tmp_hits.begin(); it != tmp_hits.end(); ++it)
				{
					it->setScore(score_default_value);
				}
				id.setHits(tmp_hits);
			}
		}
		else
		{
			if (global_intercept != 0 && global_slope != 0)
			{
				id.setScoreType("PILIS-E-value");
				vector<PeptideHit> tmp_hits = id.getHits();
				for (vector<PeptideHit>::iterator it = tmp_hits.begin(); it != tmp_hits.end(); ++it)
				{
					it->setScore(exp(global_intercept + global_slope * log(it->getScore())));
				}
				id.setHits(tmp_hits);
			}
			else
			{
				double score_default_value = (double)param_.getValue("score_default_value");
				id.setScoreType("PILIS-E-value");
				vector<PeptideHit> tmp_hits = id.getHits();
				for (vector<PeptideHit>::iterator it = tmp_hits.begin(); it != tmp_hits.end(); ++it)
				{
					it->setScore(score_default_value);
				}
				id.setHits(tmp_hits);
			}
		}
	}

	void PILISScoring::getFitParameter_(double& slope, double& intercept, const vector<double>& scores, double threshold)
	{
		slope = 0;
		intercept = 0;
		
		double survival_function_bin_size = (double)param_.getValue("survival_function_bin_size");

  	Map<UInt, double> score_dist_discrete;
  	for (vector<double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
  	{
    	UInt bin = (UInt)((*it) * survival_function_bin_size);
    	if (score_dist_discrete.has(bin))
    	{
      	score_dist_discrete[bin] += 1;
    	}
    	else
    	{
      	score_dist_discrete[bin] = 1;
    	}
  	}

  	vector<DPosition<2> > survival_function;
  	getSurvivalFunction_(score_dist_discrete, survival_function);

  	// fit the high scoring part of the survival function linearly
  	vector<double> x_values;
  	vector<double> y_values;

  	for (vector<DPosition<2> >::const_iterator sit = survival_function.begin(); sit != survival_function.end(); ++sit)
  	{
    	if (sit->getY() < threshold) 
    	{
      	x_values.push_back(log(sit->getX()));
      	y_values.push_back(log(sit->getY()));
    	}
  	}

		Math::LinearRegression lin_reg;
  	if (x_values.size() > 2)
  	{
    	lin_reg.computeRegression(0.95, x_values.begin(), x_values.end(), y_values.begin());

    	slope = lin_reg.getSlope();
    	intercept = lin_reg.getIntercept();
		}
	}

	void PILISScoring::getSurvivalFunction_(Map<UInt, double>& points, vector<DPosition<2> >& survival_function)
	{
  	// normalize the score density
  	double sum(0);
  	vector<UInt> indices;
  	for (Map<UInt, double>::ConstIterator it = points.begin(); it != points.end(); ++it)
  	{
    	sum += it->second;
    	indices.push_back(it->first);
  	}
  	for (Map<UInt, double>::Iterator it = points.begin(); it != points.end(); ++it)
  	{
    	it->second /= sum;
  	}

		double survival_function_bin_size = (double)param_.getValue("survival_function_bin_size");
  	sort(indices.begin(), indices.end());
  	for (Size i = 0; i != indices.size(); ++i)
  	{
    	//cerr << indices[i] << " ";
    	sum = 0;
    	for (Size j = i; j != indices.size(); ++j)
    	{
      	sum += points[indices[j]];
    	}
    	DPosition<2> pos;
    	pos.setX((double)indices[i] / survival_function_bin_size);
    	//cerr << (double)points[indices[i]] << endl;
    	pos.setY(sum);
    	survival_function.push_back(pos);
  	}

  	return;
	}
}

