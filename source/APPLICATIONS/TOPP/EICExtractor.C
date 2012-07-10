// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/FORMAT/EDTAFile.h>

#include <functional>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_EICExtractor EICExtractor

	@brief Extracts EICs from an MS experiment, in order to quantify analytes at a given position

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ EICExtractor \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> statistical tools, e.g., Excel, R, ... </td>
		</tr>
	</table>
</CENTER>

  Use this instead of FeatureFinder, if you have bad features  which are not recognized (much noise etc)
	or if you want to quantify non-peptides.	

	The EDTA file will specify where to search for signal. 
	Retention time is in seconds [s].
	'int' and 'charge' are ignored but need to be present. However, you MUST specify a 'rank' column. Rows with equal rank are
	summed up in intensity (e.g. useful if you have charge variants you want to sum up to enhance quantitation robustness).
	Each rank represents a so called Master Compound, which constists of one or more sub compounds.
	
Example:<br>
	<pre>
RT	m/z	int	charge	rank	
19.2	431.8599024	0	0	1	
21	678.7729237	0	0	2
25	660.7629237	0	0	2
59.2	431.8599024	0	0	3
</pre>

	Here, rows 2 and 3 will be summed up, as they have the same rank.

	As output, two files in text format are given. The detail file gives RT and m/z deltas from expected to identified signal position etc, the sum file
	represents the master compounds.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_EICExtractor.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_EICExtractor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPEICExtractor
      : public TOPPBase
{
  public:
    TOPPEICExtractor()
        : TOPPBase("EICExtractor","Extracts intensities from dedicates positions in a LC/MS map")
    {
    }

    void registerOptionsAndFlags_()
    {
	  	registerInputFileList_("in","<file>", StringList::create(""),"Input raw data file");
			setValidFormats_("in",StringList::create("mzML"));
      registerInputFile_("pos", "<file>", "", "Input config file stating where to find signal");
      setValidFormats_("pos",StringList::create("edta"));
      registerDoubleOption_("rt_tol", "", 3, "RT tolerance in [s] for finding max peak (whole RT range around RT middle)", false, false);
      registerDoubleOption_("mz_tol", "", 10, "m/z tolerance in [ppm] for finding a peak", false, false);
      registerIntOption_("rt_collect", "", 1, "# of scans up & down in RT from highest point for ppm estimation in result", false, false);
			registerOutputFile_("out","<file>", "", "Output quantitation file (summed intensities by master compounds)");
      //setValidFormats_("out", StringList::create("txt")); // 'txt' not supported yet
      registerOutputFile_("out_detail","<file>", "", "Output quantitation file");
      //setValidFormats_("out_detail", StringList::create("txt")); // 'txt' not supported yet
    }

    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      StringList in = getStringList_("in");
      String edta = getStringOption_("pos");
      String out = getStringOption_("out");
      String out_detail = getStringOption_("out_detail");
      
      DoubleReal rttol = getDoubleOption_("rt_tol");
      DoubleReal mztol = getDoubleOption_("mz_tol");
      Size rt_collect = getIntOption_("rt_collect");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------
      MzMLFile mzml_file;
      mzml_file.setLogType(log_type_);
      MSExperiment<Peak1D> exp;
      
      EDTAFile ed;
      ConsensusMap cm;
      ed.load(edta, cm);

      TextFile tf_master; // one line per master-compound, one intensity column per experiment
      TextFile tf_single; // one line for each compound, three columns for each experiment
      
      tf_master.resize(1); // for header line
      tf_single.resize(cm.size()+2); // two header lines: #1 for filenames; #2 for dRT,ppm, intensity
      tf_single[0] = "#filenames";
      tf_single[1] = "rank";
      for (Size i=0; i<cm.size(); ++i)
      {
        if (!cm[i].metaValueExists("rank"))
        {
          LOG_FATAL_ERROR << "Required column 'rank' not found in EDTA file. Aborting ...\n";
          return ILLEGAL_PARAMETERS;
        }
        Size rank;
        try
        {
          rank = String(cm[i].getMetaValue("rank")).toInt();
        }
        catch (Exception::ConversionError& /*e*/)
        {
          LOG_FATAL_ERROR << "Entry in column 'rank' (line " << i << ") is not a valid integer! Aborting ...\n";
          return ILLEGAL_PARAMETERS;
        }
        tf_single[i+2] += rank; // rank column before first experiment
      }

      for (Size fi=0; fi<in.size();++fi)
      {
        mzml_file.load(in[fi], exp);

			  if (exp.empty())
			  {
				  LOG_WARN << "The given file does not contain any conventional peak data, but might"
					            " contain chromatograms. This tool currently cannot handle them, sorry.";
				  return INCOMPATIBLE_INPUT_DATA;
			  }


        Map<Size, DoubleReal> quant;
      
        tf_single[0] += "\t" + File::basename(in[fi]) + "\t\t";
        tf_single[1] += "\tdRT\tppm\tint";

        // search for each EIC and add up
        Int not_found(0);
        for (Size i=0; i<cm.size(); ++i)
        {
          Size rank = String(cm[i].getMetaValue("rank")).toInt();

          //std::cerr << "Rt" << cm[i].getRT() << "  mz: " << cm[i].getMZ() << " R " <<  cm[i].getMetaValue("rank") << "\n";

          DoubleReal mz_da = mztol*cm[i].getMZ()/1e6; // mz tolerance in Dalton
          MSExperiment<>::ConstAreaIterator it = exp.areaBeginConst(cm[i].getRT()-rttol/2,
                                                                     cm[i].getRT()+rttol/2,
                                                                     cm[i].getMZ()-mz_da, 
                                                                     cm[i].getMZ()+mz_da);
          Peak2D max_peak;
          max_peak.setIntensity(0);
          max_peak.setRT(cm[i].getRT());
          max_peak.setMZ(cm[i].getMZ());
          for (; it != exp.areaEndConst(); ++it)
          {
            if (max_peak.getIntensity() < it->getIntensity())
            {
              max_peak.setIntensity(it->getIntensity());
              max_peak.setRT(it.getRT());
              max_peak.setMZ(it->getMZ());
            }
          }
          DoubleReal ppm = 0; // observed m/z offset
          DoubleReal q=0; // result of quantitation to store

          if (max_peak.getIntensity() == 0)
          {
            ++not_found;
          }
          else
          {
            // take median for m/z found 
            std::vector<DoubleReal> mz;
            MSExperiment<>::Iterator itm = exp.RTBegin(max_peak.getRT());
            SignedSize low = std::min<SignedSize>(std::distance(exp.begin(), itm), rt_collect);
            SignedSize high = std::min<SignedSize>(std::distance(itm, exp.end())-1, rt_collect);
            MSExperiment<>::AreaIterator itt = exp.areaBegin( (itm - low)->getRT()-0.01, (itm + high)->getRT()+0.01, cm[i].getMZ()-mz_da, cm[i].getMZ()+mz_da);
            for (; itt != exp.areaEnd(); ++itt)
            {
              mz.push_back(itt->getMZ());
            }
          
            if ((SignedSize)mz.size() > (low+high+1)) LOG_WARN << "Compound " << i << " has overlapping peaks [" << mz.size() << "/" << low+high+1 << "]\n";
          
            if ( !mz.empty() )
            {
              DoubleReal avg_mz = std::accumulate(mz.begin(), mz.end(), 0.0) / DoubleReal(mz.size());
              ppm = (avg_mz - cm[i].getMZ())/cm[i].getMZ() * 1e6;
            }

            // intensity:
            q = max_peak.getIntensity(); // max peak
/*
            // .. + left & right shoulders
            Int rt_shape_check = 50;
            std::vector< std::pair < DoubleReal, DoubleReal > > rt_shape_left, rt_shape_right;
            low = std::min<SignedSize>(std::distance(exp.begin(), itm), rt_shape_check);
            high = std::min<SignedSize>(std::distance(itm, exp.end())-1, rt_shape_check);
            itt = exp.areaBegin( (itm - low)->getRT()-0.01, (itm + high)->getRT()+0.01, cm[i].getMZ()-mz_da, cm[i].getMZ()+mz_da);
            for (; itt != exp.areaEnd(); ++itt)
            {
              if (itt.getRT() < max_peak.getRT()) rt_shape_left.push_back(make_pair(itt.getRT(), itt->getIntensity()));
              if (itt.getRT() > max_peak.getRT()) rt_shape_right.push_back(make_pair(itt.getRT(), itt->getIntensity()));
            }

            
            if ((SignedSize)(rt_shape_left.size() + rt_shape_right.size()) > (low+high)) LOG_WARN << "Compound " << i << " has overlapping RT peaks [" << (rt_shape_left.size() + rt_shape_right.size())<< "/" << low+high << "]\n";
            if ((rt_shape_left.size() + rt_shape_right.size()) > 0)
            {
              // go from RT max to left and right until 50% threshold reached
              DoubleReal thresh = max_peak.getIntensity() * 0.1;
              // go down shoulders:
              for (std::vector< std::pair < DoubleReal, DoubleReal > >::const_reverse_iterator itp = rt_shape_left.rbegin(); itp != rt_shape_left.rend(); ++itp)
              {
                if (itp->second < thresh) break;
                q += itp->second;
                ++pcount;
              }
              for (std::vector< std::pair < DoubleReal, DoubleReal > >::const_iterator itp = rt_shape_right.begin(); itp != rt_shape_right.end(); ++itp)
              {
                if (itp->second < thresh) break;
                q += itp->second;
                ++pcount;
              }
            }
*/

          }

          quant[rank] += q;
        
          tf_single[i+2] += "\t" + String(max_peak.getRT() - cm[i].getRT()) + "\t" + String(ppm)  + "\t" + String(q);
        }

        LOG_INFO << "No peaks for " << not_found << " compounds in file '" << in[fi] << "'.\n";

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        if (fi!=0) tf_master[0] += "\t";
        tf_master[0] += "sum_" + File::basename(in[fi]);
        int line(0);
			  for (Map<Size, DoubleReal>::const_iterator it=quant.begin(); it!=quant.end(); ++it)
			  {
          String data = /*String(it->first) + "\t" +*/ String(it->second);
          if (fi==0) tf_master.push_back(data);
          else tf_master[++line] += "\t" + data;
			  }

      }

      tf_master.store(out);
      tf_single.store(out_detail);

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPEICExtractor tool;
  return tool.main(argc,argv);
}

/// @endcond
