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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>

namespace OpenMS
{

	BasePairFinder::BasePairFinder()
		: FactoryProduct("BasePairFinder"),
			element_pairs_(0)
	{
		maps_.model_ = 0;
		maps_.scene_ = 0;
		map_index_.model_ = -2;
		map_index_.scene_ = -2;
		transformation_[RawDataPoint2D::RT].setSlope(1);
		transformation_[RawDataPoint2D::RT].setIntercept(0);
		transformation_[RawDataPoint2D::MZ].setSlope(1);
		transformation_[RawDataPoint2D::MZ].setIntercept(0);
	}

	BasePairFinder::~BasePairFinder()
	{
	}

	void BasePairFinder::registerChildren()
  {
    Factory< BasePairFinder>::registerProduct(SimplePairFinder::getProductName(), &SimplePairFinder::create);
    Factory< BasePairFinder>::registerProduct(DelaunayPairFinder::getProductName(), &DelaunayPairFinder::create);
  }

	int BasePairFinder::dumpElementPairs(const String& filename)
  {
    // V_dumpElementPairs() is used for a few comments about the files being
    // written.  We are silent unless output is actually being written, so
    // it is defined here inside the "else" branch.
#define V_dumpElementPairs(bla) std::cerr << bla << std::endl;
    V_dumpElementPairs("### Writing "<<filename);
    std::ofstream dump_file(filename.c_str());
    dump_file << "# " << filename<< " generated " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
    dump_file << "# 1:number 2:quality 3:firstRT 4:firstMZ 5:firstIT 6:firstQual 7:secondRT 8:secondMZ 9:secondIT 10:secondQual\n";
    for ( UInt fp = 0; fp < getElementPairs().size(); ++fp )
    {
      dump_file << fp << ' '
								<< getElementPairs()[fp].getFirst().getRT() << ' '
								<< getElementPairs()[fp].getFirst().getMZ() << ' '
								<< getElementPairs()[fp].getFirst().getIntensity() << ' '
								<< getElementPairs()[fp].getSecond().getRT() << ' '
								<< getElementPairs()[fp].getSecond().getMZ() << ' '
								<< getElementPairs()[fp].getSecond().getIntensity() << ' '
								<< std::endl;
    }
    dump_file << "# " << filename << " EOF " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
    std::string dump_filename_gp = filename + ".gp";
    V_dumpElementPairs("### Writing "<<dump_filename_gp);
    std::ofstream dump_file_gp(dump_filename_gp.c_str());
    dump_file_gp << "# " << dump_filename_gp << " generated " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
    dump_file_gp <<
			"# Gnuplot script to view element pairs\n"
			"plot   \"" << filename <<"\" using 2:3 title \"map 1\"\n"
			"replot \"" << filename <<"\" using 5:6 title \"map 2\"\n"
			"replot \"" << filename <<"\" using 2:3:($5-$2):($6-$3) w vectors nohead title \"pairs\"\n"
			;
    dump_file_gp << "# " << dump_filename_gp << " EOF " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
    V_dumpElementPairs("### You can view `"<<filename<<"' using the command line `gnuplot "<<dump_filename_gp<<" -'");
#undef V_dumpElementPairs

    return 0;
  }
} 
