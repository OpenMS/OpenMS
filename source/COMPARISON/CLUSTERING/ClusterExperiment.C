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
// $Maintainer:  $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>

#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/FORMAT/DataSetInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/COMPARISON/CLUSTERING/helper.h>
#include <OpenMS/COMPARISON/SPECTRA/SequestCompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/LinkageCluster.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperimentXMLHandler.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>


using namespace std;

namespace OpenMS
{
	ClusterExperiment::CanNotRun::CanNotRun(const char* file, int line, const char* function) throw()
    :Base(file, line, function, "ClusterExperiment::CanNotRun","ClusterExperiment::ClusterRun is incomplete")
  {
  }

  ClusterExperiment::CanNotRun::CanNotRun(const char* message,const char* file, int line, const char* function) throw()
    :Base(file, line, function, "ClusterExperiment::CanNotRun",message)
  {
  }
  
  ClusterExperiment::CanNotRun::~CanNotRun() throw()
  {
  }
  
  ClusterExperiment::NoClusterRun::NoClusterRun(const char* file, int line, const char* function) throw()
    :Base(file, line, function, "ClusterExperiment::NoClusterRun","ClusterExperiment does not contain a ClusterRun for the requested operation")
  {
  }
  
  ClusterExperiment::NoClusterRun::~NoClusterRun() throw()
  {
  }

  ClusterExperiment::ClusterExperiment()
    :adapterp_(0), datasetname_(""), dsip_(0) ,info_user_("") , info_date_(""), info_comment_(""), runps_(),currentrun_(-1) 
  {
  }

  ClusterExperiment::~ClusterExperiment()
  {
    delete dsip_;
    for (uint i = 0; i < runps_.size(); ++i)
    {
      delete runps_[i];
    }
  }

  //todo? change owner in ClusterRun?
  ClusterExperiment::ClusterExperiment(const ClusterExperiment& source)
    :adapterp_(source.adapterp_), datasetname_(source.datasetname_),dsip_(new DataSetInfo(*source.dsip_)),info_user_(source.info_user_), info_date_(source.info_date_), info_comment_(source.info_comment_), runps_(source.runps_),currentrun_(source.currentrun_)
  {
  }

  ClusterExperiment& ClusterExperiment::operator = ( const ClusterExperiment& source)
  {
    adapterp_=source.adapterp_;
    datasetname_ = source.datasetname_;
    dsip_ = new DataSetInfo(*source.dsip_);
    info_user_=source.info_user_;
    info_date_=source.info_date_;
    info_comment_=source.info_comment_;
    runps_=source.runps_;
    currentrun_=source.currentrun_;
    return *this;
  }
  
  void ClusterExperiment::deleteContents()
  {
    for (uint i = 0; i < runps_.size() ; ++i)
    {
      runps_[i]->deleteContents();
      delete runps_[i];
    }
    delete dsip_;
    delete adapterp_;
  }

  void ClusterExperiment::save(const String& file)
  {
    ofstream document(file.c_str());
    int ind = 0;
    document << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n";
    document << "<ClusterExperiment>\n";
    document << indent(++ind);
    document << "<Data ";
    document << " Name = \"" << getDataSetInfo()->name()  << "\"";
    document << " Info = \"" << getDataSetInfo()->info() << "\"";
    document << " Size = \"" << getDataSetInfo()->size("PeakList") << "\"";
    document << "/>\n";
    if ( info_user_.length() || info_comment_.length() ) 
    {
      document << "<Info ";
      document << " User = \"" << info_user_ << "\"";
      time_t rawtime;
      struct tm * timeinfo;
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      info_date_ = String(asctime (timeinfo)).substr(0,String(asctime (timeinfo)).length() -1) ; 
      document << " Date = \"" << info_date_ << "\"";
      //delete timeinfo;
      if ( info_comment_.length() )
      {
        document << ">\n";
        document << indent(++ind);
        document << "<Comment> "<< indent(0)<<  info_comment_ << indent(ind) << "</Comment>\n";
        document << indent(--ind);
      }
      else 
      {
        document << "/>\n";
      }
    }
    for (uint i = 0; i < runps_.size(); ++i)
    {
      runps_[i]->save(document,ind);
    }
    document << "<AnalysisFunctorReferences>\n";
    document << indent(++ind);
    for ( uint ci = 0; ci < this->size(); ++ci )
    {
      for ( uint ai = 0; ai < this->runps_[ci]->size() ; ++ai )
      {
        AnalysisFunctor* afp = (*this)[ci][ai].anafuncp_;
        if ( afp->needsClusterRun() )
        {
          for ( uint ici = 0; ici < this->size(); ++ici )
          {
            if ( afp->reference() == &(*this)[ici] )
            {
              document << "<Reference ClusterRunNr=\"" << ci 
                << "\" AnalysisFunctorNr=\"" << ai 
                << "\" references=\"" << ici << "\" />\n";
            }
          }
        }
      }
    }
    document << indent(--ind);
    document << "</AnalysisFunctorReferences>\n";
    document << indent(--ind);
    document << "</ClusterExperiment>\n";
  }

  void ClusterExperiment::clear_()
  {
    //adapterp_ = 0;
    delete dsip_;
    dsip_ = 0;
    info_user_ = "";
    info_date_ = "";
    info_comment_ = "";
    runps_.clear();
    for ( uint i = 0; i < runps_.size(); ++i )
    {
      delete runps_[i];
    }
    currentrun_ = -1;
  }
 
  void ClusterExperiment::load(const String& filename) throw (Exception::FileNotFound, Exception::ParseError)
  {
  	//init
    clear_();

  	//try to open file
		if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		// initialize parser
		try 
		{
			xercesc::XMLPlatformUtils::Initialize();
		}
		catch (const xercesc::XMLException& toCatch) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	  }

		xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
		
		ClusterExperimentXMLHandler handler(*this,filename);
		
		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const xercesc::XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    catch (const xercesc::SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
  }

  const DataSetInfo* ClusterExperiment::getDataSetInfo() const
  {
    if ( !dsip_ )
    {
//TODO Persistence
//      if ( adapterp_ && datasetname_ != "" )
//      {
//        oStringstream ss;
//        ss << "SELECT id FROM DataSetInfo WHERE name = '" << datasetname_ << "'";
//        adapterp_->executeQuery(ss.str(),false);
//        QSqlQuery sqlres = adapterp_->lastResult();
//        
//        if ( sqlres.size() <=0 ) throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ClusterExperiment::DataSetInfo",datasetname_.c_str());
//        
//        //dsip_ = dynamic_cast<DataSetInfo*>(adapterp_->createObject(sqlres.value(0).toInt())); //TODO Persistence
//        if ( !dsip_ ) throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ClusterExperiment::DataSetInfo",datasetname_.c_str());
//      }
    }
    return dsip_;
  }
  
  void ClusterExperiment::setDBAdapter(DBAdapter* adapter)
  {
    adapterp_ = adapter;
  }

  void ClusterExperiment::setDataSetInfo(DataSetInfo* dsip)
  {
    dsip_ = dsip;
  }

  void ClusterExperiment::setInfoUser(String name)
  {
    info_user_ = name;
  }

  void ClusterExperiment::setInfoComment(String comment)
  {
    info_comment_ = comment;
  }

  uint ClusterExperiment::size() const
  {
    return runps_.size();
  }
 
  void ClusterExperiment::analyze(int pos)
  {
    datasetsize();
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->analyze();
    }
    else throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size() -1));
  }
  
  void ClusterExperiment::cluster(int pos)
  {
    datasetsize();
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->cluster();
    }
    else throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size() -1));
  }
  
  void ClusterExperiment::run(int pos)
  {
    datasetsize();
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->run();
    }
    else throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size() -1));
  }

  bool ClusterExperiment::canRun(int pos) const
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
     return runps_[pos]->canRun();
    }
    else return false;
  }

  void ClusterExperiment::setBinSize(double value ,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->setBinSize(value);
    } 
  }

  void ClusterExperiment::setNorm(Norm norm ,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->setNorm(norm);
    } 
  }

  void ClusterExperiment::setBinSpread(uint value,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->setBinSpread(value);
    } 
  }

  int ClusterExperiment::addMower(PreprocessingFunctor* funcp,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      return runps_[pos]->addMower(funcp);;
    } 
    else
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size());
    }
  }

  void ClusterExperiment::setSimFunc(CompareFunctor* funcp,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->setSimFunc(funcp);
    } 
  }

  void ClusterExperiment::setClusterFunc(ClusterFunctor* funcp,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      runps_[pos]->setClusterFunc(funcp);
    } 
  }

  int ClusterExperiment::addAnalysisFunctor(AnalysisFunctor* func,int pos)
  {
    if (pos == -1)
    {
      pos = currentrun_;
    }
    if ( pos < (int)runps_.size() )
    {
      return runps_[pos]->addAnalysisFunctor(func);
    } 
    else return -1;
  }
  
  /**
  new run is automatically current run (for write access) <br>
  */
  int ClusterExperiment::createrun()
  {
    runps_.push_back(new ClusterRun(*this));
    currentrun_ = runps_.size() -1;
    return currentrun_;
  }

  const ClusterExperiment::ClusterRun& ClusterExperiment::operator[] (uint pos) const throw(Exception::IndexOverflow)
  {
    if ( pos >= runps_.size() )
    {
      throw (Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size()));
    } 
    return *runps_[pos];
  }

  ClusterExperiment::ClusterRun* ClusterExperiment::getClusterRun (uint pos) throw(Exception::IndexOverflow)
  {
    if ( pos >= runps_.size() )
    {
      throw (Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,runps_.size()));
    } 
    return runps_[pos];
  }

  const String& ClusterExperiment::infodate() const
  {
    return info_date_;
  }
  
  const String& ClusterExperiment::infouser() const
  {
    return info_user_;
  }

  const String& ClusterExperiment::infocomment() const
  {
    return info_comment_;
  }
  
  long ClusterExperiment::datasetid() const
  {
    if ( getDataSetInfo() ) return getDataSetInfo()->getPersistenceId();
    else return 0;
  }
  
  String ClusterExperiment::datasetname() const
  {
    if ( getDataSetInfo() ) return getDataSetInfo()->name();
    else return "no DataSetInfo";
  }

  // sets the size and creates the dataset if needed
  // todo clean up ( exceptions are wrong etc )
  uint ClusterExperiment::datasetsize() 
  {
    if ( getDataSetInfo() ) return getDataSetInfo()->size("PeakList");
    else return 0;
  }

  void ClusterExperiment::sortbyResult(String name,String param, bool desc)
  {
    sort(runps_.begin(),runps_.end(), ClusterRunAnalysisLess(name,param));
    if (desc) reverse(runps_.begin(),runps_.end());
  }

}
