// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <set>
#include <string>
#include <iostream>
#include <stdexcept>
#include <list>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <boost/lexical_cast.hpp>

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {
    //TODO add documentation, warning and logger integration
    //TODO care for casts, switch validation on
    //TODO remodel CVTermList
    //TODO extend CVTermlist with CVCollection functionality for complete replacement??
    //TODO general id openms struct for overall parameter for one id run
    MzIdentMLDOMHandler::MzIdentMLDOMHandler(const std::vector<ProteinIdentification> & pro_id, const std::vector<PeptideIdentification> & pep_id, const String & version, const ProgressLogger & logger) :
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(0),
      pep_id_(0),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id),
      mzid_parser_()
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      try
      {
         XMLPlatformUtils::Initialize();  // Initialize Xerces infrastructure
      }
      catch( XMLException& e )
      {
         char* message = XMLString::transcode( e.getMessage() );
         cerr << "XML toolkit initialization error: " << message << endl;
         XMLString::release( &message );
         // throw exception here to return ERROR_XERCES_INIT
      }

      // Tags and attributes used in XML file.
      // Can't call transcode till after Xerces Initialize()
      TAG_root = XMLString::transcode("MzIdentML");
      TAG_CV = XMLString::transcode("cvParam");
      ATTR_name = XMLString::transcode("option_a");

    }

    MzIdentMLDOMHandler::MzIdentMLDOMHandler(std::vector<ProteinIdentification> & pro_id, std::vector<PeptideIdentification> & pep_id, const String & version, const ProgressLogger & logger) :
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(&pro_id),
      pep_id_(&pep_id),
      cpro_id_(0),
      cpep_id_(0),
      mzid_parser_()
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      try
      {
         XMLPlatformUtils::Initialize();  // Initialize Xerces infrastructure
      }
      catch( XMLException& e )
      {
         char* message = XMLString::transcode( e.getMessage() );
         cerr << "XML toolkit initialization error: " << message << endl;
         XMLString::release( &message );
         // throw exception here to return ERROR_XERCES_INIT
      }

      // Tags and attributes used in XML file.
      // Can't call transcode till after Xerces Initialize()
      TAG_root        = XMLString::transcode("MzIdentML");
      TAG_CV = XMLString::transcode("cvParam");
      ATTR_name = XMLString::transcode("name");

    }

    /**
     *  Class destructor frees memory used to hold the XML tag and
     *  attribute definitions. It als terminates use of the xerces-C
     *  framework.
     */
    MzIdentMLDOMHandler::~MzIdentMLDOMHandler()
    {
      //std::cout << "destroying dom tree" << std::endl;
      try
      {
         XMLString::release( &TAG_root );
         XMLString::release( &TAG_CV );
         XMLString::release( &ATTR_name );
//         if(m_name)   XMLString::release( &m_name ); //releasing you here is releasing you twice, dunno yet why?!
      }
      catch( ... )
      {
         cerr << "Unknown exception encountered in TagNamesdtor" << endl;
      }

      // Terminate Xerces
      //std::cout << "terminate xerces" << std::endl;
      try
      {
         XMLPlatformUtils::Terminate();  // Terminate after release of memory
      }
      catch( xercesc::XMLException& e )
      {
         char* message = xercesc::XMLString::transcode( e.getMessage() );

         cerr << "XML ttolkit teardown error: " << message << endl;
         XMLString::release( &message );
      }
    }

    /**
     *  This function:
     *  - reads mzid files
     *
     *  @param file.
     */
    void MzIdentMLDOMHandler::readMzIdentMLFile(const string& mzid_file)
            throw( std::runtime_error )
    {
       // Test to see if the file is ok.
       struct stat fileStatus;

       errno = 0;
       if(stat(mzid_file.c_str(), &fileStatus) == -1) // ==0 ok; ==-1 error
       {
           if( errno == ENOENT )      // errno declared by include file errno.h
              throw ( std::runtime_error("Path file_name does not exist, or path is an empty string.") );
           else if( errno == ENOTDIR )
              throw ( std::runtime_error("A component of the path is not a directory."));
           else if( errno == ELOOP )
              throw ( std::runtime_error("Too many symbolic links encountered while traversing the path."));
           else if( errno == EACCES )
              throw ( std::runtime_error("Permission denied."));
           else if( errno == ENAMETOOLONG )
              throw ( std::runtime_error("File can not be read\n"));
       }

       // Configure DOM parser.
       mzid_parser_.setValidationScheme( XercesDOMParser::Val_Never );
       mzid_parser_.setDoNamespaces( false );
       mzid_parser_.setDoSchema( false );
       mzid_parser_.setLoadExternalDTD( false );

       try
       {
          mzid_parser_.parse( mzid_file.c_str() );

          // no need to free this pointer - owned by the parent parser object
          DOMDocument* xmlDoc = mzid_parser_.getDocument();

          // 0. AnalysisSoftware {1,unbounded}
          DOMNodeList * analysisSoftwareElements = xmlDoc->getElementsByTagName(XMLString::transcode("AnalysisSoftware"));
          if( !analysisSoftwareElements ) throw(std::runtime_error( "No AnalysisSoftware nodes" ));
          parseAnalysisSoftwareList_(analysisSoftwareElements);

          // 1. AnalysisSampleCollection ??? contact stuff

          // 2. SpectrumIdentification  {1,unbounded} ! identification runs
          DOMNodeList * spectrumIdentificationElements = xmlDoc->getElementsByTagName(XMLString::transcode("SpectrumIdentification"));
          if( !spectrumIdentificationElements ) throw(std::runtime_error( "No SpectrumIdentification nodes" ));
          parseSpectrumIdentificationElements_(spectrumIdentificationElements);

          // 3. AnalysisProtocolCollection {1,1} SpectrumIdentificationProtocol  {1,unbounded} ! identification run parameters
          DOMNodeList * spectrumIdentificationProtocolElements = xmlDoc->getElementsByTagName(XMLString::transcode("SpectrumIdentificationProtocol"));
          if( !spectrumIdentificationProtocolElements ) throw(std::runtime_error( "No SpectrumIdentificationProtocol nodes" ));
          parseSpectrumIdentificationProtocolElements_(spectrumIdentificationProtocolElements);

          // 4. SequenceCollection nodes {0,1} DBSequenceElement {1,unbounded} Peptide {0,unbounded} PeptideEvidence {0,unbounded}
          DOMNodeList * dbSequenceElements = xmlDoc->getElementsByTagName(XMLString::transcode("DBSequence"));
          if( !dbSequenceElements ) throw(std::runtime_error( "No SequenceCollection/DBSequence nodes" ));
          parseDBSequenceElements_(dbSequenceElements);

          DOMNodeList * peptideElements = xmlDoc->getElementsByTagName(XMLString::transcode("Peptide"));
          if( !peptideElements ) throw(std::runtime_error( "No SequenceCollection/Peptide nodes" ));
          parsePeptideElements_(peptideElements);

          DOMNodeList * peptideEvidenceElements = xmlDoc->getElementsByTagName(XMLString::transcode("PeptideEvidence"));
          if( !peptideEvidenceElements ) throw(std::runtime_error( "No SequenceCollection/PeptideEvidence nodes" ));
          parsePeptideEvidenceElements_(peptideEvidenceElements);
//          mzid_parser_.resetDocumentPool(); //segfault prone: do not use!

          // 5. DataCollection {1,1}
          DOMNodeList * spectraDataElements = xmlDoc->getElementsByTagName(XMLString::transcode("SpectraData"));
          if( !spectraDataElements ) throw(std::runtime_error( "No SpectraData nodes" ));
          parseInputElements_(spectraDataElements);

          DOMNodeList * searchDatabaseElements = xmlDoc->getElementsByTagName(XMLString::transcode("SearchDatabase"));
          if( !searchDatabaseElements ) throw(std::runtime_error( "No SearchDatabase nodes" ));
          parseInputElements_(searchDatabaseElements);

          DOMNodeList * sourceFileElements = xmlDoc->getElementsByTagName(XMLString::transcode("SourceFile"));
          if( !sourceFileElements ) throw(std::runtime_error( "No SourceFile nodes" ));
          parseInputElements_(sourceFileElements);

          // 6. AnalysisCollection {1,1} - build final structures
          DOMNodeList * spectrumIdentificationResultElements = xmlDoc->getElementsByTagName(XMLString::transcode("SpectrumIdentificationResult"));
          if( !spectrumIdentificationResultElements ) throw(std::runtime_error( "No SpectrumIdentificationItem nodes" ));
          parseSpectrumIdentificationResultElements_(spectrumIdentificationResultElements);

          DOMNodeList * parseProteinDetectionListElements = xmlDoc->getElementsByTagName(XMLString::transcode("ProteinDetectionList"));
          if( !parseProteinDetectionListElements ) throw(std::runtime_error( "No ProteinDetectionList nodes" ));
          parseProteinDetectionListElements_(parseProteinDetectionListElements);

          //std::cout << "I'm done. " << std::endl << "You have " << pep_id_->size() << " PeptideIdentifications and " << pro_id_->size() << " ProteinIdentifications" << std::endl; //wtf?! 8458800 PeptideIdentifications in idxml vs. Have 5936PeptideIdentifications and 1425ProteinIdentifications

       }
       catch( xercesc::XMLException& e )
       {
          char* message = xercesc::XMLString::transcode( e.getMessage() );
          ostringstream errBuf;
          errBuf << "Error parsing file: " << message << flush;
          XMLString::release( &message );
       }
    }

    void MzIdentMLDOMHandler::writeMzIdentMLFile(const string& mzid_file)
            throw( std::runtime_error )
    {
       //TODO
    }

    std::pair<CVTermList, std::map<String,DataValue> > MzIdentMLDOMHandler::parseParamGroup_(DOMNodeList * paramGroup)
    {
      CVTermList ret_cv;
      std::map<String,DataValue> ret_up;
      const  XMLSize_t cv_node_count = paramGroup->getLength();
      for( XMLSize_t cvi = 0; cvi < cv_node_count; ++cvi )
      {
        DOMNode* current_cv = paramGroup->item(cvi);
        if( current_cv->getNodeType() &&  // true is not NULL
            current_cv->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          DOMElement* element_param = dynamic_cast< xercesc::DOMElement* >( current_cv );
          if ((std::string)XMLString::transcode(element_param->getTagName()) == "cvParam")
          {
            ret_cv.addCVTerm(parseCvParam_(element_param));
          }
          else if ((std::string)XMLString::transcode(element_param->getTagName()) == "userParam")
          {
            ret_up.insert(parseUserParam_(element_param));
          }
          else
          {}// something else
        }
      }
      return std::make_pair(ret_cv,ret_up);
    }

    CVTerm MzIdentMLDOMHandler::parseCvParam_(DOMElement* param)
    {
      if ( param )
      {
        //      <cvParam accession="MS:1001469" name="taxonomy: scientific name" cvRef="PSI-MS"  value="Drosophila melanogaster"/>
        String accession = XMLString::transcode(param->getAttribute(XMLString::transcode("accession")));
        String name = XMLString::transcode(param->getAttribute(XMLString::transcode("name")));
        String cvRef = XMLString::transcode(param->getAttribute(XMLString::transcode("cvRef")));
        String value = XMLString::transcode(param->getAttribute(XMLString::transcode("value")));

        String unitAcc = XMLString::transcode(param->getAttribute(XMLString::transcode("unitAccession")));
        String unitName = XMLString::transcode(param->getAttribute(XMLString::transcode("unitName")));
        String unitCvRef = XMLString::transcode(param->getAttribute(XMLString::transcode("unitCvRef")));

        CVTerm::Unit u; // TODO @mths : make DataValue usage safe!
        if (!unitAcc.empty() && unitCvRef.empty() && unitName.empty())
        {
          u = CVTerm::Unit(unitAcc,unitCvRef,unitName);
        }
        return CVTerm(accession, name, cvRef, value, u);
      }
      else
        throw "no cv here.";
    }

    std::pair<String,DataValue> MzIdentMLDOMHandler::parseUserParam_(DOMElement* param)
    {
      if ( param )
      {
        //      <userParam name="Mascot User Comment" value="Example Mascot MS-MS search for PSI mzIdentML"/>
        String name = XMLString::transcode(param->getAttribute(XMLString::transcode("name")));
        String value = XMLString::transcode(param->getAttribute(XMLString::transcode("value")));

        String unitAcc = XMLString::transcode(param->getAttribute(XMLString::transcode("unitAccession")));
        String unitName = XMLString::transcode(param->getAttribute(XMLString::transcode("unitName")));
        String unitCvRef = XMLString::transcode(param->getAttribute(XMLString::transcode("unitCvRef")));
        String type = XMLString::transcode(param->getAttribute(XMLString::transcode("type")));

        DataValue dv;
        dv.setUnit(unitAcc+":"+unitName);

        if ( type == "xsd:float" || type == "xsd:double")
        {
          dv = value.toDouble();
        }
        else if ( type == "xsd:int" || type == "xsd:unsignedInt")
        {
          dv = value.toInt();
        }
        else
        {
          dv = value;
        }

        return std::make_pair(name,dv);
      }
      else
        throw "no up here.";
    }

    void MzIdentMLDOMHandler::parseAnalysisSoftwareList_(DOMNodeList * analysisSoftwareElements)
    {
      const  XMLSize_t as_node_count = analysisSoftwareElements->getLength();
      for( XMLSize_t swni = 0; swni < as_node_count; ++swni )
      {
        DOMNode* current_as = analysisSoftwareElements->item(swni);
        if( current_as->getNodeType() &&  // true is not NULL
            current_as->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_AnalysisSoftware = dynamic_cast< xercesc::DOMElement* >( current_as );
          //std::cout << "analysis software found: " << element_AnalysisSoftware->getAttribute(XMLString::transcode("version")) << std::endl;

          DOMElement* child = element_AnalysisSoftware->getFirstElementChild();
          while ( child )
          {
            //std::cout << "as child: " << XMLString::transcode(child->getTagName()) << std::endl;
            if ((std::string)XMLString::transcode(child->getTagName()) == "SoftwareName")
            {
              DOMElement* element_cv = child->getFirstElementChild();
              if( element_cv )
              {
                  if ((std::string)XMLString::transcode(element_cv->getTagName()) == "cvParam") // cave: might also be a UserParam
                  {
                    CVTerm swcv = parseCvParam_(element_cv);
                    //std::cout << "as child cv: " << swcv.getName() << std::endl;
                    if (search_engine_ == "" && search_engine_version_ == "") // TODO @mths check if cv is search engine cv!
                    {
                      search_engine_ = swcv.getName();
                      search_engine_version_ = XMLString::transcode(element_AnalysisSoftware->getAttribute(XMLString::transcode("version")));
                      //std::cout << "engine found: " << search_engine_ << std::endl;
                      //std::cout << "version found: " << search_engine_version_ << std::endl;
                    }
                    //else what?! refactor software in OpenMS!!
                  }
              }
            }
            child = child->getNextElementSibling();
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseDBSequenceElements_(DOMNodeList * dbSequenceElements)
    {
      const  XMLSize_t dbs_node_count = dbSequenceElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < dbs_node_count; ++c )
      {
        DOMNode * current_dbs = dbSequenceElements->item(c);
        if( current_dbs->getNodeType() &&  // true is not NULL
            current_dbs->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_dbs = dynamic_cast< xercesc::DOMElement* >( current_dbs );
          ++count;
          String id = XMLString::transcode(element_dbs->getAttribute(XMLString::transcode("id")));
          String seq = "";
          String dbref = XMLString::transcode(element_dbs->getAttribute(XMLString::transcode("searchDatabase_ref")));
          String acc = XMLString::transcode(element_dbs->getAttribute(XMLString::transcode("accession")));
          CVTermList cvs;

          DOMElement* child = element_dbs->getFirstElementChild();
          while ( child )
          {
            //std::cout << "DBSequences child" << std::endl;
            if ((std::string)XMLString::transcode(child->getTagName()) == "Seq")
            {
              seq = XMLString::transcode(child->getNodeValue()); //http://mail-archives.apache.org/mod_mbox/xerces-c-users/200610.mbox/%3C7.0.1.0.2.20061004133628.046911a8@datadirect.com%3E -> dynamic_cast< xercesc::DOMNode* >
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "cvParam")
            {
              cvs.addCVTerm( parseCvParam_(child) );
            }
            child = child->getNextElementSibling();
          }
          if (acc != "")
          {
            db_sq_map_.insert(std::make_pair(id, DBSequence{seq,dbref,acc,cvs}));
          }
        }
      }
      //std::cout << "DBSequences found: " << count  << " / " << db_sq_map_.size() << std::endl;
    }

    void MzIdentMLDOMHandler::parsePeptideElements_(DOMNodeList * peptideElements)
    {
      const  XMLSize_t pep_node_count = peptideElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < pep_node_count; ++c )
      {
        DOMNode * current_pep = peptideElements->item(c);
        if( current_pep->getNodeType() &&  // true is not NULL
            current_pep->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pep = dynamic_cast< xercesc::DOMElement* >( current_pep );
          ++count;
          String id = XMLString::transcode(element_pep->getAttribute(XMLString::transcode("id")));


          DOMNodeList* pep_sib = element_pep->getChildNodes();
          AASequence aas = parsePeptideSiblings_(pep_sib);

          pep_map_.insert(std::make_pair(id,aas));
        }
      }
      //std::cout << "Peptides found: " << count << std::endl;
      //std::cout << "example: " << pep_map_.begin()->first << " -> " << pep_map_.begin()->second.toString() << std::endl;
    }

    void MzIdentMLDOMHandler::parsePeptideEvidenceElements_(DOMNodeList * peptideEvidenceElements)
    {
      const  XMLSize_t pev_node_count = peptideEvidenceElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < pev_node_count; ++c )
      {
        DOMNode * current_pev = peptideEvidenceElements->item(c);
        if( current_pev->getNodeType() &&  // true is not NULL
            current_pev->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pev = dynamic_cast< xercesc::DOMElement* >( current_pev );
          ++count;

//          <PeptideEvidence peptide_ref="peptide_1_1" id="PE_1_1_HSP70_ECHGR_0" start="161" end="172" pre="K" post="I" isDecoy="false" dBSequence_ref="DBSeq_HSP70_ECHGR"/>

          String id = XMLString::transcode(element_pev->getAttribute(XMLString::transcode("id")));
          String peptide_ref = XMLString::transcode(element_pev->getAttribute(XMLString::transcode("peptide_ref")));
          String dBSequence_ref = XMLString::transcode(element_pev->getAttribute(XMLString::transcode("dBSequence_ref")));
          //rest is optional !!
          int start = String(XMLString::transcode(element_pev->getAttribute(XMLString::transcode("start")))).toInt();
          int end = String(XMLString::transcode(element_pev->getAttribute(XMLString::transcode("end")))).toInt();
          char pre = *XMLString::transcode(element_pev->getAttribute(XMLString::transcode("pre")));
          char post = *XMLString::transcode(element_pev->getAttribute(XMLString::transcode("post")));

          pe_ev_map_.insert(std::make_pair(id,PeptideEvidence{start,end,pre,post}));
          p_pv_map_.insert(std::make_pair(peptide_ref,id));
          pv_db_map_.insert(std::make_pair(id,dBSequence_ref));

        }
      }
      //std::cout << "PeptideEvidences found: " << count << std::endl;
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationElements_(DOMNodeList * spectrumIdentificationElements)
    {
      const  XMLSize_t si_node_count = spectrumIdentificationElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < si_node_count; ++c )
      {
        DOMNode * current_si = spectrumIdentificationElements->item(c);
        if( current_si->getNodeType() &&  // true is not NULL
            current_si->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_si = dynamic_cast< xercesc::DOMElement* >( current_si );
          ++count;
          String id = XMLString::transcode(element_si->getAttribute(XMLString::transcode("id")));
          String spectrumIdentificationProtocol_ref = XMLString::transcode(element_si->getAttribute(XMLString::transcode("spectrumIdentificationProtocol_ref")));
          String spectrumIdentificationList_ref = XMLString::transcode(element_si->getAttribute(XMLString::transcode("spectrumIdentificationList_ref")));

          String searchDatabase_ref = "";
          String spectra_data_ref = "";
          DOMElement* child = element_si->getFirstElementChild();
          while ( child )
          {
            if ((std::string)XMLString::transcode(child->getTagName()) == "InputSpectra")
            {
              spectra_data_ref = XMLString::transcode(element_si->getAttribute(XMLString::transcode("spectraData_ref")));
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "SearchDatabaseRef")
            {
              searchDatabase_ref = XMLString::transcode(element_si->getAttribute(XMLString::transcode("searchDatabase_ref")));
            }
            child = child->getNextElementSibling();
          }

          si_map_.insert(std::make_pair(id,SpectrumIdentification{spectra_data_ref, searchDatabase_ref, spectrumIdentificationProtocol_ref, spectrumIdentificationList_ref}));

          pro_id_->push_back(ProteinIdentification());
          si_pro_map_.insert(std::make_pair(spectrumIdentificationList_ref,&pro_id_->back()));

        }
      }
      //std::cout << "SpectrumIdentification found: " << count << std::endl;
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationProtocolElements_(DOMNodeList * spectrumIdentificationProtocolElements)
    {
      const  XMLSize_t si_node_count = spectrumIdentificationProtocolElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < si_node_count; ++c )
      {
        DOMNode * current_sip = spectrumIdentificationProtocolElements->item(c);
        if( current_sip->getNodeType() &&  // true is not NULL
            current_sip->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_sip = dynamic_cast< xercesc::DOMElement* >( current_sip );
          ++count;
          String id = XMLString::transcode(element_sip->getAttribute(XMLString::transcode("id")));

          CVTerm searchtype;
          String enzyme;
          CVTermList param_cv;
          std::map<String, DataValue> param_up;
          CVTermList modparam;
          long double p_tol = 0;
          long double f_tol = 0;
          CVTermList tcv;
          std::map<String, DataValue> tup;

          DOMElement* child = element_sip->getFirstElementChild();
          while ( child )
          {
            if ((std::string)XMLString::transcode(child->getTagName()) == "SearchType")
            {
              searchtype = parseCvParam_(child->getFirstElementChild());
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "AdditionalSearchParams")
            {
              std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(child->getChildNodes());

            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "ModificationParams") // TODO @all where to store the specificities?
            {
              DOMElement* sm = child->getFirstElementChild();
              while ( sm )
              {
                  // TODO @mths: check availability before cast
//                String fixedMod = XMLString::transcode(child->getAttribute(XMLString::transcode("fixedMod")));
//                long double massDelta = String(XMLString::transcode(child->getAttribute(XMLString::transcode("massDelta")))).toDouble();
//                String residues = XMLString::transcode(child->getAttribute(XMLString::transcode("residues")));
                CVTermList specificities;

                DOMElement* rule = sm->getFirstElementChild();
                while ( rule )
                {
                  specificities.consumeCVTerms(parseParamGroup_(rule->getChildNodes()).first.getCVTerms());
                  rule = rule->getNextElementSibling();
                }
                sm = sm->getNextElementSibling();
              }
              modparam = parseParamGroup_(child->getChildNodes()).first;
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "Enzymes") // TODO @all : where store multiple enzymes for one identificationrun?
            {
              if ((std::string)XMLString::transcode(child->getAttribute(XMLString::transcode("independent"))) == "0")
              {
                DOMElement* enzy = child->getFirstElementChild();
                while ( enzy )
                {
                  String id = XMLString::transcode(enzy->getAttribute(XMLString::transcode("id")));
                  String name = XMLString::transcode(enzy->getAttribute(XMLString::transcode("name")));
                  int missedCleavages = String(XMLString::transcode(enzy->getAttribute(XMLString::transcode("missedCleavages")))).toInt();
                  String semiSpecific = XMLString::transcode(enzy->getAttribute(XMLString::transcode("semiSpecific"))); //xsd:boolean
                  String cTermGain = XMLString::transcode(enzy->getAttribute(XMLString::transcode("cTermGain")));
                  String nTermGain = XMLString::transcode(enzy->getAttribute(XMLString::transcode("nTermGain")));
                  int minDistance = String(XMLString::transcode(enzy->getAttribute(XMLString::transcode("minDistance")))).toInt();
                  String enzymename = "UNKNOWN";

                  DOMElement* sub = enzy->getFirstElementChild();
                  if ( sub )
                  {
                    //SiteRegex unstorable just now
                    if ((std::string)XMLString::transcode(sub->getTagName()) != "EnzymeName")
                      sub = sub->getNextElementSibling();
                    if ( sub )
                    {
                      //take the first param for name
                      DOMElement* pren = sub->getFirstElementChild();
                      if ((std::string)XMLString::transcode(pren->getTagName()) == "userParam")
                      {
                        CVTerm param = parseCvParam_(pren->getFirstElementChild());
                        enzymename = param.getValue();
                      }
                      else if ((std::string)XMLString::transcode(pren->getTagName()) == "cvParam")
                      {
                        std::pair<String,DataValue> param = parseUserParam_(pren->getFirstElementChild());
                        enzymename = param.second.toString();
                      }
                    }
                  }
                  enzyme = enzymename;
                  enzy = enzy->getNextElementSibling();
                }
              } // else uhoh?! what nao?
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "FragmentTolerance")
            {
              std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(child->getChildNodes());
              //+- take the numerically greater
              for (std::map<String,std::vector<CVTerm> >::const_iterator it=params.first.getCVTerms().begin(); it!=params.first.getCVTerms().end(); ++it)
              {
                f_tol = std::max(f_tol, (long double)it->second.front().getValue());
              }
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "ParentTolerance")
            {
              std::pair<CVTermList,std::map<String,DataValue> > params= parseParamGroup_(child->getChildNodes());
              //+- take the numerically greater
              for (std::map<String,std::vector<CVTerm> >::const_iterator it=params.first.getCVTerms().begin(); it!=params.first.getCVTerms().end(); ++it)
              {
                p_tol = std::max(p_tol, (long double)it->second.front().getValue());
              }
            }
            else if ((std::string)XMLString::transcode(child->getTagName()) == "Threshold")
            {
              std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(child->getChildNodes());
              tcv = params.first;
              tup = params.second;
            }
            child = child->getNextElementSibling();
            //      <DatabaseFilters> omitted for now, not reflectable by our member structures
            //      <DatabaseTranslation> omitted for now, not reflectable by our member structures
            //      <Masstable> omitted for now, not reflectable by our member structures
          }
          sp_map_.insert(std::make_pair(id,SpectrumIdentificationProtocol{searchtype,enzyme,param_cv,param_up,modparam,p_tol,f_tol,tcv,tup}));

          for (Map< String, SpectrumIdentification >::ConstIterator si_it = si_map_.begin(); si_it != si_map_.end(); ++si_it)
          {
            if (si_it->second.spectrum_identification_protocol_ref == id)
            {
              // TODO @mths get these from AnalysisSoftware!
              si_pro_map_[si_it->first]->setSearchEngine(search_engine_);
              si_pro_map_[si_it->first]->setSearchEngineVersion(search_engine_version_);
              si_pro_map_[si_it->first]->setIdentifier(search_engine_); // TODO @mths: name/date of search
              // TODO @mths set SearchParameters 	search_parameters_
            }
          }
        }
      }
      //std::cout << "SpectrumIdentificationProtocol found: " << count << std::endl;
    }

    void MzIdentMLDOMHandler::parseInputElements_(DOMNodeList * inputElements)
    {
      const  XMLSize_t node_count = inputElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
        DOMNode * current_in = inputElements->item(c);
        if( current_in->getNodeType() &&  // true is not NULL
            current_in->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_in = dynamic_cast< xercesc::DOMElement* >( current_in );
          ++count;

          String id = XMLString::transcode(element_in->getAttribute(XMLString::transcode("id")));
          String location = XMLString::transcode(element_in->getAttribute(XMLString::transcode("location")));

          if ((std::string)XMLString::transcode(element_in->getTagName()) == "SpectraData")
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            //      <SpectrumIDFormat> omitted for now, not reflectable by our member structures
            input_spectra_data_.insert(std::make_pair(id,location));
          }
          else if ((std::string)XMLString::transcode(element_in->getTagName()) == "SourceFile")
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            input_source_.insert(std::make_pair(id,location));
          }
          else if ((std::string)XMLString::transcode(element_in->getTagName()) == "SearchDatabase")
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            DateTime releaseDate;
//            TODO check befor cast
//            releaseDate.set(String(XMLString::transcode(element_in->getAttribute(XMLString::transcode("releaseDate")))));
            String version = XMLString::transcode(element_in->getAttribute(XMLString::transcode("version")));
            //assumed that <DatabaseName> is the first child, following cv omitted for now
            String dbname = "";
            DOMElement* pren = element_in->getFirstElementChild();
            if ((std::string)XMLString::transcode(pren->getTagName()) == "userParam")
            {
              CVTerm param = parseCvParam_(pren->getFirstElementChild());
              dbname = param.getValue();
            }
            else if ((std::string)XMLString::transcode(pren->getTagName()) == "cvParam")
            {
              std::pair<String,DataValue> param = parseUserParam_(pren->getFirstElementChild());
              dbname = param.second.toString();
            }
            input_dbs_.insert(std::make_pair(id, DatabaseInput{dbname,location,version,releaseDate}));
          }
        }
        //std::cout << "InputFiles found: " << count << std::endl;
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationResultElements_(DOMNodeList * spectrumIdentificationResultElements)
    {
      const  XMLSize_t node_count = spectrumIdentificationResultElements->getLength();
      int count = 0;
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
        DOMNode * current_res = spectrumIdentificationResultElements->item(c);
        if( current_res->getNodeType() &&  // true is not NULL
            current_res->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_res = dynamic_cast< xercesc::DOMElement* >( current_res );
          ++count;

          String id = XMLString::transcode(element_res->getAttribute(XMLString::transcode("id")));
          String name = XMLString::transcode(element_res->getAttribute(XMLString::transcode("name")));
          String spectra_data_ref = XMLString::transcode(element_res->getAttribute(XMLString::transcode("spectraData_ref")));
          String spectrumID = XMLString::transcode(element_res->getAttribute(XMLString::transcode("spectrumID")));
          std::vector< String > substrings;
          spectrumID.split(" ", substrings);
          String scannr = substrings.back().reverse().chop(5);

          std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(element_res->getChildNodes());

          pep_id_->push_back(PeptideIdentification());
          //fill pep_id_->back() with content

          //butt ugly!
          DOMElement* parent = dynamic_cast< xercesc::DOMElement* >( element_res->getParentNode() );
          String sil = XMLString::transcode(parent->getAttribute(XMLString::transcode("id")));

          DOMElement* child = element_res->getFirstElementChild();
          while ( child )
          {
            if ((std::string)XMLString::transcode(child->getTagName()) == "SpectrumIdentificationItem")
            {
              parseSpectrumIdentificationItemElement_(child, pep_id_->back(),sil);
            }
            child = child->getNextElementSibling();
          }

//          setSignificanceThreshold
//          Origin stored in the MetaInfoInterface using the names 'MZ' and 'RT'
//          'real' origin in MetaInfoInterface 'spectrumID' as NativeID
          pep_id_->back().setIdentifier(search_engine_); // TODO @mths: set name/date of search
          pep_id_->back().setMetaValue("scannumber", scannr);

        }
      }
      //std::cout << "SpectrumIdentificationResults found: " << count << std::endl;
      //std::cout << "example: " << pep_id_->back().getHits().size() << std::endl;
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationItemElement_(DOMElement * spectrumIdentificationItemElement, PeptideIdentification& spectrum_identification, String& spectrumIdentificationList_ref)
    {
      String id = XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("id")));
      String name = XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("name")));

      // TODO @ mths : where to put calc. mz if even
//      long double calculatedMassToCharge = String(XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("calculatedMassToCharge")))).toDouble();
//      long double calculatedPI = String(XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("calculatedPI")))).toDouble();
      int chargeState = String(XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("chargeState")))).toInt();
      long double experimentalMassToCharge = String(XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("experimentalMassToCharge")))).toDouble();
      int rank = String(XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("rank")))).toInt();
      String peptide_ref = XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("peptide_ref")));
//      String sample_ref = XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("sample_ref")));
//      String massTable_ref = XMLString::transcode(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("massTable_ref")));

      XSValue::Status status;
      XSValue *val = XSValue::getActualValue(spectrumIdentificationItemElement->getAttribute(XMLString::transcode("passThreshold")), XSValue::dt_boolean, status);
      bool pass = false;
      if ( status == XSValue::st_Init )
        pass = val->fData.fValue.f_bool;
      delete val;

      long double score = 0;
      std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(spectrumIdentificationItemElement->getChildNodes());
      if (params.first.getCVTerms().has("MS:1001172"))// TODO @mths get children of MS:1001143
      {
        score = (long double)params.first.getCVTerms()["MS:1001172"].front().getValue();
      }

      String& pev = p_pv_map_[peptide_ref];
      String& dpv = pv_db_map_[pev];

      PeptideHit hit(score, rank, chargeState, pep_map_[peptide_ref]);
      for (Map< String, std::vector< CVTerm > >::ConstIterator cvs = params.first.getCVTerms().begin(); cvs != params.first.getCVTerms().end(); ++cvs)
      {
          for (std::vector< CVTerm >::const_iterator cv = cvs->second.begin(); cv != cvs->second.end(); ++cv)
          {
            hit.setMetaValue(cvs->first, cv->getValue());
          }
      }
      for (std::map<String,DataValue>::const_iterator up = params.second.begin(); up != params.second.end(); ++up)
      {
        hit.setMetaValue(up->first, up->second);
      }
      spectrum_identification.insertHit(hit);
      spectrum_identification.setMZ(experimentalMassToCharge); // TODO @ mths: why is this not in SpectrumIdentificationResult? exp. m/z for one spec should not change from one id for it to the next!

      if (pe_ev_map_.find(pev) != pe_ev_map_.end())
      {
        PeptideEvidence& pv = pe_ev_map_[pev];
        spectrum_identification.getHits().back().setAABefore(pv.pre);
        spectrum_identification.getHits().back().setAAAfter(pv.post);
      }

      DBSequence& db = db_sq_map_[dpv];
      spectrum_identification.getHits().back().addProteinAccession(db.accession);    

      if (si_pro_map_[spectrumIdentificationList_ref]->findHit(db.accession)
          != si_pro_map_[spectrumIdentificationList_ref]->getHits().end())
      { // butt ugly!
        si_pro_map_[spectrumIdentificationList_ref]->insertHit(ProteinHit());
        si_pro_map_[spectrumIdentificationList_ref]->getHits().back().setSequence(db.sequence);
        si_pro_map_[spectrumIdentificationList_ref]->getHits().back().setAccession(db.accession);
      }

      //due to redundand references this is not needed!
//      DOMElement* child = spectrumIdentificationItemElement->getFirstElementChild();
//      while ( child )
//      {
//        if ((std::string)XMLString::transcode(child->getTagName()) == "PeptideEvidenceRef")
//        {
//          ref = XMLString::transcode(element_si->getAttribute(XMLString::transcode("peptideEvidence_ref")));
//          //...
//          spectrum_identification.getHits().back().setAABefore(char acid);
//          spectrum_identification.getHits().back().setAAAfter (char acid);
//          break;
//        }
//        child = child->getNextElementSibling();
//      }

//        <Fragmentation> omitted for the time being

    }

    void MzIdentMLDOMHandler::parseProteinDetectionListElements_(DOMNodeList * proteinDetectionListElements)
    {
      const  XMLSize_t node_count = proteinDetectionListElements->getLength();
      int count = 0;int count_ag = 0;
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
        DOMNode * current_pr = proteinDetectionListElements->item(c);
        if( current_pr->getNodeType() &&  // true is not NULL
            current_pr->getNodeType() == DOMNode::ELEMENT_NODE ) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pr = dynamic_cast< xercesc::DOMElement* >( current_pr );
          ++count;

          String id = XMLString::transcode(element_pr->getAttribute(XMLString::transcode("id")));
          std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(current_pr->getChildNodes());

          // TODO @mths : this needs to be a ProteinIdentification for the ProteinDetectionListElement which is not mandatory and used in downstream analysis ProteinInference etc.
          pro_id_->push_back(ProteinIdentification());
          pro_id_->back().setSearchEngine(search_engine_);
          pro_id_->back().setSearchEngineVersion(search_engine_version_);
          pro_id_->back().setIdentifier(search_engine_); // TODO @mths: name/date of search

  //      SearchParameters 	search_parameters_
  //      DateTime 	date_
  //      String 	protein_score_type_ <- from proteindetectionprotocol
  //      DoubleReal 	protein_significance_threshold_ <- from proteindetectionprotocol

          DOMElement* child = element_pr->getFirstElementChild();
          while ( child )
          {
            if ((std::string)XMLString::transcode(child->getTagName()) == "ProteinAmbiguityGroup")
            {
              parseProteinAmbiguityGroupElement_(child, pro_id_->back());
            }
            child = child->getNextElementSibling();
            ++count_ag;
          }

        }
        //std::cout << "ProteinDetectionLists found: " << count << std::endl;
        //std::cout << "ProteinAmbiguityGroups found: " << count_ag << std::endl;
      }
    }

    void MzIdentMLDOMHandler::parseProteinAmbiguityGroupElement_(DOMElement * proteinAmbiguityGroupElement, ProteinIdentification& protein_identification)
    {
       String id = XMLString::transcode(proteinAmbiguityGroupElement->getAttribute(XMLString::transcode("id")));
       std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(proteinAmbiguityGroupElement->getChildNodes());

       //fill pro_id_->back() with content,
       DOMElement* child = proteinAmbiguityGroupElement->getFirstElementChild();
       while ( child )
       {
         if ((std::string)XMLString::transcode(child->getTagName()) == "ProteinDetectionHypothesis")
         {
           parseProteinDetectionHypothesisElement_(child, protein_identification);
         }
         child = child->getNextElementSibling();
       }
    }

    void MzIdentMLDOMHandler::parseProteinDetectionHypothesisElement_(DOMElement * proteinDetectionHypothesisElement, ProteinIdentification& protein_identification)
    {
      String dBSequence_ref = XMLString::transcode(proteinDetectionHypothesisElement->getAttribute(XMLString::transcode("dBSequence_ref")));

      std::pair<CVTermList,std::map<String,DataValue> > params = parseParamGroup_(proteinDetectionHypothesisElement->getChildNodes());

      DBSequence& db = db_sq_map_[dBSequence_ref];

      protein_identification.insertHit(ProteinHit());
      protein_identification.getHits().back().setSequence(db.sequence);
      protein_identification.getHits().back().setAccession(db.accession);
//      protein_identification.getHits().back().setCoverage((long double)params.first.getCVTerms()["MS:1001093"].front().getValue()); //TODO @ mths: calc percent
      protein_identification.getHits().back().setScore((long double)params.first.getCVTerms()["MS:1001171"].front().getValue()); //or any other score

    }

    AASequence MzIdentMLDOMHandler::parsePeptideSiblings_(DOMNodeList * peptideSiblings)
    {
      const  XMLSize_t node_count = peptideSiblings->getLength();
      String as;
      //1. Sequence
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
        DOMNode * current_sib = peptideSiblings->item(c);
        if( current_sib->getNodeType() &&  // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE )
        {
          DOMElement* element_sib = dynamic_cast< xercesc::DOMElement* >( current_sib );
          if ((std::string)XMLString::transcode(element_sib->getTagName()) == "PeptideSequence")
          {
            DOMNode* tn = element_sib->getFirstChild();
            if(tn->getNodeType() == DOMNode::TEXT_NODE)
            {
              DOMText* data = dynamic_cast<DOMText*>(tn);
              const XMLCh* val = data->getWholeText();
              as = String(XMLString::transcode(val));
            }
            else
            {
              //std::cout << "going to throw up from <PeptideSequence>" << std::endl;
              throw "ERROR : Non Text Node";
            }
          }
        }
      }
      //2. Substitutions
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
//        std::cout << "as: " << as;
        DOMNode * current_sib = peptideSiblings->item(c);
        if( current_sib->getNodeType() &&  // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE )
        {
          DOMElement* element_sib = dynamic_cast< xercesc::DOMElement* >( current_sib );
          if ((std::string)XMLString::transcode(element_sib->getTagName()) == "SubstitutionModification")
          {

            String location = XMLString::transcode(element_sib->getAttribute(XMLString::transcode("location")));
            char originalResidue = std::string(XMLString::transcode(element_sib->getAttribute(XMLString::transcode("originalResidue"))))[0];
            char replacementResidue = std::string(XMLString::transcode(element_sib->getAttribute(XMLString::transcode("replacementResidue"))))[0];

            if (!location.empty())
            {
              as[location.toInt()-1] = replacementResidue;
//              std::cout << as[location.toInt()-1] << ": " << originalResidue << "->" << replacementResidue << std::endl;
            }
            else if (as.hasSubstring(originalResidue)) //no location - every occurrence will be replaced
            {
              as.substitute(originalResidue, replacementResidue);
//              std::cout << originalResidue << "->" << replacementResidue << std::endl;
            }
            else
            {
              std::cout << "going to throw up from <PeptideSequence>" << std::endl;
              throw "ERROR : Non Text Node";
            }
          }
        }
//        std::cout << " as_subst: " << as << std::endl;
      }
      //3. Modifications
      AASequence aas = AASequence::fromString(as);
      for( XMLSize_t c = 0; c < node_count; ++c )
      {
        DOMNode * current_sib = peptideSiblings->item(c);
        if( current_sib->getNodeType() &&  // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE )
        {
          DOMElement* element_sib = dynamic_cast< xercesc::DOMElement* >( current_sib );
          if ((std::string)XMLString::transcode(element_sib->getTagName()) == "Modification")
          {
            // all attributes optional :( !!!
            int index = String(XMLString::transcode(element_sib->getAttribute(XMLString::transcode("location")))).toInt();
            //double monoisotopicMassDelta = XMLString::transcode(element_dbs->getAttribute(XMLString::transcode("monoisotopicMassDelta")));
//            std::cout << "index: " << index << std::endl;
            DOMElement* cvp = element_sib->getFirstElementChild();
            while ( cvp )
            {
              CVTerm cv = parseCvParam_(cvp);
              if (cv.getCVIdentifierRef() != "UNIMOD")
              {
    //                 e.g.  <cvParam accession="MS:1001524" name="fragment neutral loss" cvRef="PSI-MS" value="0" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO"/>
                 cvp = cvp->getNextElementSibling();
                 continue;
              }
              if (index == 0)
                aas.setNTerminalModification(cv.getName());
              else if ( index == aas.size()+1 )
                aas.setCTerminalModification(cv.getName());
              else
              {
                try
                {
                  aas.setModification(index-1, cv.getName()); //TODO @mths,Timo : do this via UNIMOD accessions
                }
                catch (int e)
                {
                    std::cout << " res: " << aas.getResidue((SignedSize)index).getName() << String(index) << " mod: " << cv.getName() << std::endl;
                    std::cout << "!!!as: " << as << "aas: " << aas.toString() <<  " res: " << aas.toString()[index] << String(index) << " mod: " << cv.getName() << std::endl;
                }
              }
              cvp = cvp->getNextElementSibling();
            }
          }
        }
      }
      return aas;
    }


  }   //namespace Internal
} // namespace OpenMS
