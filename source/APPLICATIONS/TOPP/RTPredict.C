// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page RTPredict RTPredict
	
	@brief This application is used to predict retention times for peptides
	
	The input of this application is a svm model and an analysisXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the RTModel application. 
	<br>
	The peptide sequences are extracted from the analysisXML inputfile 
	and passed to the svm. The svm then predicts retention times
	according to the trained model. The predicted retention times
	are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" /> 
	@endcode inside the peptide entities in the analysisXML output file.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTPredict
	: public TOPPBase
{
	public:
		TOPPRTPredict()
			: TOPPBase("RTPredict","predicts retention times for peptides via the svm_model that is trained by RTModel")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>",""," input file in analysisXML format");
			registerStringOption_("out","<file>","","output file in analysisXML format");
			registerStringOption_("svm_model","<file>","","svm model in libsvm format (can be produced by RTModel)");
			registerDoubleOption_("total_gradient_time","<time>",0.0,"the time (in seconds) of the gradient");
		}

		ExitCodes main_(int , char**)
		{
			AnalysisXMLFile analysisXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<IdentificationData> identifications;
			vector< String > peptides;
			vector< DoubleReal > training_retention_times;
			vector<PeptideHit> temp_peptide_hits;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			vector< vector< pair<SignedInt, DoubleReal> > >* encoded_composition_vectors;
			vector<svm_node*>* encoded_LibSVM_vectors;
			vector<DoubleReal>* predicted_retention_times;
			map< String, DoubleReal > predicted_data;
  
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
			String inputfile_name = getStringOption_("in");			
			inputFileReadable_(inputfile_name);
			String outputfile_name = getStringOption_("out");	
			outputFileWritable_(outputfile_name);
			String svmfile_name = getStringOption_("svm_model");
			inputFileReadable_(svmfile_name);			
			Real total_gradient_time = getDoubleOption_("total_gradient_time");

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			analysisXML_file.load(inputfile_name, protein_identifications, identifications);
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for(UnsignedInt i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].id.getPeptideHits();
				for(UnsignedInt j = 0; j < temp_peptide_hits.size(); j++)
				{
					peptides.push_back(temp_peptide_hits[j].getSequence());
				}
			}
			
			encoded_composition_vectors =  encoder.encodeCompositionVectors(peptides,  allowed_amino_acid_characters);
			encoded_LibSVM_vectors = encoder.encodeLibSVMVectors(*encoded_composition_vectors);
		
			svm.loadModel(svmfile_name);	
			predicted_retention_times = svm.predict(*encoded_LibSVM_vectors);
		
			delete encoded_composition_vectors;
			delete encoded_LibSVM_vectors;
		
			for(UnsignedInt i = 0; i < peptides.size(); i++)
			{
				predicted_data.insert(make_pair(peptides[i], ((*predicted_retention_times)[i] * total_gradient_time)));
			}
		
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			analysisXML_file.store(outputfile_name,
														 protein_identifications,
														 identifications,
														 predicted_data,
														 svm.getSVRProbability());
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPRTPredict tool;
	return tool.main(argc,argv);
}
  
/// @endcond





