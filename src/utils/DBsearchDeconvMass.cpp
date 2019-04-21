// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <unordered_set>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;


class DBsearchDeconvMass:
        public TOPPBase{

public:
    DBsearchDeconvMass():
    TOPPBase("DBsearchDeconvMass", "DB search on deconv masses", false)
    {}

protected:
    void registerOptionsAndFlags_() override{
        registerInputFile_("in", "<FeatureFinderIntact feature result file>", "", "Input FeatureFinderIntact result file");
        setValidFormats_("in", ListUtils::create<String>("tsv"));

        registerInputFile_("d", "<Protein DB file>", "", "Input protein DB file");
        setValidFormats_("d", ListUtils::create<String>("fasta"));

        registerDoubleOption_("tol", "<tolerance>", 0, "proteoform mass tolerance (Da)", false, false);
        registerIntOption_("maxptm", "<max number of PTMs>", 4, "maximum number of PTMs per proteoform", false, false);

    }

    struct PTM{ // struct for PTM
        string name;
        double mass;
        int numOfmod; // for multimod (1 if singlemod)

        PTM(string n, double m, int nm):
                name(n), mass(m), numOfmod(nm) {}
        PTM():
                name(""), mass(0.0), numOfmod(0) {}
    };

    struct PTMnode{ // Node for PTM
        PTM ptm; // change this to pointer?
        double mass; // masses of modifications up to this ptm
        int level; // # of ptm up to this ptm.
        PTMnode* prev;

        PTMnode(PTM tPtm, double m, int l, PTMnode* p) :
                ptm(tPtm), mass(m), level(l), prev(p) {}
        PTMnode() :
                ptm({"", 0.0, 0}), mass(0.0), level(0), prev(nullptr)  {}
    };

    static bool compareNodeByMass(const PTMnode *a, const PTMnode *b) {
        return a->mass < b->mass;
    }

    struct Proteoform{ // struct for proteoforms
        string accessions;
        int mass; // nominal mass of proteoform
        double protein_mass;
        Size PTM_index;

        Proteoform(string protein_name, int theo_mass, double prot_mass, Size ptm_index):
                accessions(protein_name), mass(theo_mass), protein_mass(prot_mass), PTM_index(ptm_index){}
    };

    ExitCodes main_(int, const char **) override {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
        String modfilePath = "CHEMISTRY/unimod_most_freq.xml";

        string deconvfilePath = getStringOption_("in");
        String fastaFilePath = getStringOption_("d");
//        String fastaFilePath = "/Users/jeek/Documents/A4B/FLASHDeconv_dy/CYS_BOVIN.fasta";
//        String deconvfilePath = "/Users/jeek/Documents/A4B/FLASHDeconv_dy/190215_Cyto_untreated_30V_ISCID.tsv";
        int max_num_ptm = getIntOption_("maxptm");
        double tolerance = getDoubleOption_("tol");

        string outfilePath = deconvfilePath.substr(0, deconvfilePath.find_last_of(".")) + "_dbmass.tsv" ; // based on deconvfilePath
        cout << outfilePath << endl;

        fstream fsout;
        fsout.open(outfilePath, fstream::out);


        //-------------------------------------------------------------
        // db search on deconvoluted masses
        //-------------------------------------------------------------
        //--- read FeatureFinderIntact results---
        vector<double> mass_vec;
        vector<string> dfile_vec; // rows of deconvfilePath
        read_FeatureFinderIntact_result_file(deconvfilePath, mass_vec, dfile_vec, "ExactMass");
        double min_prot_mass = *min_element(mass_vec.begin(), mass_vec.end()) - 1000.0;
        double max_prot_mass = *max_element(mass_vec.begin(), mass_vec.end()) + 1000.0;

        vector<PTM> ptm_list;
        getMassesFromUNIMOD(modfilePath, ptm_list, max_num_ptm);

        vector<PTMnode*> ptm_linkedlist;
        generatePTM_linkedlist(ptm_list, max_num_ptm, ptm_linkedlist);
        LOG_INFO << "ptm linked list # : " << ptm_linkedlist.size() << endl;

        sort(ptm_linkedlist.begin(), ptm_linkedlist.end(), compareNodeByMass);

        LOG_INFO << "Building PTM lists for proteoforms" << endl;
        vector<double> PTM_masses_vec;
        vector<string> PTM_lists_vec;
        get_PTM_vectors(ptm_linkedlist, PTM_masses_vec, PTM_lists_vec);
        ptm_linkedlist.clear();

        //-------------------------------------------------------------
        // generate theoretical proteoform tree
        //-------------------------------------------------------------
        LOG_INFO << "Building a proteoform tree..." << endl;
        vector<Proteoform> tree;
        multimap<int, Size> multimap_proteoform;
        buildProteoformMap(fastaFilePath, PTM_masses_vec, tree, multimap_proteoform, min_prot_mass, max_prot_mass);
        LOG_INFO << "\n# proteoforms in map : " << multimap_proteoform.size() << endl;

        //-------------------------------------------------------------
        // Search tree with deconv masses & write output
        //-------------------------------------------------------------
        LOG_INFO << "Writing final results..." << endl;
        fsout << dfile_vec[0] << "\tDeconvMass\tProteoformMass\tProteinAccession\tProteinMass\tTotalPTMmass\tPTMs" << endl;

//        cout << multimap_proteoform.size() << "\t" << PTM_masses_vec.size() << "\t" << PTM_lists_vec.size() << "\t" << mass_vec.size() << endl;
        search_map(tree, multimap_proteoform, fsout, tolerance, PTM_masses_vec, PTM_lists_vec, mass_vec, dfile_vec);
        fsout.close();
        LOG_INFO << "\n" << endl;
        return EXECUTION_OK;
    }

    void get_PTM_vectors(vector<PTMnode*>& ptmlist, vector<double>& PTM_masss_vec, vector<string>& PTM_str_vec){
        // create PTM vectors using PTM linked lists.
        Size ptmlen = ptmlist.size();
        PTM_masss_vec.reserve(ptmlen);
        PTM_str_vec.reserve(ptmlen);
        for (auto &ptm : ptmlist) {
//            double ptm_total_weight = ptm->mass;
//            string ptmlists = write_PTMlist_fromLinkedList(ptm);
            PTM_masss_vec.push_back(ptm->mass);
            PTM_str_vec.push_back(write_PTMlist_fromLinkedList(ptm));
        }
    }

    string write_PTMlist_fromLinkedList(PTMnode*& ptm_node) {
        string out = "";
        string ptms_total_mass = to_string(ptm_node->mass);
        PTMnode* tmp = ptm_node;
        while(tmp != nullptr){
            out = tmp->ptm.name + "(" + to_string(tmp->ptm.numOfmod) + ")," + out;
            tmp = tmp->prev;
        }
        delete tmp;
        return out.substr(0, out.size()-1);
    }


    void read_FeatureFinderIntact_result_file(string deconvfilePath, vector<double>& mass_vec, vector<string>& file_lines, string col_name){

        std::ifstream infile(deconvfilePath);
        string header;
        std::getline(infile, header);
        file_lines.push_back(header);
        int index_of_col = get_index_of_tsv_header(header, col_name);

        string line;
        while (std::getline(infile, line)) {
            vector<string> results;
            boost::split(results, line, boost::is_any_of("\t"));
            mass_vec.push_back(String(results[index_of_col]).toDouble());
            file_lines.push_back(line);
        }
    }

    int get_index_of_tsv_header(string header_str, string col_name){

        vector<string> header_vec;
        std::stringstream lineStream(header_str);
        string tmp;
        int index = -1;
        while (std::getline(lineStream, tmp, '\t'))
        {
            header_vec.push_back(tmp);
            if(tmp==col_name){
                index = header_vec.size()-1;
                break;
            }
        }
        lineStream.clear();
        return index;
    }

    String write_proteoform_result(Proteoform& proteoform_node, double deconv_mass, vector<double>& ptm_mass, vector<string>& ptm_str){
        String out;
        if( proteoform_node.PTM_index == Size(-1)){ // unmodified proteoform
            out = to_string(deconv_mass) + "\t" + to_string(proteoform_node.mass) + "\t"
                  + proteoform_node.accessions + "\t" + to_string(proteoform_node.protein_mass) + "\t"
                  + "0\t";
            return out;
        }

        // header : DeconvMass	ProteoformMass	ProteinMass	TotalPTMmass	PTMs
        out = to_string(deconv_mass) + "\t" + to_string(proteoform_node.mass) + "\t"
                + proteoform_node.accessions + "\t" + to_string(proteoform_node.protein_mass) + "\t"
                + to_string(ptm_mass[proteoform_node.PTM_index]) + "\t" + ptm_str[proteoform_node.PTM_index];
        return out;
    }

    int calNomialMass(double &mass){
        return lrint(mass * 0.999497);
    }

    void buildProteoformMap(string DBfilepath, vector<double>& ptm_masses, vector<Proteoform>& tree, multimap<int, Size>& multimap_proteoform, double minM, double maxM){
        tree.clear();

        // get proteoform
        vector<FASTAFile::FASTAEntry> fentry;
        FASTAFile ff;
        ff.load(DBfilepath, fentry);
        LOG_INFO <<  "# DB seq : " <<fentry.size() << endl;

        tree.reserve( fentry.size() * (ptm_masses.size()+1) );
        Size index = 0;
        float prevProgress = .0;
        // build tree
        for(auto entry = fentry.begin(); entry != fentry.end(); ++entry) {
            if ( AASequence::fromString(entry->sequence).toString().find('X') != std::string::npos ) continue; // protein contains AA 'X' are ignored

            double seq_weight = AASequence::fromString(entry->sequence).getMonoWeight();
            if(seq_weight < minM || seq_weight > maxM) continue;

            float progress = (float) (entry - fentry.begin()) / fentry.size();
            if (progress > prevProgress + .01)
            {
                printProgress(progress); //
                prevProgress = progress;
            }
            Proteoform* tmp = new Proteoform(entry->identifier, seq_weight, seq_weight, Size(-1));
            tree.push_back(*tmp); // unmodified proteoform
            delete tmp;

            multimap_proteoform.insert(make_pair(calNomialMass(seq_weight), index++));

            for (Size i=0; i<ptm_masses.size(); i++) {
                double total_weight = ptm_masses[i] + seq_weight;
                int nominal_mass = calNomialMass(total_weight) ; // nominal mass

                tmp = new Proteoform(entry->identifier, nominal_mass, seq_weight, i);
                tree.push_back(*tmp);
                delete tmp;
                multimap_proteoform.insert(make_pair(nominal_mass, index++));
            }

        }
        printProgress(1);
    }

    static void printProgress(float progress)
    {
        int barWidth = 70;
        cout << "[";
        int pos = (int) (barWidth * progress);
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
            {
                std::cout << "=";
            }
            else if (i == pos)
            {
                std::cout << ">";
            }
            else
            {
                std::cout << " ";
            }
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();
    }

    void search_map(vector<Proteoform> &tree, multimap<int, Size> &multimap_proteoform, fstream &fs, double tolerance,
                           vector<double> &ptm_mass, vector<string> &ptm_str, vector<double> &mass_vec, const vector<string> &dfile_vec) {
        SignedSize count = 0; // keep track of the row number of deconv. result file

        float prevProgress = .0;
        for(auto mass = mass_vec.begin(); mass!=mass_vec.end(); ++mass) {
            int mass_nm = calNomialMass(*mass);

            multimap<int, Size>::const_iterator low_it;
            multimap<int, Size>::const_iterator up_it;

            low_it = multimap_proteoform.lower_bound(mass_nm - tolerance);
            up_it = multimap_proteoform.upper_bound(mass_nm - tolerance);

            count++;
            // no matching in data
            if (low_it == up_it) { continue; }

            float progress = (float) (mass - mass_vec.begin()) / mass_vec.size();
            if (progress > prevProgress + .01)
            {
                printProgress(progress); //
                prevProgress = progress;
            }

            // for all matching results
            for (; low_it != up_it; ++low_it) {
                fs << dfile_vec[count] << "\t" << write_proteoform_result(tree[low_it->second], mass_nm, ptm_mass, ptm_str) << endl;
            }
        }
        fs.flush();
        printProgress(1);
    }


    void getMassesFromUNIMOD(string unimod_filepath, vector<PTM>& ptm_list, int multimod_limit){
        multimod_limit+=1;

        // File handling (should change when transferring this to OPENMS)
        auto* file_ptr = new UnimodXMLFile();

        vector<ResidueModification*> modifications;
        file_ptr->load(unimod_filepath, modifications);

        // generate list
        unordered_set<double> ptmMhash;
        int cnt = 0;
        for(auto& mod : modifications) {
            if (ptmMhash.count(mod->getDiffMonoMass()) > 0) continue;
//            if (cnt>10) break; // tmp purpose

            double mass = mod->getDiffMonoMass();
            string name = mod->getId();

            /*
             * mod->getId() : Acetyl
             * mod->getName() : not working
             * mod->getFullName() : Acetylation
             * mod->getUniModAccession() : UniMod:1
             * mod->getMonoMass() : 0
             * mod->getDiffMonoMass() : 42.010565
             * */
            PTM tmp = {name, mass, 1};
            ptm_list.push_back(tmp);
            ptmMhash.insert(mass);

            // multipmod masses (multimod_limit per one ptm)
            double multi_mass = mass;
            for (int i = 2; i < multimod_limit; ++i) {
                multi_mass += mass;
                ptm_list.push_back({name, multi_mass , i});
            }
            ++cnt;
        }
        LOG_INFO << "# ptm types : " << cnt << endl;
    }

    void generatePTM_linkedlist(vector<PTM>& ptm_list, int max_num_ptm, vector<PTMnode*>& result){

        for(auto& ptm: ptm_list){
            PTMnode* new_node = new PTMnode(ptm, ptm.mass, ptm.numOfmod, nullptr);

            vector<PTMnode*> tmp_result;
            tmp_result.push_back(new_node);

            for (auto& res: result){
                if (res->ptm.name == ptm.name) // avoid linking same modification
                    continue;
                if (res->level+ptm.numOfmod >= max_num_ptm){ // avoid the number of ptm exceeding 10
                    continue;
                }

                double curr_mass = res->mass + ptm.mass;

                new_node = new PTMnode(ptm, curr_mass, res->level+ptm.numOfmod, res);
                tmp_result.push_back(new_node);
            }

            for (auto&tmp : tmp_result){
                result.push_back(tmp);
            }
            tmp_result.clear();
        }
    }
};
    int main(int argc, const char **argv){
        DBsearchDeconvMass tool;
        return tool.main(argc, argv);
    }
