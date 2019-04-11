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
        registerInputFile_("in", "<FeatureFinderIntact result file>", "", "Input FeatureFinderIntact result file");
        setValidFormats_("in", ListUtils::create<String>("tsv"));

        registerInputFile_("d", "<Protein DB file>", "", "Input protein DB file");
        setValidFormats_("d", ListUtils::create<String>("fasta"));

        registerDoubleOption_("tol", "<tolerance>", 3.0, "proteoform mass tolerance (Da)", false, false);
        registerIntOption_("maxptm", "<max number of PTMs>", 10, "maximum number of PTMs per proteoform", false, false);

    }

    struct PTM{ // struct for PTM
        string name;
        double mass;
        int numOfmod; // for multimod (1 if singlemod)

        PTM(const string n, const double&m, const int nm):
                name(n), mass(m), numOfmod(nm) {}
        PTM():
                name(""), mass(0.0), numOfmod(0) {}
    };

    struct PTMnode{ // Node for PTM
        PTM ptm; // change this to pointer?
        double mass; // masses of modifications up to this ptm
        int level; // # of ptm up to this ptm.
        PTMnode* prev;

        PTMnode(PTM tPtm, const double& m, int l, PTMnode* p) :
                ptm(tPtm), mass(m), level(l), prev(p) {}
        PTMnode() :
                ptm({"", 0.0, 0}), mass(0.0), level(0), prev(nullptr)  {}
    };

    static bool compareNodeByMass(const PTMnode *a, const PTMnode *b) {
        return a->mass < b->mass;
    }

    struct Proteoform{ // struct for proteoforms
        string accessions;
        PTMnode* ptm_node;
        int mass; // nominal mass of proteoform
        double protein_mass;

        Proteoform(string protein_name, PTMnode* ptms, int theo_mass, double prot_mass):
                accessions(protein_name), ptm_node(ptms), mass(theo_mass), protein_mass(prot_mass) {}
    };

    struct PrtfNode{ // Node for proteoform
        Proteoform* prtf; // key

        PrtfNode *left;
        PrtfNode *right;
        int height;

        PrtfNode(Proteoform* ptrf_struct):
                prtf(ptrf_struct), left(nullptr), right(
                nullptr), height(1) {}
    };

    static int height(PrtfNode *N)
    {
        if (N == nullptr)
            return 0;
        return N->height;
    }

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
        fsout << "DeconvMass\tProteoformMass\tProteinAccession\tProteinMass\tTotalPTMmass\tPTMs" << endl;

        //-------------------------------------------------------------
        // generate PTM nodes vector
        //-------------------------------------------------------------
        vector<PTM> ptm_list;
        getMassesFromUNIMOD(modfilePath, ptm_list, 10);

        vector<PTMnode*> ptm_linkedlist;
        generatePTM_linkedlist(ptm_list, max_num_ptm, ptm_linkedlist);
        LOG_INFO << "ptm linked list # : " << ptm_linkedlist.size() << endl;

        sort(ptm_linkedlist.begin(), ptm_linkedlist.end(), compareNodeByMass);

        //-------------------------------------------------------------
        // generate theoretical proteoform tree
        //-------------------------------------------------------------
        LOG_INFO << "Building a proteoform tree..." << endl;
        PrtfNode* tree = nullptr;
        buildProteoformTree(fastaFilePath, ptm_linkedlist, tree);

        //-------------------------------------------------------------
        // Search tree with deconv masses & write output
        //-------------------------------------------------------------
        LOG_INFO << "Writing final results..." << endl;
        searchProteoformTree(deconvfilePath, tree, tolerance, fsout);
        fsout.close();

        return EXECUTION_OK;
    }

    string write_PTMlist_fromLinkedList(PTMnode*& ptm_node) {
        string out = "";
        string ptms_total_mass = to_string(ptm_node->mass);
        PTMnode* tmp = ptm_node;
        while(tmp != nullptr){
            out = tmp->ptm.name + "(" + to_string(tmp->ptm.numOfmod) + ")," + out;
            tmp = tmp->prev;
        }
        return out.substr(0, out.size()-1);
    }

    void searchProteoformTree(string deconvfilePath, PrtfNode*& tree, double tolerance, fstream &fs){
        vector<double> mass_vec;
        read_FeatureFinderIntact_result_file(deconvfilePath, mass_vec, "ExactMass");

        for(auto& mass : mass_vec){
            if (tolerance == 0) {
                PrtfNode* tmp = search_exact_node(mass, tree);
                if(tmp->prtf!= nullptr){
                    fs << write_proteoform_result(tmp, mass) << endl;
                }
            }else{
                vector <PrtfNode*> result_vec;
                search_within_tolerance(mass, 3.0, tree, result_vec);
                for(auto& tmp : result_vec){
                    fs << write_proteoform_result(tmp, mass) << endl;
                }
            }
            fs.flush();
        }
    }

    void read_FeatureFinderIntact_result_file(string deconvfilePath, vector<double>& mass_vec, string col_name){

        std::ifstream infile(deconvfilePath);
        string header;
        std::getline(infile, header);
        int index_of_col = get_index_of_tsv_header(header, col_name);

        string line;
        while (std::getline(infile, line)) {
            vector<string> results;
            boost::split(results, line, boost::is_any_of("\t"));
            mass_vec.push_back(String(results[index_of_col]).toDouble());
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

    struct PrtfNode* search_exact_node(int value, struct PrtfNode* tree){
        if (tree == nullptr)
            return nullptr;

        if(tree->prtf->mass == value) {
            return tree;
        }
        else if(tree->prtf->mass < value) {
            search_exact_node(value, tree->right);
        }
        else if(tree->prtf->mass > value) {
            search_exact_node(value, tree->left);
        }
        return nullptr;
    }

    void search_within_tolerance(double deconv_mass, double tol, struct PrtfNode *tree, vector<PrtfNode*>& result_vec){ // default value:3 Da
        if (tree == nullptr)
            return;

        result_vec.clear();
        double minM = deconv_mass-tol;
        double maxM = deconv_mass+tol;

        // get result by traversing from left-most to right-most
        traversal_tol(minM, maxM, tree, result_vec);
    }

    String write_proteoform_result(PrtfNode*& proteoform_node, double deconv_mass){
        String out;
        if( proteoform_node->prtf->ptm_node == nullptr){ // unmodified proteoform
            out = to_string(deconv_mass) + "\t" + to_string(proteoform_node->prtf->mass) + "\t"
                  + proteoform_node->prtf->accessions + "\t" + to_string(proteoform_node->prtf->protein_mass) + "\t"
                  + "0\t";
            return out;
        }
        string ptmlists = write_PTMlist_fromLinkedList(proteoform_node->prtf->ptm_node);
        // header : DeconvMass	ProteoformMass	ProteinMass	TotalPTMmass	PTMs
        out = to_string(deconv_mass) + "\t" + to_string(proteoform_node->prtf->mass) + "\t"
                + proteoform_node->prtf->accessions + "\t" + to_string(proteoform_node->prtf->protein_mass) + "\t"
                + to_string(proteoform_node->prtf->ptm_node->mass) + "\t" + ptmlists;
        return out;
    }

    void traversal_tol(double min, double max, struct PrtfNode *node, vector<PrtfNode*>& result){

        if (node == nullptr)
            return;

        double curmass = node->prtf->mass;
        if ( min < curmass )
            traversal_tol(min, max, node->left, result);

        if ( min <= curmass && max >= curmass ) {
            result.push_back(node);
        }

        if ( max > curmass )
            traversal_tol(min, max, node->right, result);
    }

    struct PrtfNode *rightRotate(struct PrtfNode *y) // A utility function to right rotate subtree rooted with y
    {
        struct PrtfNode *x = y->left;
        struct PrtfNode *t = x->right;

        // Perform rotation
        x->right = y;
        y->left = t;

        // Update heights
        y->height = max(height(y->left), height(y->right))+1;
        x->height = max(height(x->left), height(x->right))+1;

        // Return new root
        return x;
    }

    struct PrtfNode *leftRotate(struct PrtfNode *x) // A utility function to left rotate subtree rooted with x
    {
        struct PrtfNode *y = x->right;
        struct PrtfNode *t = y->left;

        // Perform rotation
        y->left = x;
        x->right = t;

        // Update heights
        x->height = max(height(x->left), height(x->right))+1;
        y->height = max(height(y->left), height(y->right))+1;

        // Return new root
        return y;
    }

    // Get Balance factor of node N
    int getBalance(struct PrtfNode *N)
    {
        if (N == nullptr)
            return 0;
        return height(N->left) - height(N->right);
    }

    struct PrtfNode* insert(struct PrtfNode* node, struct Proteoform* proteoform_ptr)
    {
        /* 1. Perform the normal BST insertion */
        if (node == nullptr){
            return( new PrtfNode(proteoform_ptr) );
        }

        int node_mass = node->prtf->mass;
        int new_mass = proteoform_ptr->mass;

        if (new_mass < node_mass)
            node->left = insert(node->left, proteoform_ptr);
        else if (new_mass > node_mass)
            node->right = insert(node->right, proteoform_ptr);
        else // Equal keys are not allowed in BST
            return node;

        /* 2. Update height of this ancestor node */
        node->height = 1 + max(height(node->left), height(node->right));

        /* 3. Get the balance factor of this ancestor
            node to check whether this node became
            unbalanced */
        int balance = getBalance(node);

        // If this node becomes unbalanced, then
        // there are 4 cases

        // Left Left Case
        if (balance > 1 && new_mass < node->left->prtf->mass)
            return rightRotate(node);

        // Right Right Case
        if (balance < -1 && new_mass > node->right->prtf->mass)
            return leftRotate(node);

        // Left Right Case
        if (balance > 1 && new_mass > node->left->prtf->mass)
        {
            node->left = leftRotate(node->left);
            return rightRotate(node);
        }

        // Right Left Case
        if (balance < -1 && new_mass < node->right->prtf->mass)
        {
            node->right = rightRotate(node->right);
            return leftRotate(node);
        }

        /* return the (unchanged) node pointer */
        return node;
    }

    void preOrder(struct PrtfNode *root) // A utility function to print preorder traversal of the tree // The function also prints height of every node
    {
        if(root != nullptr)
        {
            printf("%d ", root->prtf->mass);
            preOrder(root->left);
            preOrder(root->right);
        }
    }

    int calNomialMass(double &mass){
        return lrint(mass * 0.999497);
    }

    void buildProteoformTree(string DBfilepath, vector<PTMnode*>& ptmlist, PrtfNode*& tree){
        tree = nullptr;

        // get proteoform
        vector<FASTAFile::FASTAEntry> fentry;
        FASTAFile ff;
        ff.load(DBfilepath, fentry);

        // build tree
        for(auto& entry : fentry) {
            double seq_weight = AASequence::fromString(entry.sequence).getMonoWeight();

//            cout << entry.identifier << " : " << to_string(seq_weight) << endl;
//            double seq_weight_avg = AASequence::fromString(entry.sequence).getAverageWeight();
            tree = insert(tree, new Proteoform(entry.identifier, nullptr ,seq_weight, seq_weight)); // unmodified proteoform

            for (auto &ptm : ptmlist) {
                double total_weight = ptm->mass + seq_weight;
                int nominal_mass = calNomialMass(total_weight) ; // nominal mass
//                if (nominal_mass < minM) continue;
//                else if (nominal_mass > maxM) break; // ptmlist is sorted in ascending order

                Proteoform *prtfrm = new Proteoform(entry.identifier, ptm, nominal_mass, seq_weight);
                tree = insert(tree, prtfrm);
            }
        }
//        preOrder(tree);
    }

    void getMassesFromUNIMOD(string unimod_filepath, vector<PTM>& ptm_list, int multimod_limit){
        multimod_limit+=1;

        // File handling (should change when transferring this to OPENMS)
        UnimodXMLFile* file_ptr = new UnimodXMLFile();

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

            // multipmod masses (10 masses per one ptm)
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
            PTMnode* new_node;
            new_node = new PTMnode(ptm, ptm.mass, ptm.numOfmod, nullptr);

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
        }
    }
};
    int main(int argc, const char **argv){
        DBsearchDeconvMass tool;
        return tool.main(argc, argv);
    }
