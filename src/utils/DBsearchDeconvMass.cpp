//
// Created by JihyungKim on 2019-02-14.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>

using namespace OpenMS;
using namespace std;
using namespace Internal;


class DBsearchDeconvMass:
        public TOPPBase{

public:
    DBsearchDeconvMass():
    TOPPBase("DBsearchDeconvMass", "DB search on deconv masses", false)
    {}

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
//        string proteinSeq;
        PTMnode* ptm_node;
        double mass;

        Proteoform(string protein_name, PTMnode* ptms, double theo_mass):
                accessions(protein_name), ptm_node(ptms), mass(theo_mass) {}
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

protected:
    void registerOptionsAndFlags_() override{
        registerInputFile_("in", "<file>", "", "Input file");
    }

    ExitCodes main_(int, const char **) override {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
//        String infilePath = getStringOption_("in");
//        String outfilePath = getStringOption_("out");
        String modfilePath = "CHEMISTRY/unimod.xml";
        String fastaFilePath = "/Users/jeek/Documents/A4B/FLASHDeconv_dy/CYS_BOVIN_Homolog.fasta";
//        String deconvfilePath = "";
        double minMass = 2000;
        double maxMass = 50000;

        //-------------------------------------------------------------
        // generate PTM nodes vector
        //-------------------------------------------------------------
        vector<PTM> ptm_list;
        getMassesFromUNIMOD(modfilePath, ptm_list, 2);
        cout << "ptm # : " << ptm_list.size() << endl;
//        for (auto& mod : ptm_list){
//            cout << mod.name << "\t" << mod.mass << endl;
//        }

        vector<PTMnode*> ptm_linkedlist;
        generatePTM_linkedlist(ptm_list, maxMass, ptm_linkedlist);
        cout << "ptm # : " << ptm_linkedlist.size() << endl;

        sort(ptm_linkedlist.begin(), ptm_linkedlist.end(), compareNodeByMass);

        for(auto& ends : ptm_linkedlist){
            string out = "";
            PTMnode* tmp = ends;
            while(tmp != nullptr){
                out = tmp->ptm.name + " * " + to_string(tmp->ptm.numOfmod) +" + " + out;

                tmp = tmp->prev;
            }
            cout << out << endl;
        }

        //-------------------------------------------------------------
        // generate theoretical proteoform tree
        //-------------------------------------------------------------
        PrtfNode* tree = nullptr;
        getProteoformTree(fastaFilePath, ptm_linkedlist, minMass, maxMass, tree);


        //-------------------------------------------------------------
        // generate PTM nodes vector
        //-------------------------------------------------------------

        return EXECUTION_OK;
    }

    struct PrtfNode *rightRotate(struct PrtfNode *y) // A utility function to right rotate subtree rooted with y
    {
        struct PrtfNode *x = y->left;
        struct PrtfNode *t = x->right;

        // Perform rotation
        x->right = y;
        y->left = t;

        // Update heights
        y->height = max(y->left->height, y->right->height)+1;
        x->height = max(x->left->height, x->right->height)+1;

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
        x->height = max(x->left->height, x->right->height)+1;
        y->height = max(y->left->height, y->right->height)+1;

        // Return new root
        return y;
    }

    // Get Balance factor of node N
    int getBalance(struct PrtfNode *N)
    {
        if (N == nullptr)
            return 0;
        return N->left->height - N->right->height;
    }

    struct PrtfNode* insert(struct PrtfNode* node, struct Proteoform* proteoform_ptr)
    {
        /* 1. Perform the normal BST insertion */
        if (node == nullptr)
            return( new PrtfNode(proteoform_ptr) );

        double node_mass = node->prtf->mass;
        double new_mass = proteoform_ptr->mass;

        if (new_mass < node_mass)
            node->left = insert(node->left, proteoform_ptr);
        else if (new_mass > node_mass)
            node->right = insert(node->right, proteoform_ptr);
        else // Equal keys are not allowed in BST
            return node;

        /* 2. Update height of this ancestor node */
        node->height = 1 + max(node->left->height,
                               node->right->height);

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
            printf("%f ", root->prtf->mass);
            preOrder(root->left);
            preOrder(root->right);
        }
    }

    int calNomialMass(double &mass){
        return lrint(mass * 0.999497);
    }

    void getProteoformTree(string DBfilepath, vector<PTMnode*>& ptmlist, double minM, double maxM, PrtfNode* tree){
        tree = nullptr;

        // get proteoform
        vector<FASTAFile::FASTAEntry> fentry;
        FASTAFile ff;
        ff.load(DBfilepath, fentry);

        // build tree
        for(auto& entry : fentry) {
            double seq_weight = AASequence::fromString(entry.sequence).getMonoWeight();

            for (auto &ptm : ptmlist) {
                double total_weight = ptm->mass + seq_weight;
                if (total_weight < minM)
                    continue;
                else if (total_weight > maxM)
                    break; // ptmlist is sorted in ascending order

                Proteoform *prtfrm = new Proteoform(entry.identifier, ptm, total_weight);
                tree = insert(tree, prtfrm);
            }
        }
        preOrder(tree);
    }

    void getMassesFromUNIMOD(string unimod_filepath, vector<PTM>& ptm_list, int multimod_limit){
        multimod_limit+=1;

        // File handling (should change when transferring this to OPENMS)
        UnimodXMLFile xml_file;
        UnimodXMLFile* ptr = new UnimodXMLFile();

        vector<ResidueModification*> modifications;
        ptr->load(unimod_filepath, modifications);

        cout << modifications.size() << endl;

        exit(0);

        // generate list
        vector <PTM> ptm_mass_cand;
        int cnt = 0;
        for(auto& mod : modifications)
        {
            if (cnt>3) break; // tmp purpose

            double mass = mod->getMonoMass();
            string name = mod->getId();

            std::cout.precision(5);
            cout << mod->getId() << "\t" <<  mod->getMonoMass() << mod->getAverageMass() << endl;

            PTM tmp = {name, mass, 1};
            ptm_list.push_back(tmp);

            // multipmod masses (10 masses per one ptm)
            for (int i=2; i<multimod_limit; ++i)
                ptm_list.push_back({name, mass * i, i });

            ++cnt;
        }


    }

    void generatePTM_linkedlist(vector<PTM>& ptm_list, double maxM, vector<PTMnode*>& result){

//        vector<Node*> result; // end element of linked list
        // generate list
        for(auto& ptm: ptm_list){
            PTMnode* new_node;
            new_node = new PTMnode(ptm, ptm.mass, 1, nullptr);

            vector<PTMnode*> tmp_result;
            tmp_result.push_back(new_node);

            for (auto& res: result){
                if (res->ptm.name == ptm.name) // avoid linking same modification
                    continue;

                double curr_mass = res->mass + ptm.mass;
                if (curr_mass > maxM)
                    continue;

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
