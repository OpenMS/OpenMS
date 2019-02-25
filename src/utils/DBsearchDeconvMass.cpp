//
// Created by JihyungKim on 2019-02-14.
//

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

class DBsearchDeconvMass:
        public TOPPBase{

public:
    DBsearchDeconvMass():
    TOPPBase("DBsearchDeconvMass", "DB search on deconv masses", false)
    {}

    struct DB{
        vector<string> accessions;
        vector<string> proteinSeq;
        vector<double> mass;
    };

    struct Node{
        double mass;
        Node* prev;

        Node(const double& m, Node* p) :
                mass(m), prev(p) {}
        Node() :
                mass(0.0), prev(nullptr) {}
    };

protected:
    void registerOptionsAndFlags_() override{
        registerInputFile_("in", "<file>", "", "Input file");
    }

    ExitCodes main_(int, const char **) override {
        //-------------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
        String infilePath = getStringOption_("in");
        String outfilePath = getStringOption_("out");



        return EXECUTION_OK;
    }

    int calNomialMass(double &mass){
        return lrint(mass * 0.999497);
    }

    void getMassesFromDBseq(vector<string>& AAseq, vector<double>& masses){

    }

    void dynamicprogramming(vector<double>& dbmasses, vector<double>& ptmlist, vector<double>& candidatelist, double minM, double maxM){
        double tmp_mass = 0.0;

        vector<Node*> EndNodes;

        for(auto& mass : dbmasses){
            Node* head, tmp, newn;  // head node (protein seq node)
            head = new Node(mass, nullptr);

        }

    }

    void generatePTM_linkedlist(vector<double>& ptm_masslist, double minM, double maxM, vector<Node*>& ptm_linkedlist){
        vector <double>::iterator ptm_ptr = ptm_masslist.begin();

        for(auto& ptm: ptm_masslist){


            ptm_ptr++;
        }
    }
};