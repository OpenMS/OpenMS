#ifndef __CONFIG_SPECANNOTATE_H__
#define __CONFIG_SPECANNOTATE_H__

#include <OpenMS/config.h>

//QT-usage in method "annotate()"
#define ANNOTATE_QT

//No database: just the peakwise-cormen method using information saved in an xml file (also specified here,
// but has to be the same name as the executable of this program, only with extension .ini. and has to be in same directory as the exec.
// the specification here is used by the functional classes.
#ifndef DB_DEF   
   #define ANNOTATE_XML
   #define XML_FILE "TOPPView.ini"
#endif			

//DATABASE:
//QT-Database-Driver
#define QTDATABASEDRIVER "QMYSQL3"

//common:
#define DATABASE            "specannotate"

//Enzyme:
#define ENZ_TABLE           "enzyme"

//AminoAcid:
#define AMINO_TABLE         "aminoacid"

//MyElement:
#define ELEM_TABLE          "element"
#define ISO_TABLE           "isotope"

//AminoAcidModification:
#define MOD_TABLE           "modification"

//ProteinDigest
#define PROTEIN_TABLE       "protein"
#define SEQUENCE_TABLE      "sequence"
#define FRAGMENT_TABLE      "digest_fragment"
#define REALIZED_MOD_TABLE  "realized_modification"
#define MOD_COMB_TABLE      "modification_combination"
#define REAL_MOD_PLESS_TAB  "realized_modification_positionless"
#define MOD_COMB_PLESS_TAB  "modification_combination_positionless"

//Sample
#define PROT_MOD_SCEN_TABLE "protein_modification_scenario"
#define SAMPLE_TABLE        "sample"
#define ANNOTATION_TABLE    "annotation"




#endif
