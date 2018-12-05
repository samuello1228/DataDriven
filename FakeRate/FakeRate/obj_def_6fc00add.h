#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>
using std::string;

struct PAR0{
  float pt;
  float eta;
  float phi;
};
const string PAR0_s = "pt/F:eta/F:phi/F";

struct R_PAR:PAR0{
  float m;
};
const string R_PAR_s = PAR0_s+":m/F";

struct L_PAR:PAR0{
  float mT;
  unsigned int lFlag;
  int ID;
  int truthI;
  int truthType;
  int truthOrig;
  int firstEgMotherPdgId;
  int lepTruth;
};
const string L_PAR_s = PAR0_s+":mT/F:lFlag/i:ID/I:truthI/I:truthType/I:truthOrig/I:firstEgMotherPdgId/I:lepTruth/I";

struct EL_Par:L_PAR{
  int elID;
};
const string EL_PAR_s = L_PAR_s+":elID/i";

struct MU_Par:L_PAR{
  int muID;
};
const string MU_PAR_s = L_PAR_s+":muID/i";

struct J_PAR:PAR0{
  float m;
  unsigned int jFlag ;
};
const string J_PAR_s = PAR0_s+":m/F:jFlag/I";

struct EVT{
  unsigned long int run;
  unsigned long int event;
  unsigned int isMC;
  int flag;
  float weight;
  float averageMu;
  float pwt;
  float ElSF;
  float MuSF;
  float BtagSF;
  float JvtSF;
  float trigSF;
  float trigSF_BL;
  float qFwt;
  float fLwt;
 };
const string EVT_s = "run/l:event/l:isMC/i:flag/I:weight/F:averageMu/F:pwt/F:ElSF/F:MuSF/F:BtagSF/F:JvtSF/F:trigSF/F:trigSF_BL/F:qFwt/F:fLwt/F";

struct SIGNATURE{
  unsigned long int trigCode;//trigger info
  int nBaseLep;
  int nSigLep;
  int nJet;
  int nBJet;
  int isSS;
  int JetCut;
  int isZ;
  float Met;
  float MetX;
  float MetY;
  float MetPhi;
  float mT2;
  float HT;
  float mlj;
};
const string SIGNATURE_s = "trigCode/l:nBaseLep/I:nSigLep/I:nJet/I:nBJet/I:isSS/I:JetCut/I:isZ/I:Met/F:MetX/F:MetY/F:MetPhi/F:mT2/F:HT/F:mlj";

struct TR_PAR:PAR0{
  int pdgId; 
  int barcode;
  int motherI;
  int matchI;
  unsigned char status;
};
const string TR_PAR_s = PAR0_s+"pdgId/I:barcode/I:motherI/I:matchI/I:status/b";

enum FLAGS{
  PASS_GRL = 1<<0,
};

enum OFLAGS{
  IS_BASELINE = 1<<0, //1
  IS_SIGNAL = 1<<1, //2
  IS_PASSOR = 1<<2, //4
  IS_BAD = 1<<3, //8
  IS_LOOSEBASELINE = 1<<4, //16 
  ISO_LOOSE = 1<<8, /// isolation
  ISO_TIGHT = 1<<9,
  ISO_GRAD = 1<<10,
};

enum EFLAGS{
  EL_OQ = 1<<4, //16
};

enum MFLAGS{
  MU_COSMIC = 1<<4, //16
  TRIGGER_MATCHED = 1<<6,
};

enum JFALGS{
  JT_BJET_LOOSE = 1<<4, //16
  JT_BJET = 1<<5, //32
  JT_CJET = 1<<11,
};

typedef std::vector< L_PAR > LEPTONS;
typedef std::vector< PAR0 > PARTICLES;
typedef std::vector< J_PAR > JETS;
typedef std::vector< TR_PAR > TRUTHS;

#endif
