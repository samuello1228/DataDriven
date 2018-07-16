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

struct PAR:PAR0{
  float MET_dPhi;
};
const string PAR_s = PAR0_s+":MET_dPhi/F";

struct R_PAR:PAR{
  float m;
  float dPhi;
  float dR;
  float jet0_dPhi;
};
const string R_PAR_s = PAR_s+":m/F:dPhi:dR:jet0_dPhi";

struct L_PAR:PAR{
  float mT;
  float jet0_dR;
  float jet_dRm;
  int ID;
  int truthI;
  unsigned int lFlag;
  unsigned int isoPass;
  int truthType;
  int truthOrig;
  int firstEgMotherPdgId;
};
const string L_PAR_s = PAR_s+":mT/F:jet0_dR:jet_dRm:ID/I:truthI/I:lFlag/i:isoPass/i:truthType/I:truthOrig/I:firstEgMotherPdgId/I";

struct EL_Par:L_PAR{
  int elID;
};
const string EL_PAR_s = L_PAR_s+":elID/i";

struct MU_Par:L_PAR{
  int muID;
};
const string MU_PAR_s = L_PAR_s+":muID/i";

struct J_PAR:PAR{
  float m;
  unsigned int jFlag ;
};
const string J_PAR_s = PAR_s+":m/F:jFlag/I";

struct EVT{
  unsigned long int run;
  unsigned long int event;
  unsigned int isMC;
  unsigned int cut;
  int flag;
//   float actualMu;
  float weight;
  float averageMu;
  float pwt;
  float ElSF;
  float MuSF;
  float BtagSF;
  float JvtSF;
  float trigSF;
  float qFwt;
  float qFwt_sys_1up;
  float qFwt_sys_1dn;
  float fLwt;
  float fLwt_e_sys_1up;
  float fLwt_e_sys_1dn;
  float fLwt_u_sys_1up;
  float fLwt_u_sys_1dn;
 };
// const string EVT_s = "run/l:event/l:isMC/i:cut/i;flag/I:actualMu/F:averageMu/F:weight/F:pwt/F:ElSF/F:MuSF/F:BtagSF/F:qFwt/F:qFwt_sys_1up/F:qFwt_sys_1dn/F:fLwt/F:fLwt_e_sys_1up/F:fLwt_e_sys_1dn/F:fLwt_u_sys_1up/F:fLwt_u_sys_1dn/F";
const string EVT_s = "run/l:event/l:isMC/i:cut/i:flag/I:weight/F:averageMu/F:pwt/F:ElSF/F:MuSF/F:BtagSF/F:JvtSF/F:trigSF/F:qFwt/F:qFwt_sys_1up/F:qFwt_sys_1dn/F:fLwt/F:fLwt_e_sys_1up/F:fLwt_e_sys_1dn/F:fLwt_u_sys_1up/F:fLwt_u_sys_1dn/F";

struct SIGNATURE{
  unsigned long int trigCode;//trigger info
  unsigned long int trigMask;//trigger info
  float Met;
  float MetRel;
  float MetX;
  float MetY;
  float mT2;
  float HT;
  float mjj;
  float mlj;
};
const string SIGNATURE_s = "trigCode/l:trigMask:Met/F:MetRel/F:MetX/F:MetY/F:mT2/F:HT/F:mjj:mlj";

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
};

enum JFALGS{
  JT_BJET_LOOSE = 1<<4, //16
  JT_BJET = 1<<5, //32
  JT_CJET = 1<<11,
};

typedef std::vector< L_PAR > LEPTONS;
typedef std::vector< PAR > PARTICLES;
typedef std::vector< J_PAR > JETS;
typedef std::vector< TR_PAR > TRUTHS;

#endif
