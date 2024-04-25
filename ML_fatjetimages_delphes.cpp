
#include "SampleAnalyzer/User/Analyzer/ML_fatjetimages_delphes.h"
#include "/scratch/sj1n19/MG5/MG5_aMC_v3_3_1/HEPTools/madanalysis5/madanalysis5/tools/fastjet/include/fastjet/contrib/VariableRPlugin.hh" // to use variable R algorithm
#include "/scratch/sj1n19/MG5/MG5_aMC_v3_3_1/HEPTools/madanalysis5/madanalysis5/tools/fastjet/include/fastjet/contrib/fastjet_spectraljet/SpectralPlugin.hh"
#include "/scratch/sj1n19/MG5/MG5_aMC_v3_3_1/HEPTools/madanalysis5/madanalysis5/tools/fastjet/include/fastjet/ClusterSequence.hh"
#include "/scratch/sj1n19/MG5/MG5_aMC_v3_3_1/HEPTools/madanalysis5/madanalysis5/tools/fastjet/include/fastjet/PseudoJet.hh"
#include "/scratch/sj1n19/MG5/MG5_aMC_v3_3_1/HEPTools/madanalysis5/madanalysis5/tools/fastjet/include/fastjet/JetDefinition.hh"
#include "SampleAnalyzer/Interfaces/fastjet/ClusterAlgoFastJet.h"
#include "SampleAnalyzer/Commons/Service/LoopService.h"
#include "SampleAnalyzer/Commons/Service/Physics.h"
#include "SampleAnalyzer/Commons/Service/PDGService.h"
#include "SampleAnalyzer/Commons/Service/RandomService.h"
#include <algorithm>
#include <cmath>
#include <time.h>
using namespace MA5;
using namespace std;
using namespace fastjet;
using namespace contrib;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool ML_fatjetimages_delphes::Initialize(const MA5::Configuration& cfg,
                      const std::map<std::string,std::string>& parameters)
{
  INFO << "      <><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
  INFO << "      <>   Analysis: bjet_analysis_delphes            <>" << endmsg;
  INFO << "      <>   Recaster: S.Jain                           <>" << endmsg;
  INFO << "      <><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;

  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();

  // definition of the multiparticle "hadronic"
  PHYSICS->mcConfig().AddHadronicId(-20543);
  PHYSICS->mcConfig().AddHadronicId(-20533);
  PHYSICS->mcConfig().AddHadronicId(-20523);
  PHYSICS->mcConfig().AddHadronicId(-20513);
  PHYSICS->mcConfig().AddHadronicId(-20433);
  PHYSICS->mcConfig().AddHadronicId(-20423);
  PHYSICS->mcConfig().AddHadronicId(-20413);
  PHYSICS->mcConfig().AddHadronicId(-20323);
  PHYSICS->mcConfig().AddHadronicId(-20313);
  PHYSICS->mcConfig().AddHadronicId(-20213);
  PHYSICS->mcConfig().AddHadronicId(-10543);
  PHYSICS->mcConfig().AddHadronicId(-10541);
  PHYSICS->mcConfig().AddHadronicId(-10533);
  PHYSICS->mcConfig().AddHadronicId(-10531);
  PHYSICS->mcConfig().AddHadronicId(-10523);
  PHYSICS->mcConfig().AddHadronicId(-10521);
  PHYSICS->mcConfig().AddHadronicId(-10513);
  PHYSICS->mcConfig().AddHadronicId(-10511);
  PHYSICS->mcConfig().AddHadronicId(-10433);
  PHYSICS->mcConfig().AddHadronicId(-10431);
  PHYSICS->mcConfig().AddHadronicId(-10423);
  PHYSICS->mcConfig().AddHadronicId(-10421);
  PHYSICS->mcConfig().AddHadronicId(-10413);
  PHYSICS->mcConfig().AddHadronicId(-10411);
  PHYSICS->mcConfig().AddHadronicId(-10323);
  PHYSICS->mcConfig().AddHadronicId(-10321);
  PHYSICS->mcConfig().AddHadronicId(-10313);
  PHYSICS->mcConfig().AddHadronicId(-10311);
  PHYSICS->mcConfig().AddHadronicId(-10213);
  PHYSICS->mcConfig().AddHadronicId(-10211);
  PHYSICS->mcConfig().AddHadronicId(-5554);
  PHYSICS->mcConfig().AddHadronicId(-5544);
  PHYSICS->mcConfig().AddHadronicId(-5542);
  PHYSICS->mcConfig().AddHadronicId(-5534);
  PHYSICS->mcConfig().AddHadronicId(-5532);
  PHYSICS->mcConfig().AddHadronicId(-5524);
  PHYSICS->mcConfig().AddHadronicId(-5522);
  PHYSICS->mcConfig().AddHadronicId(-5514);
  PHYSICS->mcConfig().AddHadronicId(-5512);
  PHYSICS->mcConfig().AddHadronicId(-5503);
  PHYSICS->mcConfig().AddHadronicId(-5444);
  PHYSICS->mcConfig().AddHadronicId(-5442);
  PHYSICS->mcConfig().AddHadronicId(-5434);
  PHYSICS->mcConfig().AddHadronicId(-5432);
  PHYSICS->mcConfig().AddHadronicId(-5424);
  PHYSICS->mcConfig().AddHadronicId(-5422);
  PHYSICS->mcConfig().AddHadronicId(-5414);
  PHYSICS->mcConfig().AddHadronicId(-5412);
  PHYSICS->mcConfig().AddHadronicId(-5403);
  PHYSICS->mcConfig().AddHadronicId(-5401);
  PHYSICS->mcConfig().AddHadronicId(-5342);
  PHYSICS->mcConfig().AddHadronicId(-5334);
  PHYSICS->mcConfig().AddHadronicId(-5332);
  PHYSICS->mcConfig().AddHadronicId(-5324);
  PHYSICS->mcConfig().AddHadronicId(-5322);
  PHYSICS->mcConfig().AddHadronicId(-5314);
  PHYSICS->mcConfig().AddHadronicId(-5312);
  PHYSICS->mcConfig().AddHadronicId(-5303);
  PHYSICS->mcConfig().AddHadronicId(-5301);
  PHYSICS->mcConfig().AddHadronicId(-5242);
  PHYSICS->mcConfig().AddHadronicId(-5232);
  PHYSICS->mcConfig().AddHadronicId(-5224);
  PHYSICS->mcConfig().AddHadronicId(-5222);
  PHYSICS->mcConfig().AddHadronicId(-5214);
  PHYSICS->mcConfig().AddHadronicId(-5212);
  PHYSICS->mcConfig().AddHadronicId(-5203);
  PHYSICS->mcConfig().AddHadronicId(-5201);
  PHYSICS->mcConfig().AddHadronicId(-5142);
  PHYSICS->mcConfig().AddHadronicId(-5132);
  PHYSICS->mcConfig().AddHadronicId(-5122);
  PHYSICS->mcConfig().AddHadronicId(-5114);
  PHYSICS->mcConfig().AddHadronicId(-5112);
  PHYSICS->mcConfig().AddHadronicId(-5103);
  PHYSICS->mcConfig().AddHadronicId(-5101);
  PHYSICS->mcConfig().AddHadronicId(-4444);
  PHYSICS->mcConfig().AddHadronicId(-4434);
  PHYSICS->mcConfig().AddHadronicId(-4432);
  PHYSICS->mcConfig().AddHadronicId(-4424);
  PHYSICS->mcConfig().AddHadronicId(-4422);
  PHYSICS->mcConfig().AddHadronicId(-4414);
  PHYSICS->mcConfig().AddHadronicId(-4412);
  PHYSICS->mcConfig().AddHadronicId(-4403);
  PHYSICS->mcConfig().AddHadronicId(-4334);
  PHYSICS->mcConfig().AddHadronicId(-4332);
  PHYSICS->mcConfig().AddHadronicId(-4324);
  PHYSICS->mcConfig().AddHadronicId(-4322);
  PHYSICS->mcConfig().AddHadronicId(-4314);
  PHYSICS->mcConfig().AddHadronicId(-4312);
  PHYSICS->mcConfig().AddHadronicId(-4303);
  PHYSICS->mcConfig().AddHadronicId(-4301);
  PHYSICS->mcConfig().AddHadronicId(-4232);
  PHYSICS->mcConfig().AddHadronicId(-4224);
  PHYSICS->mcConfig().AddHadronicId(-4222);
  PHYSICS->mcConfig().AddHadronicId(-4214);
  PHYSICS->mcConfig().AddHadronicId(-4212);
  PHYSICS->mcConfig().AddHadronicId(-4203);
  PHYSICS->mcConfig().AddHadronicId(-4201);
  PHYSICS->mcConfig().AddHadronicId(-4132);
  PHYSICS->mcConfig().AddHadronicId(-4122);
  PHYSICS->mcConfig().AddHadronicId(-4114);
  PHYSICS->mcConfig().AddHadronicId(-4112);
  PHYSICS->mcConfig().AddHadronicId(-4103);
  PHYSICS->mcConfig().AddHadronicId(-4101);
  PHYSICS->mcConfig().AddHadronicId(-3334);
  PHYSICS->mcConfig().AddHadronicId(-3324);
  PHYSICS->mcConfig().AddHadronicId(-3322);
  PHYSICS->mcConfig().AddHadronicId(-3314);
  PHYSICS->mcConfig().AddHadronicId(-3312);
  PHYSICS->mcConfig().AddHadronicId(-3303);
  PHYSICS->mcConfig().AddHadronicId(-3224);
  PHYSICS->mcConfig().AddHadronicId(-3222);
  PHYSICS->mcConfig().AddHadronicId(-3214);
  PHYSICS->mcConfig().AddHadronicId(-3212);
  PHYSICS->mcConfig().AddHadronicId(-3203);
  PHYSICS->mcConfig().AddHadronicId(-3201);
  PHYSICS->mcConfig().AddHadronicId(-3122);
  PHYSICS->mcConfig().AddHadronicId(-3114);
  PHYSICS->mcConfig().AddHadronicId(-3112);
  PHYSICS->mcConfig().AddHadronicId(-3103);
  PHYSICS->mcConfig().AddHadronicId(-3101);
  PHYSICS->mcConfig().AddHadronicId(-2224);
  PHYSICS->mcConfig().AddHadronicId(-2214);
  PHYSICS->mcConfig().AddHadronicId(-2212);
  PHYSICS->mcConfig().AddHadronicId(-2203);
  PHYSICS->mcConfig().AddHadronicId(-2114);
  PHYSICS->mcConfig().AddHadronicId(-2112);
  PHYSICS->mcConfig().AddHadronicId(-2103);
  PHYSICS->mcConfig().AddHadronicId(-2101);
  PHYSICS->mcConfig().AddHadronicId(-1114);
  PHYSICS->mcConfig().AddHadronicId(-1103);
  PHYSICS->mcConfig().AddHadronicId(-545);
  PHYSICS->mcConfig().AddHadronicId(-543);
  PHYSICS->mcConfig().AddHadronicId(-541);
  PHYSICS->mcConfig().AddHadronicId(-535);
  PHYSICS->mcConfig().AddHadronicId(-533);
  PHYSICS->mcConfig().AddHadronicId(-531);
  PHYSICS->mcConfig().AddHadronicId(-525);
  PHYSICS->mcConfig().AddHadronicId(-523);
  PHYSICS->mcConfig().AddHadronicId(-521);
  PHYSICS->mcConfig().AddHadronicId(-515);
  PHYSICS->mcConfig().AddHadronicId(-513);
  PHYSICS->mcConfig().AddHadronicId(-511);
  PHYSICS->mcConfig().AddHadronicId(-435);
  PHYSICS->mcConfig().AddHadronicId(-433);
  PHYSICS->mcConfig().AddHadronicId(-431);
  PHYSICS->mcConfig().AddHadronicId(-425);
  PHYSICS->mcConfig().AddHadronicId(-423);
  PHYSICS->mcConfig().AddHadronicId(-421);
  PHYSICS->mcConfig().AddHadronicId(-415);
  PHYSICS->mcConfig().AddHadronicId(-413);
  PHYSICS->mcConfig().AddHadronicId(-411);
  PHYSICS->mcConfig().AddHadronicId(-325);
  PHYSICS->mcConfig().AddHadronicId(-323);
  PHYSICS->mcConfig().AddHadronicId(-321);
  PHYSICS->mcConfig().AddHadronicId(-315);
  PHYSICS->mcConfig().AddHadronicId(-313);
  PHYSICS->mcConfig().AddHadronicId(-311);
  PHYSICS->mcConfig().AddHadronicId(-215);
  PHYSICS->mcConfig().AddHadronicId(-213);
  PHYSICS->mcConfig().AddHadronicId(-211);
  PHYSICS->mcConfig().AddHadronicId(111);
  PHYSICS->mcConfig().AddHadronicId(113);
  PHYSICS->mcConfig().AddHadronicId(115);
  PHYSICS->mcConfig().AddHadronicId(130);
  PHYSICS->mcConfig().AddHadronicId(211);
  PHYSICS->mcConfig().AddHadronicId(213);
  PHYSICS->mcConfig().AddHadronicId(215);
  PHYSICS->mcConfig().AddHadronicId(221);
  PHYSICS->mcConfig().AddHadronicId(223);
  PHYSICS->mcConfig().AddHadronicId(225);
  PHYSICS->mcConfig().AddHadronicId(310);
  PHYSICS->mcConfig().AddHadronicId(311);
  PHYSICS->mcConfig().AddHadronicId(313);
  PHYSICS->mcConfig().AddHadronicId(315);
  PHYSICS->mcConfig().AddHadronicId(321);
  PHYSICS->mcConfig().AddHadronicId(323);
  PHYSICS->mcConfig().AddHadronicId(325);
  PHYSICS->mcConfig().AddHadronicId(331);
  PHYSICS->mcConfig().AddHadronicId(333);
  PHYSICS->mcConfig().AddHadronicId(335);
  PHYSICS->mcConfig().AddHadronicId(411);
  PHYSICS->mcConfig().AddHadronicId(413);
  PHYSICS->mcConfig().AddHadronicId(415);
  PHYSICS->mcConfig().AddHadronicId(421);
  PHYSICS->mcConfig().AddHadronicId(423);
  PHYSICS->mcConfig().AddHadronicId(425);
  PHYSICS->mcConfig().AddHadronicId(431);
  PHYSICS->mcConfig().AddHadronicId(433);
  PHYSICS->mcConfig().AddHadronicId(435);
  PHYSICS->mcConfig().AddHadronicId(441);
  PHYSICS->mcConfig().AddHadronicId(443);
  PHYSICS->mcConfig().AddHadronicId(445);
  PHYSICS->mcConfig().AddHadronicId(511);
  PHYSICS->mcConfig().AddHadronicId(513);
  PHYSICS->mcConfig().AddHadronicId(515);
  PHYSICS->mcConfig().AddHadronicId(521);
  PHYSICS->mcConfig().AddHadronicId(523);
  PHYSICS->mcConfig().AddHadronicId(525);
  PHYSICS->mcConfig().AddHadronicId(531);
  PHYSICS->mcConfig().AddHadronicId(533);
  PHYSICS->mcConfig().AddHadronicId(535);
  PHYSICS->mcConfig().AddHadronicId(541);
  PHYSICS->mcConfig().AddHadronicId(543);
  PHYSICS->mcConfig().AddHadronicId(545);
  PHYSICS->mcConfig().AddHadronicId(551);
  PHYSICS->mcConfig().AddHadronicId(553);
  PHYSICS->mcConfig().AddHadronicId(555);
  PHYSICS->mcConfig().AddHadronicId(1103);
  PHYSICS->mcConfig().AddHadronicId(1114);
  PHYSICS->mcConfig().AddHadronicId(2101);
  PHYSICS->mcConfig().AddHadronicId(2103);
  PHYSICS->mcConfig().AddHadronicId(2112);
  PHYSICS->mcConfig().AddHadronicId(2114);
  PHYSICS->mcConfig().AddHadronicId(2203);
  PHYSICS->mcConfig().AddHadronicId(2212);
  PHYSICS->mcConfig().AddHadronicId(2214);
  PHYSICS->mcConfig().AddHadronicId(2224);
  PHYSICS->mcConfig().AddHadronicId(3101);
  PHYSICS->mcConfig().AddHadronicId(3103);
  PHYSICS->mcConfig().AddHadronicId(3112);
  PHYSICS->mcConfig().AddHadronicId(3114);
  PHYSICS->mcConfig().AddHadronicId(3122);
  PHYSICS->mcConfig().AddHadronicId(3201);
  PHYSICS->mcConfig().AddHadronicId(3203);
  PHYSICS->mcConfig().AddHadronicId(3212);
  PHYSICS->mcConfig().AddHadronicId(3214);
  PHYSICS->mcConfig().AddHadronicId(3222);
  PHYSICS->mcConfig().AddHadronicId(3224);
  PHYSICS->mcConfig().AddHadronicId(3303);
  PHYSICS->mcConfig().AddHadronicId(3312);
  PHYSICS->mcConfig().AddHadronicId(3314);
  PHYSICS->mcConfig().AddHadronicId(3322);
  PHYSICS->mcConfig().AddHadronicId(3324);
  PHYSICS->mcConfig().AddHadronicId(3334);
  PHYSICS->mcConfig().AddHadronicId(4101);
  PHYSICS->mcConfig().AddHadronicId(4103);
  PHYSICS->mcConfig().AddHadronicId(4112);
  PHYSICS->mcConfig().AddHadronicId(4114);
  PHYSICS->mcConfig().AddHadronicId(4122);
  PHYSICS->mcConfig().AddHadronicId(4132);
  PHYSICS->mcConfig().AddHadronicId(4201);
  PHYSICS->mcConfig().AddHadronicId(4203);
  PHYSICS->mcConfig().AddHadronicId(4212);
  PHYSICS->mcConfig().AddHadronicId(4214);
  PHYSICS->mcConfig().AddHadronicId(4222);
  PHYSICS->mcConfig().AddHadronicId(4224);
  PHYSICS->mcConfig().AddHadronicId(4232);
  PHYSICS->mcConfig().AddHadronicId(4301);
  PHYSICS->mcConfig().AddHadronicId(4303);
  PHYSICS->mcConfig().AddHadronicId(4312);
  PHYSICS->mcConfig().AddHadronicId(4314);
  PHYSICS->mcConfig().AddHadronicId(4322);
  PHYSICS->mcConfig().AddHadronicId(4324);
  PHYSICS->mcConfig().AddHadronicId(4332);
  PHYSICS->mcConfig().AddHadronicId(4334);
  PHYSICS->mcConfig().AddHadronicId(4403);
  PHYSICS->mcConfig().AddHadronicId(4412);
  PHYSICS->mcConfig().AddHadronicId(4414);
  PHYSICS->mcConfig().AddHadronicId(4422);
  PHYSICS->mcConfig().AddHadronicId(4424);
  PHYSICS->mcConfig().AddHadronicId(4432);
  PHYSICS->mcConfig().AddHadronicId(4434);
  PHYSICS->mcConfig().AddHadronicId(4444);
  PHYSICS->mcConfig().AddHadronicId(5101);
  PHYSICS->mcConfig().AddHadronicId(5103);
  PHYSICS->mcConfig().AddHadronicId(5112);
  PHYSICS->mcConfig().AddHadronicId(5114);
  PHYSICS->mcConfig().AddHadronicId(5122);
  PHYSICS->mcConfig().AddHadronicId(5132);
  PHYSICS->mcConfig().AddHadronicId(5142);
  PHYSICS->mcConfig().AddHadronicId(5201);
  PHYSICS->mcConfig().AddHadronicId(5203);
  PHYSICS->mcConfig().AddHadronicId(5212);
  PHYSICS->mcConfig().AddHadronicId(5214);
  PHYSICS->mcConfig().AddHadronicId(5222);
  PHYSICS->mcConfig().AddHadronicId(5224);
  PHYSICS->mcConfig().AddHadronicId(5232);
  PHYSICS->mcConfig().AddHadronicId(5242);
  PHYSICS->mcConfig().AddHadronicId(5301);
  PHYSICS->mcConfig().AddHadronicId(5303);
  PHYSICS->mcConfig().AddHadronicId(5312);
  PHYSICS->mcConfig().AddHadronicId(5314);
  PHYSICS->mcConfig().AddHadronicId(5322);
  PHYSICS->mcConfig().AddHadronicId(5324);
  PHYSICS->mcConfig().AddHadronicId(5332);
  PHYSICS->mcConfig().AddHadronicId(5334);
  PHYSICS->mcConfig().AddHadronicId(5342);
  PHYSICS->mcConfig().AddHadronicId(5401);
  PHYSICS->mcConfig().AddHadronicId(5403);
  PHYSICS->mcConfig().AddHadronicId(5412);
  PHYSICS->mcConfig().AddHadronicId(5414);
  PHYSICS->mcConfig().AddHadronicId(5422);
  PHYSICS->mcConfig().AddHadronicId(5424);
  PHYSICS->mcConfig().AddHadronicId(5432);
  PHYSICS->mcConfig().AddHadronicId(5434);
  PHYSICS->mcConfig().AddHadronicId(5442);
  PHYSICS->mcConfig().AddHadronicId(5444);
  PHYSICS->mcConfig().AddHadronicId(5503);
  PHYSICS->mcConfig().AddHadronicId(5512);
  PHYSICS->mcConfig().AddHadronicId(5514);
  PHYSICS->mcConfig().AddHadronicId(5522);
  PHYSICS->mcConfig().AddHadronicId(5524);
  PHYSICS->mcConfig().AddHadronicId(5532);
  PHYSICS->mcConfig().AddHadronicId(5534);
  PHYSICS->mcConfig().AddHadronicId(5542);
  PHYSICS->mcConfig().AddHadronicId(5544);
  PHYSICS->mcConfig().AddHadronicId(5554);
  PHYSICS->mcConfig().AddHadronicId(10111);
  PHYSICS->mcConfig().AddHadronicId(10113);
  PHYSICS->mcConfig().AddHadronicId(10211);
  PHYSICS->mcConfig().AddHadronicId(10213);
  PHYSICS->mcConfig().AddHadronicId(10221);
  PHYSICS->mcConfig().AddHadronicId(10223);
  PHYSICS->mcConfig().AddHadronicId(10311);
  PHYSICS->mcConfig().AddHadronicId(10313);
  PHYSICS->mcConfig().AddHadronicId(10321);
  PHYSICS->mcConfig().AddHadronicId(10323);
  PHYSICS->mcConfig().AddHadronicId(10331);
  PHYSICS->mcConfig().AddHadronicId(10333);
  PHYSICS->mcConfig().AddHadronicId(10411);
  PHYSICS->mcConfig().AddHadronicId(10413);
  PHYSICS->mcConfig().AddHadronicId(10421);
  PHYSICS->mcConfig().AddHadronicId(10423);
  PHYSICS->mcConfig().AddHadronicId(10431);
  PHYSICS->mcConfig().AddHadronicId(10433);
  PHYSICS->mcConfig().AddHadronicId(10441);
  PHYSICS->mcConfig().AddHadronicId(10443);
  PHYSICS->mcConfig().AddHadronicId(10511);
  PHYSICS->mcConfig().AddHadronicId(10513);
  PHYSICS->mcConfig().AddHadronicId(10521);
  PHYSICS->mcConfig().AddHadronicId(10523);
  PHYSICS->mcConfig().AddHadronicId(10531);
  PHYSICS->mcConfig().AddHadronicId(10533);
  PHYSICS->mcConfig().AddHadronicId(10541);
  PHYSICS->mcConfig().AddHadronicId(10543);
  PHYSICS->mcConfig().AddHadronicId(10551);
  PHYSICS->mcConfig().AddHadronicId(10553);
  PHYSICS->mcConfig().AddHadronicId(20113);
  PHYSICS->mcConfig().AddHadronicId(20213);
  PHYSICS->mcConfig().AddHadronicId(20223);
  PHYSICS->mcConfig().AddHadronicId(20313);
  PHYSICS->mcConfig().AddHadronicId(20323);
  PHYSICS->mcConfig().AddHadronicId(20333);
  PHYSICS->mcConfig().AddHadronicId(20413);
  PHYSICS->mcConfig().AddHadronicId(20423);
  PHYSICS->mcConfig().AddHadronicId(20433);
  PHYSICS->mcConfig().AddHadronicId(20443);
  PHYSICS->mcConfig().AddHadronicId(20513);
  PHYSICS->mcConfig().AddHadronicId(20523);
  PHYSICS->mcConfig().AddHadronicId(20533);
  PHYSICS->mcConfig().AddHadronicId(20543);
  PHYSICS->mcConfig().AddHadronicId(20553);
  PHYSICS->mcConfig().AddHadronicId(100443);
  PHYSICS->mcConfig().AddHadronicId(100553);
  PHYSICS->mcConfig().AddHadronicId(9900440);
  PHYSICS->mcConfig().AddHadronicId(9900441);
  PHYSICS->mcConfig().AddHadronicId(9900443);
  PHYSICS->mcConfig().AddHadronicId(9900551);
  PHYSICS->mcConfig().AddHadronicId(9900553);
  PHYSICS->mcConfig().AddHadronicId(9910441);
  PHYSICS->mcConfig().AddHadronicId(9910551);

  // definition of the multiparticle "invisible"
  PHYSICS->mcConfig().AddInvisibleId(-16);
  PHYSICS->mcConfig().AddInvisibleId(-14);
  PHYSICS->mcConfig().AddInvisibleId(-12);
  PHYSICS->mcConfig().AddInvisibleId(12);
  PHYSICS->mcConfig().AddInvisibleId(14);
  PHYSICS->mcConfig().AddInvisibleId(16);
  PHYSICS->mcConfig().AddInvisibleId(1000022);
  PHYSICS->mcConfig().AddInvisibleId(1000039);

  // Initializing PhysicsService for RECO
  PHYSICS->recConfig().Reset();

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("myregion");

  // ===== Selections ===== //

  // ===== Histograms ===== //
/*  Manager()->AddHisto("NBQuarks", 12,0.0,11.0);
  Manager()->AddHisto("bquark pt", 50, 0.0,1000.0);
  Manager()->AddHisto("bb DelR", 30,0.0,6.0);
  Manager()->AddHisto("LHiggs pt", 50, 0.0, 1000.0);
  Manager()->AddHisto("LHiggs DelR", 30,0.0,6.0);*/
  Manager()->AddHisto("PassSelection", 2,0.0,1.0);
  Manager()->AddHisto("Nbjets", 6,1.0,6.0);
  Manager()->AddHisto("bb mass", 50,0.0,1700.0);
  Manager()->AddHisto("bjet pt all", 40, 0.0, 1100.0);
  Manager()->AddHisto("Two bjet mass", 40, 0.0, 1250.0);
  Manager()->AddHisto("bjet mass1", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet mass2", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet mass", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet DelR", 30,0.0,6.0);
  Manager()->AddHisto("bjet pt1", 40,0.0,1000.0);
  Manager()->AddHisto("bjet pt2", 40,0.0,1000.0);
  Manager()->AddHisto("bjet pt3", 40,0.0,1000.0);
  Manager()->AddHisto("bjet pt4", 40,0.0,1000.0);
  Manager()->AddHisto("bjet_mass1", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet_mass2", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet_mass3", 40, 0.0, 250.0);
  Manager()->AddHisto("bjet_mass4", 40, 0.0, 250.0);
  //Manager()->AddHisto("NBQuarks", 6,1.0,6.0);
  // No problem during initialization
  return true;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool ML_fatjetimages_delphes::Execute(SampleFormat& sample, const EventFormat& event)
{
  MAfloat32 __event_weight__ = 1.0;
  if (weighted_events_ && event.mc()!=0) __event_weight__ = event.mc()->weight();

  if (sample.mc()!=0) sample.mc()->addWeightedEvents(__event_weight__);
  Manager()->InitializeForNewEvent(__event_weight__);
/*
  // Parton level particles
  // Light Higgs container
  std::vector<const MCParticleFormat*> LHiggs;
  std::vector<const MCParticleFormat*> pair1;
  std::vector<const MCParticleFormat*> pair2;
  for(MAuint32 i=0;i<event.mc()->particles().size();i++)
  {
    const MCParticleFormat *ThisParticle = &event.mc()->particles()[i];
    if(abs(ThisParticle->pdgid()) == 25)
    {
      int n = 0; // to track higgs in mothers list
      for(MAuint32 j=0;j<ThisParticle->mothers().size();j++)
      {
        if(abs(ThisParticle->mothers()[j]->pdgid()) == 25) // avoid overcounting, check mothers for higgs
        {
          n++;
        }
      }
      if(n==0)
      {
        LHiggs.push_back(ThisParticle); // store higgs if no higgs' found in mothers
      }
      else continue;
    }
  }


  // b and c quarks for truth tagging
  //std::vector<const MCParticleFormat*> Cquarks;
  std::vector<const MCParticleFormat*> AllBQuarks;
  std::vector<const MCParticleFormat*> Bquarks;
  std::vector<const MCParticleFormat*> AntiBQuarks;
  std::vector<const MCParticleFormat*> BPair1;
  std::vector<const MCParticleFormat*> BPair2;
  for(MAuint32 i=0;i<event.mc()->particles().size();i++)
  {
    const MCParticleFormat *ThisParticle = &event.mc()->particles()[i];
    if(abs(ThisParticle->pdgid()) == 5)
    {
      int nb = 0; // to track b in mothers list
      int nh = 0; // to track higgs in b mothers list
      for(MAuint32 j=0;j<ThisParticle->mothers().size();j++)
      {
        if(abs(ThisParticle->mothers()[j]->pdgid()) == 5) // avoid overcounting, check mothers for b quarks
        {
          nb++;
        }
        if(abs(ThisParticle->mothers()[j]->pdgid()) == 25) // check mothers for higgs
        {
          nh++;
        }
      }
      if((nb==0) & (ThisParticle->pdgid() == 5))
      {
        AllBQuarks.push_back(ThisParticle); // store b quark if no b's in mothers
        Bquarks.push_back(ThisParticle);
      }
      if((nb==0) & (ThisParticle->pdgid() == -5))
      {
        AllBQuarks.push_back(ThisParticle); // store anti b quark if no b's in mothers
        AntiBQuarks.push_back(ThisParticle);
      }
      else continue;
    }
  }

  // Remove additional b quarks if >4
  // Sort BQuark container such that they're pt ordered
  std::vector<const MCParticleFormat*> SortedBQuarks;
  SORTER->sort(AllBQuarks, PTordering);

  if(AllBQuarks.size() > 4)
  {
    AllBQuarks.resize(4);
    //while(AllBQuarks.size() > 4)
    //{
    //  AllBQuarks.erase(AllBQuarks.begin());
    //}
  }

  // Gather b pairs together - messy as ttbar can contain 2,3 or 4 b quarks
  // 2 b quarks is trivial
  if(AllBQuarks.size() == 2)
  {
    BPair1.push_back(AllBQuarks[0]);
    BPair1.push_back(AllBQuarks[1]);
  }
  // for 3, we construct closest pair
  if(AllBQuarks.size() == 3)
  {
    double dist12 = AllBQuarks[0]->dr(AllBQuarks[1]);
    double dist13 = AllBQuarks[0]->dr(AllBQuarks[2]);
    double dist23 = AllBQuarks[1]->dr(AllBQuarks[2]);
    double mindist1 = min({dist12,dist13,dist23});
    if(mindist1 == dist12)
    {
      BPair1.push_back(AllBQuarks[0]);
      BPair1.push_back(AllBQuarks[1]);
    }
    if(mindist1 == dist13)
    {
      BPair1.push_back(AllBQuarks[0]);
      BPair1.push_back(AllBQuarks[2]);
    }
    if(mindist1 == dist23)
    {
      BPair1.push_back(AllBQuarks[1]);
      BPair1.push_back(AllBQuarks[2]);
    }
  }
  // for 4, we find closest partner for one, then pair other two by default
  if(AllBQuarks.size() == 4)
  {
    double dist12 = AllBQuarks[0]->dr(AllBQuarks[1]);
    double dist13 = AllBQuarks[0]->dr(AllBQuarks[2]);
    double dist14 = AllBQuarks[0]->dr(AllBQuarks[3]);
    double mindist2 = min({dist12,dist13,dist14});
    if(mindist2 == dist12)
    {
      BPair1.push_back(AllBQuarks[0]);
      BPair1.push_back(AllBQuarks[1]);
      BPair2.push_back(AllBQuarks[2]);
      BPair2.push_back(AllBQuarks[3]);
    }
    if(mindist2 == dist13)
    {
      BPair1.push_back(AllBQuarks[0]);
      BPair1.push_back(AllBQuarks[2]);
      BPair2.push_back(AllBQuarks[1]);
      BPair2.push_back(AllBQuarks[3]);
    }
    if(mindist2 == dist14)
    {
      BPair1.push_back(AllBQuarks[0]);
      BPair1.push_back(AllBQuarks[3]);
      BPair2.push_back(AllBQuarks[1]);
      BPair2.push_back(AllBQuarks[2]);
    }
  }

  for(MAuint32 i=0;i<AllBQuarks.size();i++)
  {
    const MCParticleFormat *ThisQuark = AllBQuarks[i];
    for(MAuint32 j=0;j<ThisQuark->mothers().size();j++)
    {
      if((ThisQuark->mothers()[j]) == (LHiggs[0]))
      {
        pair1.push_back(ThisQuark);
      }
      else if((ThisQuark->mothers()[j]) == (LHiggs[1]))
      {
        pair2.push_back(ThisQuark);
      }
    }
  }

*/
  // Filling track/towers for clustering
  std::vector<PseudoJet> ClusterParticles;
  // EFlow Tracks
  for(MAuint32 i=0;i<event.rec()->EFlowTracks().size();i++)
  {
    const RecTrackFormat* thisTrack = & event.rec()->EFlowTracks()[i];
    // Sanity check for kinematics
    //if (abs(thisTrack->eta()) > 2.5) continue;
    //if (thisTrack->pt() < 0.5) continue;
    ClusterParticles.push_back(PseudoJet(thisTrack->px(),
                                   thisTrack->py(),
                                   thisTrack->pz(),
                                   thisTrack->e()));
  }


  // EFlow Neutral Hadrons
  for(MAuint32 i=0;i<event.rec()->EFlowNeutralHadrons().size();i++)
  {
    const RecParticleFormat* thisHadron = & event.rec()->EFlowNeutralHadrons()[i];
    // Sanity check for kinematics
    //if (abs(thisHadron->eta()) > 2.5) continue;
    //if (thisHadron->pt() < 0.5) continue;
    ClusterParticles.push_back(PseudoJet(thisHadron->px(),
                                   thisHadron->py(),
                                   thisHadron->pz(),
                                   thisHadron->e()));
  }

  // EFlow Photons
  for(MAuint32 i=0;i<event.rec()->EFlowPhotons().size();i++)
  {
    const RecParticleFormat* thisPhoton = & event.rec()->EFlowPhotons()[i];

    // Sanity check for kinematics
    //if (abs(thisPhoton->eta()) > 2.5) continue;
    //if (thisPhoton->pt() < 0.5) continue;
    ClusterParticles.push_back(PseudoJet(thisPhoton->px(),
                                   thisPhoton->py(),
                                   thisPhoton->pz(),
                                   thisPhoton->e()));
  }
/* // Use these when you have acesss to jets directly after running delphes
  std::vector<PseudoJet> Raw_Jets;

  for (MAuint32 i=0;i<event.rec()->jets().size();i++)
  {
    Raw_Jets.push_back(PseudoJet(event.rec()->jets()[i].px(),
                                event.rec()->jets()[i].py(),
                                event.rec()->jets()[i].pz(),
                                event.rec()->jets()[i].e()
                                ));
  }

*/
/*

// NOTE: WHEN SELECTING VARR/ANTIKT, ALSO NEED TO CHANGE THE Reff DEFINITIONS IN THE
// TAGGING CODE BELOW, FOR BOTH B AND C TAGGING
// Comment out below block uf using a fixed cone

  // VARIABLE R CLUSTERING
  // Jet clustering using fastjet with variable R plugin
  // Set varR parameters and higgs mass
  double H_mass = 700.;
  double h_mass = 125.;
  double rho = 300.;
  double Rmin = 0.4;
  double Rmax = 2.0;

  VariableRPlugin lvjet_pluginAKT(rho, Rmin, Rmax, VariableRPlugin::AKTLIKE);
  JetDefinition jet_defAKT(&lvjet_pluginAKT);
  ClusterSequence clust_seqAKT(ClusterParticles, jet_defAKT);

  // Do clustering with chosen algorithm
  std::vector<PseudoJet> Raw_Jets = sorted_by_pt(clust_seqAKT.inclusive_jets()); // edit to choose which clustering to perform!
  std::vector<PseudoJet> Jets; // Jet container for after cuts have been applied

// End of block to comment out
*/
// Comment out block below if using varR
  // ANTI KT CLUSTERING
  double Reff=1.2;
  double H_mass = 700.;
  double h_mass = 125.;
  JetDefinition jet_def(antikt_algorithm, Reff); // jet definition
  ClusterSequence clust_seq(ClusterParticles, jet_def); // run clustering
  vector<fastjet::PseudoJet> Raw_Jets = sorted_by_pt(clust_seq.inclusive_jets());
  std::vector<PseudoJet> Jets;
// end of block to comment out

// Only keep jets with pt > 200GeV and pt < 400GeV

  for(MAuint32 i=0;i<Raw_Jets.size();i++)
  {
    if(Raw_Jets[i].pt() > 200.) Jets.push_back(Raw_Jets[i]);
  }

  // Grab B and C Quarks - use for tagging above jets
  std::vector<const MCParticleFormat*> BQuarks;
  std::vector<const MCParticleFormat*> CQuarks;
  for(MAuint32 i=0;i<event.mc()->particles().size();i++)
  {
    const MCParticleFormat *ThisParticle = &event.mc()->particles()[i];
    if(abs(ThisParticle->pdgid()) == 5)
    {
      int n = 0; // to track b in mothers list
      for(MAuint32 j=0;j<ThisParticle->mothers().size();j++)
      {
        if(abs(ThisParticle->mothers()[j]->pdgid()) == 5) // avoid overcounting, check mothers for b quarks
        {
          n++; // store b quark if no b's in momthers
        }
      }
      if(n==0)
      {
        BQuarks.push_back(ThisParticle);
      }
      else continue;
    }
    if(abs(ThisParticle->pdgid()) == 4)
    {
      int m = 0; // to track c in mothers list
      for(MAuint32 j=0;j<ThisParticle->mothers().size();j++)
      {
        if(abs(ThisParticle->mothers()[j]->pdgid()) == 4) // avoid overcounting, check mothers for c quarks
        {
          m++;
        }
      }
      if(m==0)
      {
        CQuarks.push_back(ThisParticle); // store c quark if no b's in mothers
      }
      else continue;
    }
  }

  //srand (time(NULL)); // seed for tagging efficiencies

  // BTagger using MCBQuarks
  std::vector<PseudoJet> BJets;
  //std::vector<PseudoJet> BCands;
  //vector<PseudoJet> BQCands;
  std::vector<const MCParticleFormat*> BQCands;
  for(MAuint32 i=0;i<Jets.size();i++)
  {
    for(MAuint32 j=0;j<BQuarks.size();j++)
    {
      //const MCParticleFormat *ThisBQuark = BQuarks[j];
      double phi_diff = abs(BQuarks[j]->phi() - Jets[i].phi());
      if(abs(phi_diff) > M_PI) phi_diff = 2*M_PI - phi_diff;
      double dist = sqrt(pow(BQuarks[j]->eta() - Jets[i].eta(),2.0) + pow(phi_diff,2.0));
      //double Reff = rho / Jets[i].pt(); //use for varr
      // Set Reff value

      if (rho / Jets[i].pt() < Rmin)
        {
          Reff = Rmin;
        }
      else if (rho / Jets[i].pt() > Rmax)
        {
          Reff = Rmax;
        }
      else Reff = rho / Jets[i].pt();

      // End of block to comment out
      if(dist <= Reff)
      {
        BQCands.push_back(BQuarks[j]); //bquarks within Reff are listed as bquarks candidates for fatjets
        continue;
      }
      if(dist > Reff)
      {
        continue;
      }
    }

    if(BQCands.size() == 2) // if exactly 2 bquarks were found for the fatjet candidate then we can keep the jet as double b-tagged fatjet
    {
      // Apply btagging efficiency
      //double btag_rand = rand() % 101;
      //double btag_eff = 0.85*tanh(0.0025*Jets[i].pt())*(25.0/(1+0.063*Jets[i].pt()));
      //if(btag_rand > 100-(100*btag_eff)) BJets.push_back(Jets[i]);
      BJets.push_back(Jets[i]);
      for (MAuint32 k=0;k<BQCands.size();k++)
      {
        for(MAuint32 j=0;j<BQuarks.size();j++)
        {
          if (BQCands[k] == BQuarks[j]) BQuarks.erase(BQuarks.begin() + j);
        }
      }
    }
    BQCands.clear();
  }
/*
  // CTagger using MCBQuarks
  std::vector<PseudoJet> CJets;
  std::vector<const MCParticleFormat*> CQCands;
  for(MAuint32 i=0;i<Jets.size();i++)
  {
    for(MAuint32 j=0;j<CQuarks.size();j++)
    {
      //const MCParticleFormat *ThisBQuark = BQuarks[j];
      double phi_diff = abs(CQuarks[j]->phi() - Jets[i].phi());
      if(abs(phi_diff) > M_PI) phi_diff = 2*M_PI - phi_diff;
      double dist = sqrt(pow(CQuarks[j]->eta() - Jets[i].eta(),2.0) + pow(phi_diff,2.0));
      double Reff = rho / Jets[i].pt();
      // Set Reff value

      if (rho / Jets[i].pt() < Rmin)
        {
          Reff = Rmin;
        }
      else if (rho / Jets[i].pt() > Rmax)
        {
          Reff = Rmax;
        }
      else Reff = rho / Jets[i].pt();

      // End of block to comment out
      if(dist <= Reff)
      {
        CQCands.push_back(CQuarks[j]); //bquarks within Reff are listed as bquarks candidates for fatjets
        continue;
      }
      if(dist > Reff)
      {
        continue;
      }
    }

    if(CQCands.size() == 2) // if exactly 2 bquarks were found for the fatjet candidate then we can keep the jet as double b-tagged fatjet
    {
      // Apply btagging efficiency
      //double c_rand = rand() % 101;
      //double c_bmistag = 0.25*tanh(0.018*Jets[i].pt())*(1/(1+ 0.0013*Jets[i].pt()));
      //if(c_rand > 100-(100*c_bmistag)) BJets.push_back(Jets[i]);
      else CJets.push_back(Jets[i]);
      for (MAuint32 k=0;k<CQCands.size();k++)
      {
        for(MAuint32 j=0;j<CQuarks.size();j++)
        {
          if (CQCands[k] == CQuarks[j]) CQuarks.erase(CQuarks.begin() + j);
        }
      }
    }
    CQCands.clear();
  }


  for(MAuint32 i=0;i<Jets.size();i++)
  {
    double l_rand = rand() % 101;
    double l_bmistag = 0.01+0.000038*Jets[i].pt();
    if(l_rand > 100-(100*l_bmistag)) BJets.push_back(Jets[i]);
  }

 */

  std::vector<PseudoJet> SortedBJets = sorted_by_pt(BJets);
  //std::vector<PseudoJet> SortedJets = sorted_by_pt(Raw_Jets);
  //if(SortedJets.size() > 2) SortedJets.resize(2);

  // Count number of bjets following pt cut to leading btagged jet
  int numbjets = 0;

  std::vector<PseudoJet> Leading_BJet, SubLeading_BJet, Sub2Leading_BJet, Sub3Leading_BJet;

  // Filling containers based on number of bjets found
  if(BJets.size() > 0)
  {
    Leading_BJet.push_back(SortedBJets[0]);
    numbjets ++;

    if(BJets.size() >= 2)
    {
      SubLeading_BJet.push_back(SortedBJets[1]);
      numbjets++;
    }

    if(BJets.size() >= 3)
    {
      Sub2Leading_BJet.push_back(SortedBJets[2]);
      numbjets++;
    }

    if(BJets.size() >= 4)
    {
      Sub3Leading_BJet.push_back(SortedBJets[3]);
      numbjets++;
    }

  }

  std::vector<PseudoJet> Bpair1;

  if(SortedBJets.size() == 2)
  {
    Bpair1.push_back(SortedBJets[0]);
    Bpair1.push_back(SortedBJets[1]);
  }

  PseudoJet Dijets;

  if(SortedBJets.size() == 2)
  {
    Dijets = PseudoJet(SortedBJets[0].px()+SortedBJets[1].px(),
                                  SortedBJets[0].py()+SortedBJets[1].py(),
                                  SortedBJets[0].pz()+SortedBJets[1].pz(),
                                  SortedBJets[0].e()+SortedBJets[1].e()
                                  );
  }

  // Histogram declarations
  //   * Plot: PassSelection
  {
    if(SortedBJets.size() == 2 && SortedBJets[0].m() > 100. && SortedBJets[0].m() < 150. && SortedBJets[1].m() > 100. && SortedBJets[1].m() < 150. &&// select events with 4 b jets, and m_bbbb close to m_H
      Dijets.m() > 650. && Dijets.m() < 750.)
    {
      Manager()->FillHisto("PassSelection", 0.9);
    }
    else Manager()->FillHisto("PassSelection",0);
  }
/*
  {
    if(SortedJets.size() == 2 && SortedJets[0].m() > 100. && SortedJets[0].m() < 150. && SortedJets[1].m() > 100. && SortedJets[1].m() < 150. &&// select events with 4 b jets, and m_bbbb close to m_H
      JDijets.m() > 650. && JDijets.m() < 750.)
    {
      Manager()->FillHisto("PassSelection1", 0.9);
    }
    else Manager()->FillHisto("PassSelection1",0);
  }


  // BQuark pt
  {
    for(MAuint32 i=0;i<AllBQuarks.size();i++)
    {
      Manager()->FillHisto("bquark pt", AllBQuarks[i]->pt());
    }
  }

  // bb DelR
  {
    if(BPair1.size() > 0)
    {
      MAfloat32 dist1 = BPair1[0]->dr(BPair1[1]);
      Manager()->FillHisto("bb DelR", dist1);
    }
    if(BPair2.size() > 0)
    {
      MAfloat32 dist2 = BPair2[0]->dr(BPair2[1]);
      Manager()->FillHisto("bb DelR", dist2);
    }
  }

  // LHiggs pt
  {
    for(MAuint32 i=0;i<LHiggs.size();i++)
    {
      if(LHiggs.size() > 0) Manager()->FillHisto("LHiggs pt", LHiggs[i]->pt());
    }
    if(LHiggs.size() == 0) Manager()->FillHisto("LHiggs pt", 0.);
  }

  // LHiggs DelR
  {
    if(LHiggs.size() > 0)
    {
      MAfloat32 dist = LHiggs[0]->dr(LHiggs[1]);
      Manager()->FillHisto("LHiggs DelR", dist);
    }
    if(LHiggs.size() == 0) Manager()->FillHisto("LHiggs pt", 0.);
  }

*/
  //   * Plot: Nbjets
  {
    if(SortedBJets.size() > 0) // don't plot zeros (makes histos look nicer)
    {
      Manager()->FillHisto("Nbjets", SortedBJets.size());
    }
  }
/*
  {
    if(SortedJets.size() > 0) // don't plot zeros (makes histos look nicer)
    {
      Manager()->FillHisto("Njets", SortedJets.size());
    }
  }
*/
  //   * Plot: bb mass
  {
    if(SortedBJets.size() == 2)
    {
      PseudoJet bdijet = PseudoJet(SortedBJets[0].px()+SortedBJets[1].px(),
                                    SortedBJets[0].py()+SortedBJets[1].py(),
                                    SortedBJets[0].pz()+SortedBJets[1].pz(),
                                    SortedBJets[0].e()+SortedBJets[1].e()
                                    );
      Manager()->FillHisto("bb mass", bdijet.m());
    }
  }
/*
  {
    if(SortedJets.size() == 2)
    {
      PseudoJet Jbdijet = PseudoJet(SortedJets[0].px()+SortedJets[1].px(),
                                    SortedJets[0].py()+SortedJets[1].py(),
                                    SortedJets[0].pz()+SortedJets[1].pz(),
                                    SortedJets[0].e()+SortedJets[1].e()
                                    );
      Manager()->FillHisto("jj mass", Jbdijet.m());
    }
  }
*/
  //BJet mass
  {
    for(MAuint32 j=0;j<SortedBJets.size();j++)
    {
      Manager()->FillHisto("bjet mass", SortedBJets[j].m());
    }
  }
/*
  {
    for(MAuint32 j=0;j<SortedJets.size();j++)
    {
      Manager()->FillHisto("jet mass", SortedJets[j].m());
    }
  }
*/
  //   * Plot: bjet pt all
  {
    for(MAuint32 j=0;j<SortedBJets.size();j++)
    {
      Manager()->FillHisto("bjet pt all", SortedBJets[j].pt());
    }
  }
/*
  {
    for(MAuint32 j=0;j<SortedJets.size();j++)
    {
      Manager()->FillHisto("jet pt all", SortedJets[j].pt());
    }
  }
*/
  // *Plot:twobjet mass
  {
    if(SortedBJets.size() == 2)
    {
      for(MAuint32 j=0;j<SortedBJets.size();j++)
      {
        Manager()->FillHisto("Two bjet mass", SortedBJets[j].m());
      }
    }
  }

  // *plot:bjet1 mass1
  {
    if(SortedBJets.size() == 2)
    {
      for(MAuint32 j=0;j<SortedBJets.size();j++)
      {
        Manager()->FillHisto("bjet mass1", SortedBJets[0].m());
        Manager()->FillHisto("bjet mass2", SortedBJets[1].m());
      }
    }
  }
/*
  {
    if(SortedJets.size() == 2)
    {
      for(MAuint32 j=0;j<SortedJets.size();j++)
      {
        Manager()->FillHisto("jet mass1", SortedJets[0].m());
        Manager()->FillHisto("jet mass2", SortedJets[1].m());
      }
    }
  }
*/
  // bjet DelR
  {
    if(Bpair1.size() > 0)
    {
      double bphi_diff = abs( Bpair1[0].phi()- Bpair1[1].phi());
      if(abs(bphi_diff) > M_PI) bphi_diff = 2*M_PI - bphi_diff;
      double bdist = sqrt(pow(Bpair1[0].eta() - Bpair1[1].eta(),2.0) + pow(bphi_diff,2.0));
      Manager()->FillHisto("bjet DelR", bdist);
    }
  }

/*
  // NLHiggs
  {
    Manager()->FillHisto("NLHiggs", LHiggs.size());
  }

  {
    Manager()->FillHisto("pair1", pair1.size());
  }

  {
    Manager()->FillHisto("pair2", pair2.size());
  }


  //bb from same LHiggs DelR
  {
    if(pair1.size() > 0)
    {
      MAfloat32 dist1 = pair1[0]->dr(pair1[1]);
      Manager()->FillHisto("bb from same LHiggs DelR1", dist1);
    }
    if(pair2.size() > 0)
    {
      MAfloat32 dist2 = pair2[0]->dr(pair2[1]);
      Manager()->FillHisto("bb from same LHiggs DelR2", dist2);
    }
  }
*/

  //   * Plot: bjet pt1
  {
    if(Leading_BJet.size() > 0)
    {
      Manager()->FillHisto("bjet pt1", Leading_BJet[0].pt());
    }
  }

  //   * Plot: bjet pt2
  {
    if(SubLeading_BJet.size() > 0)
    {
      Manager()->FillHisto("bjet pt2", SubLeading_BJet[0].pt());
    }
  }

  //   * Plot: bjet pt3
  {
    if(Sub2Leading_BJet.size() > 0)
    {
      Manager()->FillHisto("bjet pt3", Sub2Leading_BJet[0].pt());
    }
  }

  //   * Plot: bjet pt4
    {
      if(Sub3Leading_BJet.size() > 0)
      {
        Manager()->FillHisto("bjet pt4", Sub3Leading_BJet[0].pt());
      }
    }

    //   * Plot: bjet pt1
    {
      if(Leading_BJet.size() > 0)
      {
        Manager()->FillHisto("bjet_mass1", Leading_BJet[0].m());
      }
    }

    //   * Plot: bjet pt2
    {
      if(SubLeading_BJet.size() > 0)
      {
        Manager()->FillHisto("bjet_mass2", SubLeading_BJet[0].m());
      }
    }

    //   * Plot: bjet pt3
    {
      if(Sub2Leading_BJet.size() > 0)
      {
        Manager()->FillHisto("bjet_mass3", Sub2Leading_BJet[0].m());
      }
    }

    //   * Plot: bjet pt4
      {
        if(Sub3Leading_BJet.size() > 0)
        {
          Manager()->FillHisto("bjet_mass4", Sub3Leading_BJet[0].m());
        }
      }



  // Only store events for leading fat b tagged jets
    if(Leading_BJet.size() > 0)
    {
      // reclustering the constituents of leading fat b-jet using anti-kt algo with R=0.2
      double R = 0.4;
      JetDefinition jet_def(kt_algorithm, R); // jet DEFINITIONS
      ClusterSequence cs(Leading_BJet[0].constituents(), jet_def); // run CLUSTERING
      vector<fastjet::PseudoJet> sJets = sorted_by_pt(cs.inclusive_jets());
      vector<fastjet::PseudoJet> sub_Jets;
      for(MAuint32 i=0;i<sJets.size();i++)
      {
        if(sJets[i].pt() > 0.05*Leading_BJet[0].pt())
           sub_Jets.push_back(sJets[i]);
      }

      ofstream myfile1;
      myfile1.open("delphesfatjet_ma5_signal_doublebtag_subjet04_new_constituents_withisr_mpi_r12_100k_info.csv",ios::app);
      for(MAuint32 i=0;i<sub_Jets.size();i++)
      {
        for(MAuint32 j=0;j<sub_Jets[i].constituents().size();j++)
        {
          myfile1 << i+1 << "," << sub_Jets[i].pt() << "," << sub_Jets[i].m() << "," << sub_Jets[i].eta() << "," << sub_Jets[i].phi() << "," <<
                 j+1 << "," << sub_Jets[i].constituents()[j].user_index() << "," << sub_Jets[i].constituents()[j].pt() << "," <<
                 sub_Jets[i].constituents()[j].m() << "," << sub_Jets[i].constituents()[j].eta() << "," << sub_Jets[i].constituents()[j].phi() << "\n";
        }
      }
      myfile1 << "0,0,0,0,0,0,0,0,0,0,0\n"; // line of zeroes to represent new event
      myfile1.close();


      // Write constituent information to text file
      ofstream myfile2;
      myfile2.open("delphesfatjet_ma5_signal_doublebtag_leadingbjet_new_constituent_withisr_mpi_r12_100k_info.csv",ios::app);
      for(MAuint32 i=0;i<Leading_BJet[0].constituents().size();i++)
      {
        myfile2 << i+1 << "," << Leading_BJet[0].constituents()[i].user_index() << "," << Leading_BJet[0].constituents()[i].pt() << "," <<
                  Leading_BJet[0].constituents()[i].m() << "," << Leading_BJet[0].constituents()[i].eta() << "," << Leading_BJet[0].constituents()[i].phi()<< "\n";
      }

      myfile2 << "0,0,0,0,0,0\n"; // line of zeroes to represent new event
      myfile2.close();
      
    }

  return true;

}

void ML_fatjetimages_delphes::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
}
