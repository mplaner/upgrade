// Package:    Ntuplizer
// Class:      Ntuplizer
// 
/**\class Ntuplizer Ntuplizer.cc Ntuple/Ntuplizer/src/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael David Planer
//         Created:  Fri Jul 19 03:51:27 EDT 2013
// $Id$
//
//


// system include files
#include <memory>

//for conversion vertex truth
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include <exception>

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/JetReco/interface/GenJet.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerSimple.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"

#include <stdint.h>
#include "TTree.h"
#include "TBranch.h"

using namespace std;
using namespace edm;

#define MAX_PHOTON 2
#define MAX_GEN 2
//
// class declaration
//

class Diphoton : public edm::EDAnalyzer {
public:
  explicit Diphoton(const edm::ParameterSet&);
  ~Diphoton();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob(const edm::EventSetup&, const edm::ParameterSet&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::EDGetTokenT<View<pat::Photon> > photonToken_;
  edm::EDGetTokenT<View<pat::PackedGenParticle> > genToken_;
  //  edm::EDGetTokenT<View<reco::Photon> > photonToken_;
  //edm::EDGetTokenT<View<reco::GenParticle> > genToken_;
  //  edm::Handle<std::vector<SimTrack>> simTrack;
  edm::EDGetTokenT<View<reco::Vertex> > vertexToken_;
  edm::EDGetTokenT<View<SimTrack> >  simTrackToken_;
  edm::EDGetTokenT<View<SimVertex> > simVertexToken_;
  edm::EDGetTokenT<std::vector<PCaloHit> > simCaloHitToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalRecHitsToken_;

  //  edm::InputTag  g4_simTk_Token_;
  //  edm::InputTag g4_simVtx_Token_;       
  
  edm::Service<TFileService> fs;

  TTree *t_photon;


  //photon branches
  TBranch *b_p_e12;
  TBranch *b_p_e15;
  TBranch *b_p_eAll;
  TBranch *b_p1_e12;
  TBranch *b_p1_e15;
  TBranch *b_p1_eAll;
  TBranch *b_p_r9;
  TBranch *b_p_full5x5r9;
  TBranch *b_p_ecal_energy;
  TBranch *b_p_photon_energy;
  TBranch *b_p_dr;
  TBranch *b_p_eta;
  TBranch *b_p_phi;
  TBranch *b_p_SCeta;
  TBranch *b_p_SCphi;
  TBranch *b_p_energy_1x5;
  TBranch *b_p_energy_2x5;
  TBranch *b_p_energy_3x3;
  TBranch *b_p_energy_full3x3;
  TBranch *b_p_energy_5x5;
  TBranch *b_p_energy_full5x5;
  TBranch *b_p_raw_energy;
  TBranch *b_p_raw_pt;
  TBranch *b_p_es_energy;
  TBranch *b_p_pt;
  TBranch *b_p_nRecHits;
  TBranch *b_p_eRecHits;
  TBranch *b_p_size;
  TBranch *b_p_isGap;
  TBranch *b_p_sigmaIetaIeta;
  TBranch *b_p_hOverE;
  TBranch *b_p_ecalRecHitSumEtConeDR03;
  TBranch *b_p_hcalTowerSumEtConeDR03;
  TBranch *b_p_trkSumPtSolidConeDR03;
  TBranch *b_p_neutralIsoConeDR03;
  TBranch *b_p_photonIsoConeDR03;
  TBranch *b_p_chargedIsoConeDR03;

  TBranch *b_p_mass;
  TBranch *b_p_rawmass;
  TBranch *b_p_genVrawmass;
  TBranch *b_p_genErawmass;
  TBranch *b_p_5x5mass;
  TBranch *b_p_full5x5mass;
  TBranch *b_p_3x3mass;
  TBranch *b_p_full3x3mass;

  
  float p_e12[MAX_PHOTON];
  float p_e15[MAX_PHOTON];
  float p_eAll[MAX_PHOTON];
  float p1_e12[MAX_PHOTON];
  float p1_e15[MAX_PHOTON];
  float p1_eAll[MAX_PHOTON];
  
  float p_photonenergy[MAX_PHOTON];
  float p_dr[MAX_PHOTON];
  double p_esenergy[MAX_PHOTON];
  double p_rawenergy[MAX_PHOTON];
  float p_rawpt[MAX_PHOTON];
  float p_ecalenergy[MAX_PHOTON];
  float p_energy_1x5[MAX_PHOTON];
  float p_energy_2x5[MAX_PHOTON];
  float p_energy_3x3[MAX_PHOTON];
  float p_energy_full3x3[MAX_PHOTON];
  float p_energy_5x5[MAX_PHOTON];
  float p_energy_full5x5[MAX_PHOTON];
  float p_pt_3x3[MAX_PHOTON];
  float p_pt_full3x3[MAX_PHOTON];
  float p_pt_5x5[MAX_PHOTON];
  float p_pt_full5x5[MAX_PHOTON];
  float p_pt[MAX_PHOTON];
  float p_nRecHits[MAX_PHOTON];
  float p_r9[MAX_PHOTON];
  float p_full5x5r9[MAX_PHOTON];
  float p_eRecHits[85][360][2];
  int   p_RecHitInCluster[85][360][2];
  
  bool p_isGap[MAX_PHOTON];
  float p_sigmaIetaIeta[MAX_PHOTON];
  float p_hOverE[MAX_PHOTON];
  float p_ecalRecHitSumEtConeDR03[MAX_PHOTON];
  float p_hcalTowerSumEtConeDR03[MAX_PHOTON];
  float p_trkSumPtSolidConeDR03[MAX_PHOTON];

  float p_neutralIsoConeDR03[MAX_PHOTON];
  float p_photonIsoConeDR03[MAX_PHOTON];
  float p_chargedIsoConeDR03[MAX_PHOTON];
  float p_mass;
  float p_rawmass;
  float p_genErawmass;
  float p_genVrawmass;
  float p_5x5mass;
  float p_full5x5mass;
  float p_3x3mass;
  float p_full3x3mass;

  
  double p_eta[MAX_PHOTON];
  double p_phi[MAX_PHOTON];
  double p_SCeta[MAX_PHOTON];
  double p_SCphi[MAX_PHOTON];
  int    p_size;



  //generator level branches
  TBranch *b_g_isConv;
  TBranch *b_g_calHits;
  TBranch *b_g_status;
  TBranch *b_g_ePt;
  TBranch *b_g_pPt;  
  TBranch* b_g_hdaughter;
  TBranch* b_g_zdaughter;
  TBranch* b_g_pt;
  TBranch* b_g_eta;
  TBranch *b_g_phi;
  TBranch *b_g_pdgid;
  TBranch *b_g_size;
  TBranch *b_g_mass;
  TBranch *b_g_run;
  TBranch *b_g_lumi;
  TBranch *b_g_event;
  
  //pileup branches
  TBranch *b_pu_size;
  TBranch *b_pu_bx;
  TBranch *b_vtx_size;
  TBranch *b_vtx_dz;
  
  float  g_calHits[85][360][2];
  int    g_pdgid[MAX_GEN];
  int    g_isConv[MAX_GEN];
  int    g_status[MAX_GEN];
  float  g_pPt[MAX_GEN];
  float  g_ePt[MAX_GEN];
  bool   g_hdaughter[MAX_GEN];
  bool   g_zdaughter[MAX_GEN];
  float  g_pt[MAX_GEN];
  double g_eta[MAX_GEN];
  double g_phi[MAX_GEN];
  float  g_mass;
  int    g_size;
  
  int    pu_size[16]; //number of vertices per bx
  int    pu_bx[16];   //bx index
  int    vtx_size;
  float    vtx_dz;
  
  int n_event;
  int n_lumi;
  int n_run;

};

Diphoton::Diphoton(const edm::ParameterSet& iConfig):
  photonToken_( consumes<View<pat::Photon> >( iConfig.getParameter<edm::InputTag> ( "photonTag" ) ) ),
  genToken_( consumes<View<pat::PackedGenParticle> >( iConfig.getParameter<edm::InputTag>( "genTag" ) ) ),
  //photonToken_( consumes<View<reco::Photon> >( iConfig.getParameter<edm::InputTag> ( "photonTag" ) ) ),
  //genToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>( "genTag" ) ) ),
  vertexToken_(    consumes<View<reco::Vertex> >( iConfig.getParameter<edm::InputTag> ( "vertexTag" ) ) ),
  simTrackToken_(  consumes<View<SimTrack>     >( iConfig.getParameter<edm::InputTag> ( "simTrackTag" ) ) ),
  simVertexToken_( consumes<View<SimVertex>    >( iConfig.getParameter<edm::InputTag> ( "simVertexTag" ) ) ),
  simCaloHitToken_( consumes<std::vector<PCaloHit>>( iConfig.getParameter<edm::InputTag> ( "simCaloHitTag" ) ) ),
  ecalRecHitsToken_( consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalRecHits") ) )
{
  t_photon = fs->make<TTree>("photon","photon");
  
  /*********branch initialization for pileup******/
  b_pu_size          = t_photon->Branch("pusize",pu_size, "pusize[16]/I");
  b_pu_bx            = t_photon->Branch("puBx",pu_bx,"puBx[16]/I");
  b_vtx_size            = t_photon->Branch("vtxsize",&vtx_size,"vtxsize/I");
  b_vtx_dz            = t_photon->Branch("vtxdz",&vtx_dz,"vtxdz/F");
  
/*****branch initialization for photons****/
  b_p_size          = t_photon->Branch("psize",&p_size, "psize/I");
  b_p_nRecHits      = t_photon->Branch("pNRecHits"       ,p_nRecHits,       "pNRecHits[psize]/F");
  b_p_5x5mass       = t_photon->Branch("p5x5mass",&p_5x5mass, "p5x5mass/F");
  b_p_full5x5mass       = t_photon->Branch("pfull5x5mass",&p_full5x5mass, "pfull5x5mass/F");
  b_p_3x3mass       = t_photon->Branch("p3x3mass",&p_3x3mass, "p3x3mass/F");
  b_p_full3x3mass       = t_photon->Branch("pfull3x3mass",&p_full3x3mass, "pfull3x3mass/F");
  b_p_dr            = t_photon->Branch("pdr",      p_dr,      "pdr[psize]/F");
  b_p_e12           = t_photon->Branch("p_e12", p_e12,"p_e12[psize]/F");
  b_p_e15           = t_photon->Branch("p_e15", p_e15,"p_e15[psize]/F");
  b_p_eAll          = t_photon->Branch("p_eAll",p_eAll, "p_eAll[psize]/F");
  b_p1_e12           = t_photon->Branch("p1_e12", p1_e12,"p1_e12[psize]/F");
  b_p1_e15           = t_photon->Branch("p1_e15", p1_e15,"p1_e15[psize]/F");
  b_p1_eAll          = t_photon->Branch("p1_eAll",p1_eAll, "p1_eAll[psize]/F");
  b_p_photon_energy   = t_photon->Branch("pPhoton_energy", p_photonenergy,"pPhoton_energy[psize]/F");
  b_p_ecal_energy   = t_photon->Branch("pSC_energy", p_ecalenergy,"pSC_energy[psize]/F");
  b_p_eta           = t_photon->Branch("pEta",p_eta,"pEta[psize]/D");
  b_p_phi           = t_photon->Branch("pPhi",p_phi,"pPhi[psize]/D");
  b_p_SCeta           = t_photon->Branch("pSCeta",p_SCeta,"pSCeta[psize]/D");
  b_p_SCphi           = t_photon->Branch("pSCphi",p_SCphi,"pSCphi[psize]/D");
  b_p_es_energy     = t_photon->Branch("pEsenergy", p_esenergy, "pEsenergy[psize]/D");
  b_p_raw_energy    = t_photon->Branch("pRawenergy",p_rawenergy,"pRawenergy[psize]/D");
  b_p_raw_energy    = t_photon->Branch("pRawpt",p_rawpt,"pRawpt[psize]/F");
  b_p_pt            = t_photon->Branch("pPt"       ,p_pt,       "pPt[psize]/F");
  b_p_mass          = t_photon->Branch("pmass",&p_mass, "pmass/F");
  b_p_genVrawmass   = t_photon->Branch("pgenVrawmass",&p_genVrawmass, "pgenVrawmass/F");
  b_p_genErawmass   = t_photon->Branch("pgenErawmass",&p_genErawmass, "pgenErawmass/F");
  b_p_rawmass       = t_photon->Branch("prawmass",&p_rawmass, "prawmass/F");
  b_p_energy_1x5    = t_photon->Branch("p1x5_energy",p_energy_1x5,"p1x5_energy[psize]/F");
  b_p_energy_2x5    = t_photon->Branch("p2x5_energy",p_energy_2x5,"p2x5_energy[psize]/F");
  b_p_energy_3x3    = t_photon->Branch("p3x3_energy",p_energy_3x3,"p3x3_energy[psize]/F");
  b_p_energy_full3x3    = t_photon->Branch("pfull3x3_energy",p_energy_full3x3,"pfull3x3_energy[psize]/F");
  b_p_energy_5x5    = t_photon->Branch("p5x5_energy",p_energy_5x5,"p5x5_energy[psize]/F");
  b_p_energy_full5x5    = t_photon->Branch("pfull5x5_energy",p_energy_full5x5,"p5x5_fullenergy[psize]/F");
  b_p_r9            = t_photon->Branch("pR9",p_r9,"pR9[psize]/F");
  b_p_full5x5r9            = t_photon->Branch("pfull5x5R9",p_full5x5r9,"pfull5x5R9[psize]/F");
  b_p_isGap         = t_photon->Branch("pGap",p_isGap,"pGap[psize]/O");
  b_p_sigmaIetaIeta = t_photon->Branch("pSigmaIetaIeta",p_sigmaIetaIeta,"pSigmaIetaIeta[psize]/F");
  b_p_hOverE        =   t_photon->Branch("pHoverE",p_hOverE,"pHoverE[psize]/F");
  b_p_ecalRecHitSumEtConeDR03 = t_photon->Branch("p_ecalRecHitSumEtConeDR03",p_ecalRecHitSumEtConeDR03,"p_ecalRecHitSumEtConeDR03[psize]/F");
  b_p_hcalTowerSumEtConeDR03  = t_photon->Branch("p_hcalTowerSumEtConeDR03" ,p_hcalTowerSumEtConeDR03 ,"p_hcalTowerSumEtConeDR03[psize]/F");
  b_p_trkSumPtSolidConeDR03=t_photon->Branch("p_trkSumPtSolidConeDR03",p_trkSumPtSolidConeDR03,"p_trkSumPtSolidConeDR03[psize]/F");
  b_p_photonIsoConeDR03 = t_photon->Branch("p_photonIsoConeDR03",p_photonIsoConeDR03,"p_photonIsoConeDR03[psize]/F");
  b_p_neutralIsoConeDR03  = t_photon->Branch("p_neutralIsoConeDR03" ,p_neutralIsoConeDR03 ,"p_neutralIsoConeDR03[psize]/F");
  b_p_chargedIsoConeDR03=t_photon->Branch("p_chargedIsoConeDR03",p_chargedIsoConeDR03,"p_chargedIsoConeDR03[psize]/F");
  
  /*****branch initialization for gen particles ****/  
  b_g_size =       t_photon->Branch("gsize",&g_size,"gsize/I");
  b_g_mass =       t_photon->Branch("gmass",&g_mass,"gmass/F");
  b_g_pt =         t_photon->Branch("gPt", g_pt,"gPt[gsize]/F");
  b_g_isConv         = t_photon->Branch("gConv",g_isConv,"gConv[gsize]/I");
  b_g_status         = t_photon->Branch("gStatus",g_status,"gStatus[gsize]/I");
  b_g_ePt         = t_photon->Branch("gePt",g_ePt,"gePt[gsize]/F");
  b_g_pPt         = t_photon->Branch("gpPt",g_pPt,"gpPt[gsize]/F");
  b_g_hdaughter =  t_photon->Branch("gHdaughter", g_hdaughter,"gHdaughter[gsize]/O");
  b_g_zdaughter =  t_photon->Branch("gZdaughter", g_zdaughter,"gZdaughter[gsize]/O");
  b_g_eta =        t_photon->Branch("gEta",g_eta,"gEta[gsize]/D");
  b_g_phi =        t_photon->Branch("gPhi",g_phi,"gPhi[gsize]/D");  
  b_g_pdgid =      t_photon->Branch("pdgid",g_pdgid,"pdgid[gsize]/I");
  b_g_run =        t_photon->Branch("run_number",&n_run,"run_number/I");
  b_g_lumi =       t_photon->Branch("lumi_number",&n_lumi,"lumi_number/I");
  b_g_event =      t_photon->Branch("event_number",&n_event,"event_number/I");
}


Diphoton::~Diphoton()
{
}

void
Diphoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

      /*********for conv truth***********/
   /*
   //   std::cout << "in ntuplizer " << std::endl;
   
   //edm::Handle<std::vector<SimTrack>> simTrack;
   //   Handle<View<SimTrack>> simTracks;
   Handle<SimTrackContainer> simTracks;
   if(!iEvent.getByToken(simTrackToken_,simTracks))
     {
       std::cout << "simTracks not found" << std::endl;
       return;
     }
   
   //edm::Handle<std::vector<SimVertex>> simVertex;
   //Handle<View<SimVertex>> simVertices;
   Handle<SimVertexContainer> simVertices;
   if(!iEvent.getByToken(simVertexToken_,simVertices))
     {
       std::cout << "simVertices not found" << std::endl;
       return;
     }
   std::vector<SimTrack> theSimTracks;
   std::vector<SimVertex> theSimVertices;
   
   theSimTracks.insert(    theSimTracks.end(), simTracks->begin(), simTracks->end());
   theSimVertices.insert(theSimVertices.end(),simVertices->begin(),simVerices->end());

   PhotonMCTruthFinder*  thePhotonMCTruthFinder_;
   thePhotonMCTruthFinder_ = new PhotonMCTruthFinder();
   //std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find(&simTracks,  &simVertices);
   std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find(theSimTracks,  theSimVertices);
   */
   /*********end conv truth***********/
   Handle<View<SimTrack>> simTracks;
   iEvent.getByToken(simTrackToken_ ,simTracks);
   Handle<View<SimVertex>> simVertices;
   iEvent.getByToken(simVertexToken_ ,simVertices);
   Handle<EcalRecHitCollection> ecalRecHits;
   iEvent.getByToken(ecalRecHitsToken_, ecalRecHits);
   
   std::vector<SimTrack> theSimTracks;     
   std::vector<SimVertex> theSimVertices; 
   for(unsigned int i=0; i<simTracks->size() ; i++)
     {
       Ptr<SimTrack> track_it= simTracks->ptrAt( i );
       theSimTracks.push_back(*track_it);
     }
   for(unsigned int i=0; i<simVertices->size() ; i++)
     {
       Ptr<SimVertex> vert_it= simVertices->ptrAt( i );
       theSimVertices.push_back(*vert_it);
     }
   PhotonMCTruthFinder*  thePhotonMCTruthFinder_;
   thePhotonMCTruthFinder_ = new PhotonMCTruthFinder();
   std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find(theSimTracks,  theSimVertices);
   //   edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
   //if(!iEvent.getByToken("addPileupInfo",PupInfo))
   // {
   //   std::cout << "Pileup Summary not found" << std::endl;
   //   return;
   // }
   Handle<View<reco::Vertex>> recoVertices;
   iEvent.getByToken(vertexToken_,recoVertices);
   //Handle<View<reco::Photon>> photons;
   Handle<View<pat::Photon>> photons;
   iEvent.getByToken(photonToken_ ,photons);
   //Handle<View<reco::GenParticle>> genCollection;
   Handle<View<pat::PackedGenParticle>> genCollection;
   iEvent.getByToken(genToken_ ,genCollection);
   std::cout.flush();   
      
     /*************start Pileup**************************************************/
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   for(int i=0;i< 16; i++)
     {
       pu_size[i] =  0;
       pu_bx[i]   = -1;
     }
   /*
   int nPU=0; //used to assign index for pileup events
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) //should be 16 pileup bxs per event
     {
       pu_size[nPU]=PVI->getPU_NumInteractions();
       pu_bx[nPU] = PVI->getBunchCrossing();   
       //std::cout << "nPU " << nPU << " pu size " << pu_size[nPU] << std::endl;
       nPU++;
     }
   */
   //std::cout << nPU << " PU bx per event" << std::endl;
   

   n_event = iEvent.id().event();
   n_run   = iEvent.id().run();
   n_lumi  = iEvent.id().luminosityBlock();
   if (n_event ==24802)
     return;
   
   /***************************start photons *********************************************/
   for(int i=0;i<MAX_PHOTON;i++)
     {
       p_e12[i]=-1;
       p_e15[i]=-1;
       p_eAll[i]=-1;
       p1_e12[i]=-1;
       p1_e15[i]=-1;
       p1_eAll[i]=-1;
       p_ecalenergy[i]=-100;
       p_photonenergy[i]=-100;
       p_dr[i]=.2;
       p_eta[i]=-100;
       p_phi[i]=-100;
       p_SCeta[i]=-100;
       p_SCphi[i]=-100;
       p_r9[i]=-100;
       p_full5x5r9[i]=-100;
       p_energy_1x5[i]=-100;
       p_energy_2x5[i]=-100;
       p_energy_3x3[i]=-100;
       p_energy_full3x3[i]=-100;
       p_energy_5x5[i]=-100;
       p_energy_full5x5[i]=-100;
       p_pt_3x3[i]=-100;
       p_pt_full3x3[i]=-100;
       p_pt_5x5[i]=-100;
       p_pt_full5x5[i]=-100;
       p_esenergy[i]=-100;
       p_rawenergy[i]=-100;
       p_rawpt[i]=-100;
       p_pt[i]=-100;
       p_nRecHits[i]=0;
       p_isGap[i]=0;
       p_hOverE[i]=-100;
       p_sigmaIetaIeta[i]=-100;
       p_ecalRecHitSumEtConeDR03[i]=-100;
       p_hcalTowerSumEtConeDR03[i]= -100;
       p_trkSumPtSolidConeDR03[i]= -100;
       p_chargedIsoConeDR03[i] = -100;
       p_neutralIsoConeDR03[i] = -100;
       p_photonIsoConeDR03[i] = -100;
       
       g_eta[i]=-100;
       g_phi[i]=-100;
       g_pt[i]=-100;
       g_status[i] = -1;
       g_ePt[i]= 0;
       g_pPt[i]= 0;
       g_pdgid[i]=0;
       g_isConv[i]=-1;
       g_hdaughter[i]=0;
       g_zdaughter[i]=0;
     }

   for(int ix=0; ix<85; ix++){
     for(int iy=0; iy<360; iy++){
       for(int iz=0; iz<2; iz++){
	 g_calHits[ix][iy][iz] = 0.;
	 p_eRecHits[ix][iy][iz] = 0.;
	 p_RecHitInCluster[ix][iy][iz] = 0;
       }
     }
   }
   p_size=0;
   g_size=0;
   p_mass=0;
   p_rawmass=0;
   p_genVrawmass=0;
   p_genErawmass=0;
   p_5x5mass=0;  
   p_full5x5mass=0;
   p_3x3mass=0;
   p_full3x3mass=0;
   g_mass=0;
   //std::cout << "start gen collection" << std::endl;
   //   for(reco::GenParticleCollection::const_iterator gen_it=genCollection->begin(); gen_it!=genCollection->end() ; ++gen_it)
   reco::Vertex::Point hardVertex( 0, 0, 0 );
   for(int j=0;j<2;j++)
     {
       for(unsigned int i=0; i<genCollection->size() ; i++)
	 {
	   //Ptr<reco::GenParticle> gen_it= genCollection->ptrAt( i );
	   Ptr<pat::PackedGenParticle> gen_it= genCollection->ptrAt( i );
	   //       std::cout << "in gen particle collection" << std::endl;
	   /*	   if(gen_it->pdgId()==25 ) //higgs vertex
	     {
	       std::cout << "reco PV z: " << recoVertices->ptrAt(0)->z() << std::endl;
	       std::cout << "gen H   z: "  << gen_it->vertex().z() << std::endl;
	       }*/
	   if(gen_it->pdgId()!=22) //check for photons
	     continue;
	   //std::cout << "pt : " << gen_it->pt() << std::endl;
	   if(gen_it->pt()<15)
	     continue;
	   if(fabs(gen_it->eta())>1.444)
	     continue;
	   if(g_size>999)
	     break;
	   if(gen_it->status()!=1)
	     continue;
	   if(gen_it->pt()<=g_pt[j]) //only fill if current pT is higher than threshold
	     continue;
	   if(gen_it->pt()==g_pt[0]) //don't refill with highest pT (second highest only)
	     continue;
	   g_hdaughter[j]=1;
	   //	   vtx_dz =   recoVertices->ptrAt(0)->z() - hardVertex.z();
	   //vtx_dz =   recoVertices->ptrAt(0)->z();
	   g_eta[j] = gen_it->eta();
	   g_phi[j] = gen_it->phi();
	   g_pt[j] =  gen_it->pt();
	   g_pdgid[j] = 22;
	   g_status[j] = gen_it->status();
	   //std::cout << "g_pt("<<j<<"):  " << g_pt[j] << std::endl;
	   for(std::vector<PhotonMCTruth>::const_iterator it=mcPhotons.begin(); it!=mcPhotons.end() ; ++it)                              
	     {                                                                                                                           
	       //std::cout << it - mcPhotons.begin() << std::endl;                                                                       
	       if(fabs(gen_it->eta()-(*it).fourMomentum().pseudoRapidity())<.1)                
		 {                                                                                                                       
		   //std::cout << " geneta: "  << gen_it->eta() << " truth eta: " << (*it).fourMomentum().pseudoRapidity() << std::endl;
		   if(fabs(gen_it->phi()-(*it).fourMomentum().phi())<.1)                                                               
		     {                                                                                                                   
		       //std::cout << " genphi: "  << gen_it->phi() << " truth phi: " << (*it).fourMomentum().phi() << std::endl; 
		       if(fabs(gen_it->pt()-(*it).fourMomentum().et())<.1)          
			 {                                                                                                               
			   //  std::cout << " genpt: "  << gen_it->pt() << " truth pt: " << (*it).fourMomentum().et() << std::endl;    
			   g_isConv[j] = 0;                                                                                         
			   if(fabs(g_eta[j])<1.566)  //check for conversion in front of barrel                                      
			     if((*it).vertex().perp()<80)                                                                                
			       if((*it).isAConversion())                                                                                 
				 g_isConv[j] = 1;        
			   // std::cout << "g_isConv : " << g_isConv[g_size] << std::endl;                                               
			 }                                                                                                               
		     }                                                                                                                   
		 }                                                                                                                       
	     }                                     
	 }
     }
   if(g_pt[0]<20)
     g_size=0;
   else if(g_pt[1]<15)
     g_size=1;
   else
     {
       g_size=2;
       g_mass = sqrt(2*g_pt[0]*g_pt[1]*(cosh(g_eta[0]-g_eta[1]) - cos(g_phi[0]-g_phi[1])));
     }

   PCaloHitContainer::const_iterator ihit;
   Handle<PCaloHitContainer> pcaloeb;
   iEvent.getByToken(simCaloHitToken_,pcaloeb);
   //if(pcaloeb.isValid())
     {
       for (ihit=pcaloeb->begin(); ihit!=pcaloeb->end(); ihit++){
	 EBDetId  id     = EBDetId(ihit->id());
	 int    ix     = id.ietaAbs();
	 int    iy     = id.iphi();
	 int    iz     = id.zside();
	 if(iz==0) continue;
	 if(iz<0) iz=0;
	 g_calHits[ix-1][iy-1][iz] +=  ihit->energy();
	 //std::cout << "ix: " << ix << " iy " << iy << " iz " << iz << "   --ihit energy: " << ihit->energy();
       }
     }
   
   for(int j=0;j<2;j++)
     {
       for(unsigned int i=0; i<photons->size() ; i++)
	 {
	   Ptr<pat::Photon> it = photons->ptrAt( i );
	   //Ptr<reco::Photon> it = photons->ptrAt( i );
	   
	   /*
	   if(it->pt()<p_pt[j]) //only fill with highest pt
	     continue;
	   if(it->pt()==p_pt[0]) //don't refill with highest pT (second highest only)
	   continue;
	   */
	   float dphi=std::abs(it->phi()-g_phi[j]); if (dphi>float(M_PI)) dphi-=float(2*M_PI);  
	   dphi = sqrt((it->eta()-g_eta[j])*(it->eta()-g_eta[j]) + dphi*dphi);
	   if(dphi<p_dr[j])
	     p_dr[j] =dphi; 
	   else 
	     continue;
	   //apply cut based working points.
	   if(fabs(it->eta())<1.444) //EB
	     {
	     }
	   else
	     continue;
	   p_chargedIsoConeDR03[j] = it->chargedHadronIso();
	   p_neutralIsoConeDR03[j] = it->neutralHadronIso();
	   p_photonIsoConeDR03[j] = it->photonIso();
	   
	   p_eta[j] = it->eta();
	   p_phi[j] = it->phi();
	   p_SCeta[j] = it->superCluster()->eta();
	   p_SCphi[j] = it->superCluster()->phi();
	   p_ecalenergy[j] = it->superCluster()->energy();
	   p_photonenergy[j] = it->getCorrectedEnergy(it->getCandidateP4type());
	   p_r9[j] = it->r9();
	   p_full5x5r9[j] = it->full5x5_r9();
	   p_energy_1x5[j] = it->e1x5();
	   p_energy_2x5[j] = it->e2x5();
	   p_energy_3x3[j] = it->e3x3();
	   p_energy_full3x3[j] = it->full5x5_e3x3();
	   p_energy_5x5[j] = it->e5x5();
	   p_energy_full5x5[j] = it->full5x5_e5x5();
	   p_pt_3x3[j] = it->e3x3()/cosh(it->superCluster()->eta());
	   p_pt_full3x3[j] = it->full5x5_e3x3()/cosh(it->superCluster()->eta());
	   p_pt_5x5[j] = it->e5x5()/cosh(it->superCluster()->eta());
	   p_pt_full5x5[j] = it->full5x5_e5x5()/cosh(it->superCluster()->eta());
	   p_esenergy[j] = it->superCluster()->preshowerEnergy();
	   p_rawenergy[j]= it->superCluster()->rawEnergy();
	   p_rawpt[j]= it->superCluster()->rawEnergy()/cosh(it->superCluster()->eta());
	   p_pt[j]= it->pt();
	   try
	     {
	       it->superCluster()->clustersBegin();
	     }
	   catch(exception& e)
	     {
	       continue;
	     }
	   std::vector<double> cells;
	   std::vector<double> SIMcells;
	   if(p_pt[j]<20)
	     continue;
	   //for(reco::CaloClusterPtrVector::const_iterator bcItr = it->superCluster()->clustersBegin();
	     //bcItr != it->superCluster()->clustersEnd(); bcItr++)
	     //{
	     //	 p_nRecHits[j]+=(*bcItr)->size();
	     //}
	   //p_nRecHits[j]=it->superCluster()->hitsAndFractions().size();
	   for( auto&& detid_frac : it->superCluster()->hitsAndFractions() ) 
	     {
	       auto hit = ecalRecHits->find(detid_frac.first);
	       if ( hit == ecalRecHits->end() and detid_frac.first.subdetId() == DetId::Ecal ) 
		 {
		   std::cout << "invalid rechit" << std::endl;
		   continue;
		 }
	       //grab sim hit sums for rechit
	       EBDetId  id     = EBDetId(hit->id());
	       int    ix     = id.ietaAbs();
	       int    iy     = id.iphi();
	       int    iz     = id.zside();
	       //std::cout << "ix: " << ix << std::endl;
	       if(iz==0) continue;
	       if(iz<0) iz=0;
	       auto genHit = g_calHits[ix-1][iy-1][iz];
	       double frac = detid_frac.second;
	       cells.push_back(hit->energy()*frac);
	       SIMcells.push_back(genHit*frac);
	     }
	   std::sort(cells.begin(), cells.end());
	   std::sort(SIMcells.begin(), SIMcells.end());
	   for(auto it=cells.rbegin(); it!=cells.rend(); ++it) 
	     {	   
	       if ( p_nRecHits[j] < 12 )
		 {
		   p_e12[j] += *it;
		 }
	       if ( p_nRecHits[j] < 15 )
		 {
		   p_e15[j] += *it;
		 }
	       p_eAll[j] += *it;
	       p_nRecHits[j]++;
	     }
	   int nCellTemp=0;
	   for(auto it=SIMcells.rbegin(); it!=SIMcells.rend(); ++it) 
	     {	   
	       if ( nCellTemp < 12 )
		 {
		   p1_e12[j] += *it;
		 }
	       if ( nCellTemp < 15 )
		 {
		   p1_e15[j] += *it;
		 }
	       p1_eAll[j] += *it;
	       nCellTemp++;
	     }
	   
	   p_isGap[j] = it->isEBEEGap();
	   p_hOverE[j]= it->hadronicOverEm();
	   p_sigmaIetaIeta[j]= it->sigmaIetaIeta();
	   p_ecalRecHitSumEtConeDR03[j]= it->ecalRecHitSumEtConeDR03();
	   p_hcalTowerSumEtConeDR03[j]= it->hcalTowerSumEtConeDR03();
	   p_trkSumPtSolidConeDR03[j]= it->trkSumPtSolidConeDR03();
	 }
       //std::cout << "supercluster size " << it->superCluster()->sfdize() << std::endl;
     }
   if(p_pt[0]<20)
     p_size=0;
   else if(p_pt[1]<15)
     p_size=1;
   else
     {
       p_size=2;
       p_mass = sqrt(2*p_pt[0]*p_pt[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_rawmass = sqrt(2*p_rawpt[0]*p_rawpt[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_genVrawmass = sqrt(2*p_rawpt[0]*p_rawpt[1]*(cosh(g_eta[0]-g_eta[1]) - cos(g_phi[0]-g_phi[1])));
       p_genErawmass = sqrt(2*g_pt[0]*g_pt[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_5x5mass = sqrt(2*p_pt_5x5[0]*p_pt_5x5[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_full5x5mass = sqrt(2*p_pt_full5x5[0]*p_pt_full5x5[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_full3x3mass = sqrt(2*p_pt_full3x3[0]*p_pt_full3x3[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
       p_3x3mass = sqrt(2*p_pt_3x3[0]*p_pt_3x3[1]*(cosh(p_eta[0]-p_eta[1]) - cos(p_phi[0]-p_phi[1])));
     }	
   vtx_size = recoVertices->size();
   //std::cout << " vertex size: " << recoVertices->size() << std::endl;
   /************ end match photons ***********************/
   if(g_size==2&&g_mass>124.9&&g_mass<125.1)
     t_photon->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
Diphoton::beginJob(const edm::EventSetup& es, const edm::ParameterSet& ps)
{
  std::cout << "beginning job" << std::endl;
  
    
  
  
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Diphoton::endJob() 
{
}

void
Diphoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Diphoton);
