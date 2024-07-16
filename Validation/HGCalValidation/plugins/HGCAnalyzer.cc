// -*- C++ -*-
//
// Package:    analysis/HGCAnalyzer
// Class:      HGCAnalyzer
//
/**\class HGCAnalyzer HGCAnalyzer.cc analysis/HGCAnalyzer/plugins/HGCAnalyzer.cc

 Description: create ntuple for HGCAL clusters/tracksters study
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tielige Mengke
//         Created:  Thu, 22 Oct 2020 16:50:33 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <iomanip>
#include <string>
#include <array>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/RecoAlgos/interface/MultiVectorManager.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalGeometryMode.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/Point3D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/GenVector/PxPyPzM4D.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/LorentzVector.h"

using namespace std;

//
// class declaration
//

#define numTracksterTypes 3

class HGCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCAnalyzer(const edm::ParameterSet&);
  ~HGCAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void analyzeHGCalSimHit(edm::Handle<std::vector<PCaloHit>> const& simHits,
                          string detectorType,
                          DetId rechitID,
                          double rechitEn);

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  mutable hgcal::RecHitTools rhtools_;

  std::vector<const HGCalDDDConstants*> hgcCons_;
  std::vector<const HGCalGeometry*> hgcGeometry_;
  std::vector<std::string> geometrySource_;

  edm::InputTag eeSimHitSource, fhSimHitSource, bhSimHitSource;
  edm::EDGetTokenT<std::vector<PCaloHit>> eeSimHitToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> fhSimHitToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> bhSimHitToken_;

  std::vector<edm::EDGetTokenT<HGCRecHitCollection>> hitsTokens_;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_;
  // edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterstrkem_;
  // edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersem_;
  //edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterstrk_;
  //edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstershad_;
  //  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersmip_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersmrg_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersclue3dhigh_;
  // edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersclue3dlow_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> simtracksters_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> timelayercluster_;
  edm::EDGetTokenT<ticl::LayerClusterToCaloParticleAssociator> LCAssocByEnergyScoreProducer_;
  edm::EDGetTokenT<std::vector<reco::PFJet>> ak4Jets_;
  edm::EDGetTokenT<std::vector<reco::GenJet>> ak4GenJets_;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs_;
  TTree* tree_;

  unsigned int runId_;
  unsigned int eventId_;

  // std::string s_tracksters[9] = {"tracksterEM", "tracksterHAD", "tracksterMIP","tracksterMerge", "tracksterTrk", "tracksterTrkEM",
  std::string s_tracksters[numTracksterTypes] = {"tracksterMerge", "tracksterCLUE3DHigh", "simTrackster"};
  edm::EDGetTokenT<std::unordered_map<DetId, const unsigned int>> hitMap_;

  std::vector<int> gpdgId_;
  std::vector<ROOT::Math::PxPyPzEVector> gp_;
  std::vector<ROOT::Math::XYZPointF> gpPosition_;
  std::vector<ROOT::Math::PxPyPzEVector> pfj_;
  std::vector<ROOT::Math::PxPyPzEVector> genj_;

  std::vector<ROOT::Math::XYZPointF> lcPosition_;
  std::vector<float> lcEnergy_;
  std::vector<uint32_t> lcLayer_;
  std::vector<float> lcTime_;
  std::vector<float> lcTimeError_;
  std::vector<int> lcId_;
  std::vector<std::vector<int>> lcHits_;

  std::vector<ROOT::Math::XYZPointF> rhPosition_;
  std::vector<float> rhEnergy_;
  std::vector<float> rhTime_;
  std::vector<float> rhTimeError_;
  std::vector<uint32_t> rhLayer_;

  std::vector<std::vector<unsigned int>> vertices_[numTracksterTypes];
  std::vector<std::vector<unsigned int>> vertex_multiplicity_[numTracksterTypes];
  std::vector<int> seedIndex_[numTracksterTypes];
  std::vector<float> time_[numTracksterTypes];
  std::vector<float> timeError_[numTracksterTypes];
  std::vector<float> regressed_energy_[numTracksterTypes];
  std::vector<float> raw_energy_[numTracksterTypes];
  std::vector<float> raw_em_energy_[numTracksterTypes];
  std::vector<float> raw_pt_[numTracksterTypes];
  std::vector<float> raw_em_pt_[numTracksterTypes];
  std::vector<ROOT::Math::XYZVector> barycenter_[numTracksterTypes];
  std::vector<std::vector<float>> sigmas_[numTracksterTypes];
  std::vector<std::vector<float>> sigmasPCA_[numTracksterTypes];
  std::vector<std::vector<float>> id_probabilities_[numTracksterTypes];
  std::vector<float> sig_tmp, sigPCA_tmp, iP_tmp;

  std::vector<int> temp_;
  std::vector<int> temp2_;
  std::vector<float> temp3_;
  std::vector<float> temp4_;

  std::string t_name[14] = {"vertices",
                            "vertex_multiplicity",
                            "seedIndex",
                            "time",
                            "timeError",
                            "regressed_energy",
                            "raw_energy",
                            "raw_em_energy",
                            "raw_pt",
                            "raw_em_pt",
                            "barycenter",
                            "sigmas",
                            "sigmasPCA",
                            "id_probabilities"};

  std::vector<int> cpdgId_;
  std::vector<ROOT::Math::PxPyPzEVector> cp_;
  std::vector<std::vector<int>> cpSC_;
  std::vector<int> cpG4T0evt_;
  std::vector<int> cpG4T0bx_;
  std::vector<ROOT::Math::XYZPointF> cpOrigin_;

  std::vector<float> scEnergy_;
  std::vector<float> scSimEnergy_;
  std::vector<std::vector<int>> scHits_;
  std::vector<std::vector<float>> scHitsEnergyFrac_;

  std::vector<DetId> rhcontainer_;
  std::vector<DetId>::iterator it;

  std::vector<std::vector<float>> lc2cpScore_;
  std::vector<std::vector<int>> lc2cpId_;
  std::vector<std::vector<float>> lc2cpEnergy_;  //for testing idx

  std::vector<std::vector<float>> cp2lcScore_;
  std::vector<std::vector<int>> cp2lcId_;
  std::vector<std::vector<float>> cp2lcEnergy_;  //for testing idx

  std::vector<float> rhall_energy_;
  std::vector<int> rhall_layer_;
  // PF JETS
  std::vector<float> pfJetsPt_;
  std::vector<float> pfJetsEta_;
  std::vector<float> pfJetsPhi_;
  //

  // Gen JETS
  std::vector<float> ak4genJetsPt_;
  std::vector<float> ak4genJetsEta_;
  std::vector<float> ak4genJetsPhi_;

  int setZside_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HGCAnalyzer::HGCAnalyzer(const edm::ParameterSet& iConfig)
    : geomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      setZside_(iConfig.getUntrackedParameter<int>("setZside", 1)) {
  // (skold) geometrySource_ = iConfig.getParameter<std::vector<std::string>>("geometrySource");

  eeSimHitToken_ = consumes<std::vector<PCaloHit>>(iConfig.getParameter<edm::InputTag>("eeSimHitSource"));
  fhSimHitToken_ = consumes<std::vector<PCaloHit>>(iConfig.getParameter<edm::InputTag>("fhSimHitSource"));
  bhSimHitToken_ = consumes<std::vector<PCaloHit>>(iConfig.getParameter<edm::InputTag>("bhSimHitSource"));

  auto hitsTags = iConfig.getParameter<std::vector<edm::InputTag>>("hits");
  for (const auto& tag : hitsTags) {
    hitsTokens_.push_back(consumes<HGCRecHitCollection>(tag));
  }

  clusters_ = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("layer_clusters"));
  timelayercluster_ =
      consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("time_layerclusters"));
  //tracksterstrkem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrkem"));
  //trackstersem_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersem"));
  //trackstershad_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstershad"));
  //trackstersmip_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersmip"));
  //tracksterstrk_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("tracksterstrk"));
  trackstersmrg_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersmrg"));
  trackstersclue3dhigh_ =
      consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersCLUE3DHigh"));
  //trackstersclue3dlow_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersCLUE3DLow"));
  simtracksters_ = consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("simtracksters"));
  caloParticles_ = consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"));
  ak4Jets_ = consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("ak4Jets"));
  ak4GenJets_ = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("ak4GenJets"));

  hitMap_ = consumes<std::unordered_map<DetId, const unsigned int>>(iConfig.getParameter<edm::InputTag>("hitMapTag"));

  genParticles_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_particles"));

  LCAssocByEnergyScoreProducer_ = consumes<ticl::LayerClusterToCaloParticleAssociator>(
      iConfig.getParameter<edm::InputTag>("lcAssocByEnergyScoreProducer"));
}

HGCAnalyzer::~HGCAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void HGCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  runId_ = iEvent.id().run();
  eventId_ = iEvent.id().event();

  MultiVectorManager<HGCRecHit> rechitManager;
  for (const auto& token : hitsTokens_) {
    edm::Handle<HGCRecHitCollection> hitsHandle;
    iEvent.getByToken(token, hitsHandle);
    rechitManager.addVector(*hitsHandle);
  }

  edm::Handle<std::vector<reco::CaloCluster>> clusterHandle;
  iEvent.getByToken(clusters_, clusterHandle);
  auto const& clusters = *clusterHandle;

  edm::Handle<std::vector<CaloParticle>> caloParticleHandle;
  iEvent.getByToken(caloParticles_, caloParticleHandle);
  auto const& caloParticles = *caloParticleHandle;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> timelayerclusterHandle;
  iEvent.getByToken(timelayercluster_, timelayerclusterHandle);
  auto const& clustersTime = *timelayerclusterHandle;

  //edm::Handle<std::vector<ticl::Trackster>> trackstertrkemHandle;
  //iEvent.getByToken(tracksterstrkem_, trackstertrkemHandle);

  //edm::Handle<std::vector<ticl::Trackster>> tracksteremHandle;
  //iEvent.getByToken(trackstersem_, tracksteremHandle);

  //edm::Handle<std::vector<ticl::Trackster>> tracksterhadHandle;
  //iEvent.getByToken(trackstershad_, tracksterhadHandle);

  //edm::Handle<std::vector<ticl::Trackster>> trackstermipHandle;
  //iEvent.getByToken(trackstersmip_, trackstermipHandle);

  //edm::Handle<std::vector<ticl::Trackster>> trackstertrkHandle;
  //iEvent.getByToken(tracksterstrk_, trackstertrkHandle);

  edm::Handle<std::vector<ticl::Trackster>> trackstermrgHandle;
  iEvent.getByToken(trackstersmrg_, trackstermrgHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterclue3dhighHandle;
  iEvent.getByToken(trackstersclue3dhigh_, tracksterclue3dhighHandle);

  // edm::Handle<std::vector<ticl::Trackster>> tracksterclue3dlowHandle;
  // iEvent.getByToken(trackstersclue3dlow_, tracksterclue3dlowHandle);

  edm::Handle<std::vector<ticl::Trackster>> simtracksterHandle;
  iEvent.getByToken(simtracksters_, simtracksterHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle[numTracksterTypes] = {
      trackstermrgHandle, tracksterclue3dhighHandle, simtracksterHandle};

  edm::Handle<std::vector<reco::GenParticle>> genParticleHandle;
  edm::Handle<std::vector<reco::PFJet>> pfJetHandle;
  edm::Handle<std::vector<reco::GenJet>> genJetHandle;

  edm::Handle<std::unordered_map<DetId, const unsigned int>> hitMapHandle;
  iEvent.getByToken(hitMap_, hitMapHandle);
  const auto hitmap = *hitMapHandle;
  edm::Handle<ticl::LayerClusterToCaloParticleAssociator> LCAssocByEnergyScoreHandle;
  iEvent.getByToken(LCAssocByEnergyScoreProducer_, LCAssocByEnergyScoreHandle);

  ticl::RecoToSimCollection cpsInLayerClusterMap =
      LCAssocByEnergyScoreHandle->associateRecoToSim(clusterHandle, caloParticleHandle);
  ticl::SimToRecoCollection cPOnLayerMap =
      LCAssocByEnergyScoreHandle->associateSimToReco(clusterHandle, caloParticleHandle);

  // const CaloGeometry *caloGeom = &es.getData(geomToken_);

  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(geomToken_);
  rhtools_.setGeometry(*geom);

  int layers_ = rhtools_.lastLayerBH();

  // SimHits
  edm::Handle<std::vector<PCaloHit>> eeSimHits;
  iEvent.getByToken(eeSimHitToken_, eeSimHits);

  edm::Handle<std::vector<PCaloHit>> fhSimHits;
  iEvent.getByToken(fhSimHitToken_, fhSimHits);

  edm::Handle<std::vector<PCaloHit>> bhSimHits;
  iEvent.getByToken(bhSimHitToken_, bhSimHits);

  //all HGCAL rechits
  rhall_layer_.clear();
  rhall_energy_.clear();

  std::map<int, float> elsum;

  for (const auto& hit : rechitManager) {
    DetId detId = hit.id();
    int layer = rhtools_.getLayerWithOffset(detId) + layers_ * ((rhtools_.zside(detId) + 1) >> 1) - 1;
    elsum[layer] += hit.energy();
    if (detId.det() == DetId::Detector::HGCalEE) {
      if (hit.energy() > 30.0 || (hit.energy() > 5.00 && hit.energy() < 5.2)) {
        cout << "rechitsEE-check en " << hit.energy() << "  layer " << layer << "   lastlayer " << layers_;
        cout << "  layer " << rhtools_.getLayerWithOffset(detId) << endl;
        int subdet, zside, layer, wafer, celltype, cell;
        HGCalTestNumbering::unpackHexagonIndex(detId, subdet, zside, layer, wafer, celltype, cell);
        cout << "subdet-rechit " << subdet << "  zside " << zside << "  layer " << layer << "  wafer " << wafer
             << "  celltype " << celltype << "  cell " << cell;

        if (eeSimHits.isValid()) {
          analyzeHGCalSimHit(eeSimHits, "EE", detId, hit.energy());
        }
      }
    } else if (detId.det() == DetId::Detector::HGCalHSi) {
      if (hit.energy() > 40.0 || (hit.energy() > 5.0 && hit.energy() < 5.1)) {
        cout << "rechitsFH-check en " << hit.energy() << "  layer " << layer << "   lastlayer " << layers_;
        cout << "  layer " << rhtools_.getLayerWithOffset(detId) << endl;
        if (fhSimHits.isValid()) {
          analyzeHGCalSimHit(fhSimHits, "FH", detId, hit.energy());
        }
      }
    } else if (detId.det() == DetId::Detector::HGCalHSc) {
      if (hit.energy() > 1000.0) {
        cout << "rechitsFH-check en " << hit.energy() << "  layer " << layer << "   lastlayer " << layers_;
        cout << "  layer " << rhtools_.getLayerWithOffset(detId) << endl;
        if (bhSimHits.isValid()) {
          analyzeHGCalSimHit(bhSimHits, "BH", detId, hit.energy());
        }
      }
    }
  }

  std::map<int, float>::iterator iter;
  for (iter = elsum.begin(); iter != elsum.end(); ++iter) {
    rhall_energy_.push_back(iter->second);
    rhall_layer_.push_back(iter->first);
  }

  elsum.clear();

  iEvent.getByToken(genParticles_, genParticleHandle);

  //genParticles
  gpdgId_.clear();
  gp_.clear();
  gpPosition_.clear();

  auto const& genParticles = *genParticleHandle;
  for (auto const& gp : genParticles) {
    ROOT::Math::PxPyPzEVector p4;
    p4.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy());
    gp_.push_back(p4);
    gpdgId_.push_back(gp.pdgId());
    ROOT::Math::XYZPointF v;
    v.SetXYZ(gp.vx(), gp.vy(), gp.vz());
    gpPosition_.push_back(v);
  }

  //layercluster
  int rh_idx = 0;
  lcEnergy_.clear();
  lcLayer_.clear();
  lcTime_.clear();
  lcTimeError_.clear();
  lcId_.clear();
  lcPosition_.clear();
  lcHits_.clear();
  rhEnergy_.clear();
  rhTime_.clear();
  rhTimeError_.clear();
  rhPosition_.clear();
  rhcontainer_.clear();
  rhLayer_.clear();
  lc2cpScore_.clear();
  lc2cpId_.clear();
  lc2cpEnergy_.clear();
  for (auto const& lc : clusters) {
    int lc_idx = &lc - &clusters[0];
    if (setZside_ * lc.z() < 0)
      continue;
    const auto firstHitDetId = lc.hitsAndFractions()[0].first;
    int layerId = rhtools_.getLayerWithOffset(firstHitDetId) + layers_ * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
    const auto& hits_and_fractions = lc.hitsAndFractions();

    for (const auto& it_h : hits_and_fractions) {
      const auto rh_detid = it_h.first;
      if (hitmap.find(rh_detid) == hitmap.end()) {
        continue;
      }
      if ((rh_detid.det() != DetId::Forward) && (rh_detid.det() != DetId::HGCalEE) &&
          (rh_detid.det() != DetId::HGCalHSi) && (rh_detid.det() != DetId::HGCalHSc)) {
        std::cout << "Not HGCAL detector" << std::endl;
        continue;
      }
      auto global = rhtools_.getPosition(it_h.first);
      ROOT::Math::XYZPointF rh_p;
      rh_p.SetXYZ(global.x(), global.y(), global.z());
      rhPosition_.push_back(rh_p);
      rhTime_.push_back(rechitManager[hitmap.at(rh_detid)].time());
      rhTimeError_.push_back(rechitManager[hitmap.at(rh_detid)].timeError());
      rhEnergy_.push_back(rechitManager[hitmap.at(rh_detid)].energy());
      int rhl = rhtools_.getLayerWithOffset(rh_detid) + layers_ * ((rhtools_.zside(rh_detid) + 1) >> 1) - 1;
      rhLayer_.push_back(rhl);
      temp_.push_back(rh_idx);
      rhcontainer_.push_back(rh_detid);
      rh_idx++;
    }
    lcHits_.push_back(temp_);
    temp_.clear();  //reuse container
    lcEnergy_.push_back(lc.energy());
    lcLayer_.push_back(layerId);
    lcId_.push_back(lc_idx);
    ROOT::Math::XYZPointF lc_p;
    lc_p.SetXYZ(lc.x(), lc.y(), lc.z());
    lcPosition_.push_back(lc_p);
    lcTime_.push_back(clustersTime.get(lc_idx).first);
    lcTimeError_.push_back(clustersTime.get(lc_idx).second);

    //metrics
    //std::cout<<"lc_idx "<<lc_idx<<std::endl;
    const edm::Ref<std::vector<reco::CaloCluster>> lcRef(clusterHandle, lc_idx);
    const auto& cpsIt = cpsInLayerClusterMap.find(lcRef);
    if (cpsIt == cpsInLayerClusterMap.end()) {
      lc2cpScore_.push_back(temp3_);
      lc2cpId_.push_back(temp2_);
      lc2cpEnergy_.push_back(temp4_);
    } else {
      const auto& cps = cpsIt->val;
      for (const auto& cpPair : cps) {
        if (cpPair.first->g4Tracks()[0].eventId().event() != 0 or
            cpPair.first->g4Tracks()[0].eventId().bunchCrossing() != 0)
          continue;

        auto const& cp_linked =
            std::find_if(std::begin(cPOnLayerMap[cpPair.first]),
                         std::end(cPOnLayerMap[cpPair.first]),
                         [&lcRef](const std::pair<edm::Ref<reco::CaloClusterCollection>, std::pair<float, float>>& p) {
                           return p.first == lcRef;
                         });

        if (cp_linked == cPOnLayerMap[cpPair.first].end())
          continue;

        temp3_.push_back(cpPair.second);
        temp2_.push_back(cpPair.first.index());
        temp4_.push_back(cp_linked->second.first);
      }
      lc2cpScore_.push_back(temp3_);
      lc2cpId_.push_back(temp2_);
      lc2cpEnergy_.push_back(temp4_);
      temp3_.clear();
      temp2_.clear();
      temp4_.clear();
    }
  }
  //
  //caloparticles
  int sc_idx = 0;
  cp_.clear();
  cpdgId_.clear();
  cpSC_.clear();
  cpOrigin_.clear();
  cpG4T0evt_.clear();
  cpG4T0bx_.clear();
  scEnergy_.clear();
  scHits_.clear();
  scHitsEnergyFrac_.clear();
  cp2lcScore_.clear();
  cp2lcId_.clear();
  cp2lcEnergy_.clear();
  for (const auto& cp : caloParticles) {
    if (cp.g4Tracks()[0].eventId().event() != 0 or cp.g4Tracks()[0].eventId().bunchCrossing() != 0)
      continue;  //saving caloparticles from hard scattering only
    if (setZside_ * cp.g4Tracks()[0].trackerSurfacePosition().Z() < 0)
      continue;
    int cp_idx = &cp - &caloParticles[0];
    cpG4T0evt_.push_back(cp.g4Tracks()[0].eventId().event());
    cpG4T0bx_.push_back(cp.g4Tracks()[0].eventId().bunchCrossing());
    ROOT::Math::PxPyPzEVector p4;
    p4.SetPxPyPzE(cp.px(), cp.py(), cp.pz(), cp.energy());
    cp_.push_back(p4);
    cpdgId_.push_back(cp.pdgId());
    //initial g4track vertex
    ROOT::Math::XYZPointF cp_o;
    cp_o.SetXYZ(cp.g4Tracks()[0].trackerSurfacePosition().X(),
                cp.g4Tracks()[0].trackerSurfacePosition().Y(),
                cp.g4Tracks()[0].trackerSurfacePosition().Z());
    cpOrigin_.push_back(cp_o);
    const SimClusterRefVector& simClusterRefVector = cp.simClusters();
    for (const auto& sc : simClusterRefVector) {
      const SimCluster& simCluster = (*(sc));
      const auto& hits_and_fractions = simCluster.hits_and_fractions();
      for (const auto& it_h : hits_and_fractions) {
        DetId rh_detid = (it_h.first);
        if (!hitmap.count(rh_detid))
          continue;
        if ((rh_detid.det() != DetId::Forward) && (rh_detid.det() != DetId::HGCalEE) &&
            (rh_detid.det() != DetId::HGCalHSi) && (rh_detid.det() != DetId::HGCalHSc)) {
          std::cout << "Not HGCAL detector" << std::endl;
          continue;
        }

        it = std::find(rhcontainer_.begin(), rhcontainer_.end(), rh_detid);
        if (it != rhcontainer_.end()) {
          temp_.push_back(it - rhcontainer_.begin());
          temp3_.push_back(it_h.second);
        } else {
          auto global = rhtools_.getPosition(it_h.first);
          ROOT::Math::XYZPointF rh_p;
          rh_p.SetXYZ(global.x(), global.y(), global.z());
          rhPosition_.push_back(rh_p);
          rhTime_.push_back(rechitManager[hitmap.at(rh_detid)].time());
          rhTimeError_.push_back(rechitManager[hitmap.at(rh_detid)].timeError());
          rhEnergy_.push_back(rechitManager[hitmap.at(rh_detid)].energy());
          int rhl = rhtools_.getLayerWithOffset(rh_detid) + layers_ * ((rhtools_.zside(rh_detid) + 1) >> 1) - 1;
          rhLayer_.push_back(rhl);
          temp_.push_back(rh_idx);
          temp3_.push_back(it_h.second);
          rh_idx++;
        }
      }
      scEnergy_.push_back(simCluster.energy());
      scSimEnergy_.push_back(simCluster.simEnergy());
      scHits_.push_back(temp_);
      scHitsEnergyFrac_.push_back(temp3_);
      temp_.clear();
      temp3_.clear();
      temp2_.push_back(sc_idx);
    }
    cpSC_.push_back(temp2_);
    temp2_.clear();

    const edm::Ref<CaloParticleCollection> cpRef(caloParticleHandle, cp_idx);
    const auto& lcsIt = cPOnLayerMap.find(cpRef);
    if (lcsIt == cPOnLayerMap.end()) {
      cp2lcScore_.push_back(temp4_);
      cp2lcId_.push_back(temp2_);
      cp2lcEnergy_.push_back(temp3_);
    } else {
      const auto& lcs = lcsIt->val;
      for (const auto& lcPair : lcs) {
        temp4_.push_back(lcPair.second.second);
        temp2_.push_back(lcPair.first.index());
        temp3_.push_back(lcPair.second.first);
      }
      cp2lcScore_.push_back(temp4_);
      cp2lcId_.push_back(temp2_);
      cp2lcEnergy_.push_back(temp3_);
      temp4_.clear();
      temp2_.clear();
      temp3_.clear();
    }
  }

  //tracksters
  // for (int i=0;i<7;i++){
  // for (int i=0;i<6;i++){
  for (int i = 0; i < numTracksterTypes; i++) {
    vertices_[i].clear();
    vertex_multiplicity_[i].clear();
    time_[i].clear();
    timeError_[i].clear();
    regressed_energy_[i].clear();
    raw_energy_[i].clear();
    raw_em_energy_[i].clear();
    raw_pt_[i].clear();
    raw_em_pt_[i].clear();
    auto const& tracksters = *tracksterHandle[i];
    for (auto const& tkr : tracksters) {
      if (!tkr.vertices().empty()) {
        if (tkr.barycenter().z() * setZside_ < 0)
          continue;
        vertices_[i].push_back(tkr.vertices());
        std::vector<unsigned int> temp_vec;
        for (auto& vm : tkr.vertex_multiplicity()) {
          temp_vec.push_back(unsigned(vm));
        }
        vertex_multiplicity_[i].push_back(temp_vec);
        seedIndex_[i].push_back(tkr.seedIndex());
        time_[i].push_back(tkr.time());
        timeError_[i].push_back(tkr.timeError());
        regressed_energy_[i].push_back(tkr.regressed_energy());
        raw_energy_[i].push_back(tkr.raw_energy());
        raw_em_energy_[i].push_back(tkr.raw_em_energy());
        raw_pt_[i].push_back(tkr.raw_pt());
        raw_em_pt_[i].push_back(tkr.raw_em_pt());
        //			barycenter_[i].push_back(tkr.barycenter());
        for (int j = 0; j < 3; j++) {
          sig_tmp.push_back(tkr.sigmas()[j]);
          sigPCA_tmp.push_back(tkr.sigmasPCA()[j]);
        }
        for (int j = 0; j < 8; j++)
          iP_tmp.push_back(tkr.id_probabilities()[j]);
        sigmas_[i].push_back(sig_tmp);
        sigmasPCA_[i].push_back(sigPCA_tmp);
        id_probabilities_[i].push_back(iP_tmp);
        sig_tmp.clear();
        sigPCA_tmp.clear();
        iP_tmp.clear();
      }
    }
  }
  pfJetsPt_.clear();
  pfJetsEta_.clear();
  pfJetsPhi_.clear();
  iEvent.getByToken(ak4Jets_, pfJetHandle);
  auto const& PFJets = *pfJetHandle;
  for (auto const& pfJets : PFJets) {
    cout << " pf jets Pt : " << pfJets.pt() << " eta : " << pfJets.eta() << " phi : " << pfJets.phi() << endl;
    /*" emFraction : "<<pfJets.emEnergyFraction() "<< hadEnergyInHB: "<<pfJets.hadEnergyInHB()<<endl;*/
    pfJetsPt_.push_back(pfJets.pt());
    pfJetsEta_.push_back(pfJets.eta());
    pfJetsPhi_.push_back(pfJets.phi());
    //ROOT::Math::PtEtaPhiMVector p4;
    //double pt = pfJets.pt();
    //double eta = pfJets.eta();
    //double phi = pfJets.phi();
    //double energy = pfJets.energy();
    ROOT::Math::PxPyPzEVector p4;
    p4.SetPxPyPzE(pfJets.px(), pfJets.py(), pfJets.pz(), pfJets.energy());
    pfj_.push_back(p4);
    //p4.SetPtEtaPhiM(pt,eta,phi,mass);
    //pfj_.push_back(p4);

    //TLorentzVector jet;
    //jet.SetPtEtaPhiE(pt, eta, phi, energy);
  }
  ak4genJetsPt_.clear();
  ak4genJetsEta_.clear();
  ak4genJetsPhi_.clear();
  iEvent.getByToken(ak4GenJets_, genJetHandle);
  auto const& GENJets = *genJetHandle;
  for (auto const& genJets : GENJets) {
    cout << " gen jets Pt : " << genJets.pt() << " gen eta : " << genJets.eta() << " gen phi : " << genJets.phi()
         << endl;
    /*" emFraction : "<<pfJets.emEnergyFraction() "<< hadEnergyInHB: "<<pfJets.hadEnergyInHB()<<endl;*/
    ak4genJetsPt_.push_back(genJets.pt());
    ak4genJetsEta_.push_back(genJets.eta());
    ak4genJetsPhi_.push_back(genJets.phi());
    //ROOT::Math::PtEtaPhiMVector p4;
    //double pt = genJets.pt();
    //double eta = genJets.eta();
    //double phi = genJets.phi();
    //double energy = genJets.energy();
    ROOT::Math::PxPyPzEVector p4;
    p4.SetPxPyPzE(genJets.px(), genJets.py(), genJets.pz(), genJets.energy());
    genj_.push_back(p4);
    //p4.SetPtEtaPhiM(pt,eta,phi,mass);
    //pfj_.push_back(p4);

    //TLorentzVector genjet;
    //genjet.SetPtEtaPhiE(pt, eta, phi, energy);
  }

  tree_->Fill();
  //std::cout<<"finished one event"<<std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void HGCAnalyzer::beginJob() {
  tree_ = fs_->make<TTree>("tree", "TICL objects");
  //tree_->SetAutoSave(0);

  tree_->Branch("runId", &runId_);
  tree_->Branch("eventId", &eventId_);

  tree_->Branch("genParticle", &gp_);
  tree_->Branch("genParticlePosition", &gpPosition_);
  tree_->Branch("genPdgId", &gpdgId_);

  tree_->Branch("lcEnergy", &lcEnergy_);
  tree_->Branch("lcLayer", &lcLayer_);
  tree_->Branch("lcId", &lcId_);
  tree_->Branch("lcPosition", &lcPosition_);
  tree_->Branch("lcTime", &lcTime_);
  tree_->Branch("lcTimeError", &lcTimeError_);
  tree_->Branch("lcHits", &lcHits_);
  tree_->Branch("lc2cpScore", &lc2cpScore_);
  tree_->Branch("lc2cpId", &lc2cpId_);
  tree_->Branch("lc2cpEnergy", &lc2cpEnergy_);

  //	tree_->Branch("rhTime", &rhTime_);
  //	tree_->Branch("rhTimeError", &rhTimeError_);
  //tree_->Branch("rhEnergy", &rhEnergy_);
  tree_->Branch("rhPosition", &rhPosition_);
  tree_->Branch("rhLayer", &rhLayer_);

  for (int i = 0; i < numTracksterTypes; i++) {
    std::string s_temp[14]{};
    for (int j = 0; j < 14; j++) {
      s_temp[j] = s_tracksters[i] + "_" + t_name[j];
    }
    //		std::cout<<"i= "<<i<<s_temp[0]<<std::endl;
    tree_->Branch(s_temp[0].c_str(), &vertices_[i]);
    tree_->Branch(s_temp[1].c_str(), &vertex_multiplicity_[i]);
    tree_->Branch(s_temp[2].c_str(), &seedIndex_[i]);
    tree_->Branch(s_temp[3].c_str(), &time_[i]);
    tree_->Branch(s_temp[4].c_str(), &timeError_[i]);
    tree_->Branch(s_temp[5].c_str(), &regressed_energy_[i]);
    tree_->Branch(s_temp[6].c_str(), &raw_energy_[i]);
    tree_->Branch(s_temp[7].c_str(), &raw_em_energy_[i]);
    tree_->Branch(s_temp[8].c_str(), &raw_pt_[i]);
    tree_->Branch(s_temp[9].c_str(), &raw_em_pt_[i]);
    tree_->Branch(s_temp[10].c_str(), &barycenter_[i]);
    tree_->Branch(s_temp[11].c_str(), &sigmas_[i]);
    tree_->Branch(s_temp[12].c_str(), &sigmasPCA_[i]);
    tree_->Branch(s_temp[13].c_str(), &id_probabilities_[i]);
  }
  tree_->Branch("cpdgId", &cpdgId_);
  tree_->Branch("cp", &cp_);
  tree_->Branch("cpSC", &cpSC_);
  tree_->Branch("cpOrigin", &cpOrigin_);
  tree_->Branch("cpG4T0evt", &cpG4T0evt_);
  tree_->Branch("cpG4T0bx", &cpG4T0bx_);
  tree_->Branch("cp2lcScore", &cp2lcScore_);
  tree_->Branch("cp2lcId", &cp2lcId_);
  tree_->Branch("cp2lcEnergy", &cp2lcEnergy_);

  tree_->Branch("scEnergy", &scEnergy_);
  tree_->Branch("scSimEnergy", &scSimEnergy_);
  tree_->Branch("scHits", &scHits_);
  tree_->Branch("scHitsEnergyFrac", &scHitsEnergyFrac_);

  tree_->Branch("rhall_energy", &rhall_energy_);
  tree_->Branch("rhall_layer", &rhall_layer_);

  tree_->Branch("pfJetPt", &pfJetsPt_);
  tree_->Branch("pfJetEta", &pfJetsEta_);
  tree_->Branch("pfJetPhi", &pfJetsPhi_);
  tree_->Branch("pfJet", &pfj_);

  tree_->Branch("genJetPt", &ak4genJetsPt_);
  tree_->Branch("genJetEta", &ak4genJetsEta_);
  tree_->Branch("genJetPhi", &ak4genJetsPhi_);
  tree_->Branch("genJet", &genj_);
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCAnalyzer::endJob() { std::cout << "finished" << std::endl; }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //desc.setUnknown(); removing lc and timelc for time being
  desc.setAllowAnything();
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
  //desc.add<edm::InputTag>("trackstersem", edm::InputTag("ticlTrackstersEM"));
  //desc.add<edm::InputTag>("trackstershad", edm::InputTag("ticlTrackstersHAD"));
  // desc.add<edm::InputTag>("trackstersmip", edm::InputTag("ticlTrackstersMIP"));
  desc.add<edm::InputTag>("trackstersmrg", edm::InputTag("ticlCandidate"));
  //desc.add<edm::InputTag>("tracksterstrk", edm::InputTag("ticlTrackstersTrk"));
  //desc.add<edm::InputTag>("tracksterstrkem", edm::InputTag("ticlTrackstersTrkEM"));
  desc.add<edm::InputTag>("trackstersCLUE3DHigh", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  // desc.add<edm::InputTag>("trackstersCLUE3DLow", edm::InputTag("ticlTrackstersCLUE3DLow"));
  desc.add<edm::InputTag>("simtracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("gen_particles", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("hitMapTag", edm::InputTag("recHitMapProducer", "hgcalRecHitMap"));
  desc.add<edm::InputTag>("lcAssocByEnergyScoreProducer", edm::InputTag("lcAssocByEnergyScoreProducer"));
  desc.add<std::vector<edm::InputTag>>("hits",
                                       {edm::InputTag("HGCalRecHit", "HGCEERecHits"),
                                        edm::InputTag("HGCalRecHit", "HGCHEFRecHits"),
                                        edm::InputTag("HGCalRecHit", "HGCHEBRecHits")});
  desc.add<edm::InputTag>("eeSimHitSource", edm::InputTag("g4SimHits", "HGCHitsEE"));
  desc.add<edm::InputTag>("fhSimHitSource", edm::InputTag("g4SimHits", "HGCHitsHEfront"));
  desc.add<edm::InputTag>("bhSimHitSource", edm::InputTag("g4SimHits", "HGCHitsHEback"));
  desc.add<edm::InputTag>("ak4Jets", edm::InputTag("ak4PFJets"));
  desc.add<edm::InputTag>("ak4GenJets", edm::InputTag("ak4GenJets"));

  std::vector<std::string> source = {"HGCalEESensitive", "HGCalHESiliconSensitive", "HGCalHEScintillatorSensitive"};
  desc.add<std::vector<std::string>>("geometrySource", source);

  descriptions.addDefault(desc);
}

void HGCAnalyzer::analyzeHGCalSimHit(edm::Handle<std::vector<PCaloHit>> const& simHits,
                                     string detectorType,
                                     DetId rechitId,
                                     double rechitEn) {
  // const HGCalTopology& hTopo = hgcGeometry_[idet]->topology();
  cout << "HGCAnalyzer::analyzeHGCalSimHit   detector type=" << detectorType << "   simHits.size()=" << simHits->size()
       << endl;
  double esum = 0.0;
  for (std::vector<PCaloHit>::const_iterator simHit = simHits->begin(); simHit != simHits->end(); ++simHit) {
    int subdet, zside, layer, wafer, celltype, cell;
    HGCalTestNumbering::unpackHexagonIndex(simHit->id(), subdet, zside, layer, wafer, celltype, cell);
    // std::pair<float, float> xy = hgcCons_[idet]->locateCell(cell, layer, wafer, false);
    // float zp = hgcCons_[idet]->waferZ(layer, false);
    //  if (zside < 0)  zp = -zp;
    //  float xp = (zp < 0) ? -xy.first / 10 : xy.first / 10;
    //  float yp = xy.second / 10.0;
    // if(simHit->energy()<0.005) continue;
    if (simHit->id() == rechitId) {
      cout << "subdet " << subdet << "  zside " << zside << "  layer " << layer << "  wafer " << wafer << "  celltype "
           << celltype << "  cell " << cell;
      cout << "  e " << simHit->energy() << "   time " << simHit->time();
      cout << "  " << endl;
      esum = esum + simHit->energy();
      // cout<<"   x "<<xp<<"   y "<<yp<<"  z "<<zp<<endl;
    }
  }
  if (esum > 0.0) {
    int subdet, zside, layer, wafer, celltype, cell;
    HGCalTestNumbering::unpackHexagonIndex(rechitId, subdet, zside, layer, wafer, celltype, cell);
    cout << "aaaDetType=" << detectorType;
    cout << "   rechitEn " << rechitEn << "   simhitSum  " << esum << "  rec/sim =" << rechitEn / esum;
    cout << "   layer " << layer << "  wafer " << wafer << "  celltype " << celltype << "  rawID " << rechitId.rawId();
    cout << " " << endl;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCAnalyzer);
