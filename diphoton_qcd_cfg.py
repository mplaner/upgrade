
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Demo")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
from EventFilter.EcalRawToDigi.EcalUnpackerData_cfi import ecalEBunpacker
process.ecalDigis = ecalEBunpacker.clone()
process.load("RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoEcal.Configuration.RecoEcal_cff') #superclustering
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('RecoEgamma.Configuration.RecoEgamma_cff') #recphoton+electron
process.load('RecoJets.Configuration.CaloTowersRec_cff')
process.load('MagneticField.ParametrizedEngine.parametrizedMagneticField_OAE_3_8T_cfi')
process.load('Configuration.StandardSequences.Services_cff')


from RecoEgamma.EgammaPhotonProducers.photonSequence_cff import *

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '') 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

#run over multiple files at once, put full path name in eosFileList.txt
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        "root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_3_0_pre4/RelValQCD_Pt-15To7000_Flat_14TeV/MINIAODSIM/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/B08D7B77-5C89-E711-9400-0242AC130002.root"
        #"/store/mc/PhaseIITDRSpring17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v3/110000/32F8BA11-4373-E711-AC6E-1CC1DE1CEDB2.root"
        ),
                            secondaryFileNames=cms.untracked.vstring("/store/relval/CMSSW_9_3_0_pre4/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/8A3BD9F4-2B89-E711-823A-0242AC130002.root",
                                                                     "/store/relval/CMSSW_9_3_0_pre4/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/8CD22FB9-2789-E711-B475-0242AC130002.root"
                                                                     )
                            )
#root://cms-xrd-global.cern.ch/

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService", fileName = cms.string("local_output_ntuple.root") ) #output file name
process.ntuple = cms.EDAnalyzer('QCDdiphoton')
process.ntuple.photonTag = cms.InputTag("slimmedPhotons")
process.ntuple.genTag = cms.InputTag("packedGenParticles")
process.ntuple.vertexTag = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ntuple.ecalRecHits   = cms.InputTag("reducedEgamma:reducedEBRecHits")
#process.ntuple.photonTag = cms.InputTag("photons")
#process.ntuple.photonTag = cms.InputTag("gedPhotons")
#process.ntuple.genTag = cms.InputTag("genParticles")
#process.ntuple.vertexTag = cms.InputTag("offlinePrimaryVertices")
process.p = cms.Path(process.ntuple) #run the ntuples

process.dump=cms.EDAnalyzer('EventContentAnalyzer') 
#process.p = cms.Path(process.dump) #check event content of input file
