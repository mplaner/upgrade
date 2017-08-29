
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

#run over multiple files at once, put full path name in eosFileList.txt
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/MINIAODSIM/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/D4CA8680-4C6A-E711-8545-D067E5F91B8A.root"
        #"/store/mc/PhaseIITDRSpring17DR/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v3/110000/32F8BA11-4373-E711-AC6E-1CC1DE1CEDB2.root"
        ),
                            secondaryFileNames=cms.untracked.vstring(
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/0A8247C7-F369-E711-BBBD-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/1AAE2B86-F669-E711-AC67-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/3AF42DAB-F569-E711-BE20-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/409C62D2-F569-E711-A0F4-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/5687DA30-F669-E711-9456-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/580B3067-F869-E711-A162-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/5AA09852-F469-E711-9F1A-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/5C2A3FED-F469-E711-936F-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/685B72DA-F669-E711-B6F0-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/6A7C2B43-F869-E711-960B-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/6C94A32F-F769-E711-B3D2-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/6CA675C7-F369-E711-AC17-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/788B31E1-F569-E711-BE20-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/7CDDF5E9-F569-E711-AF33-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/8A5F3741-F569-E711-8063-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/92C23393-F569-E711-873F-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/983AF7CA-EF69-E711-91FB-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/A0763FF3-F469-E711-9816-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/A4B1191B-F069-E711-973F-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/AA6A2C6C-F069-E711-873F-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/ACF9C1CF-F769-E711-A773-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/B09FE1F6-F769-E711-A6F7-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/B0F8DA83-F569-E711-93BC-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/CE95B6BA-F569-E711-8D3B-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/D09FBC46-F869-E711-90CE-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/D0E3034D-F869-E711-8863-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/D26681C7-F369-E711-9B2B-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/DE926C4A-F869-E711-9535-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/E07A0C4D-F869-E711-9AD6-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/E8F4690C-F569-E711-9D78-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/ECF57899-EF69-E711-83D0-0242AC130002.root",
        "/store/relval/CMSSW_9_1_1_patch3/RelValH125GGgluonfusion_14/GEN-SIM-DIGI-RAW/PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA300-v1/00000/FADB3143-F869-E711-93E2-0242AC130002.root"
        )
                            )
#root://cms-xrd-global.cern.ch/

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService", fileName = cms.string("local_output_ntuple.root") ) #output file name
process.ntuple = cms.EDAnalyzer('Diphoton')
process.ntuple.photonTag = cms.InputTag("slimmedPhotons")
process.ntuple.genTag = cms.InputTag("packedGenParticles")
process.ntuple.vertexTag = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ntuple.simVertexTag = cms.InputTag("g4SimHits")
process.ntuple.simTrackTag = cms.InputTag("g4SimHits")
process.ntuple.simCaloHitTag = cms.InputTag("g4SimHits:EcalHitsEB")
process.ntuple.ecalRecHits   = cms.InputTag("reducedEgamma:reducedEBRecHits")
#process.ntuple.photonTag = cms.InputTag("photons")
#process.ntuple.photonTag = cms.InputTag("gedPhotons")
#process.ntuple.genTag = cms.InputTag("genParticles")
#process.ntuple.vertexTag = cms.InputTag("offlinePrimaryVertices")
process.p = cms.Path(process.ntuple) #run the ntuples

process.dump=cms.EDAnalyzer('EventContentAnalyzer') 
#process.p = cms.Path(process.dump) #check event content of input file
