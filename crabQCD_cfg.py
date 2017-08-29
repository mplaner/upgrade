from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'upgrade_Mplaner1000_200_qcdRelVal'
config.General.workArea = 'upgrade1000_200_qcdRelVal'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'diphoton_qcd_cfg.py'

config.Data.allowNonValidInputDataset = True

config.Data.useParent = True
#config.Data.inputDataset =          '/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/PhaseIITDRSpring17MiniAOD-noPUNewCT_NewCT_91X_upgrade2023_realistic_v3-v3/MINIAODSIM'
#config.Data.inputDataset =          '/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/PhaseIITDRSpring17MiniAOD-PU200NewCT_NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM'
#config.Data.inputDataset =          ''
#config.Data.inputDataset =          ''
config.Data.inputDataset =          '/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_0_pre4-PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/MINIAODSIM'


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outputDatasetTag = 'CRAB3_test'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/mplaner/qcdRelVal_1000_200' 
config.Data.publication = True

config.Site.storageSite = 'T3_US_NotreDame'
