from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'FinalB_upgrade_Mplaner4500_200'
config.General.workArea = 'FinalB_upgrade4500_200'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'diphoton_simhitTest.py'

config.Data.useParent = True
#config.Data.inputDataset =          '
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA300-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA1000-v1/MINIAODSIM'
                                     #/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000-v1/MINIAODSIM
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000Ult-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA4500Ult-v1/MINIAODSIM'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA300-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA1000-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA3000Ult-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA4500Ult-v1/MINIAODSIM'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA300-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA1000-v1/MINIAODSIM'
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000Ult-v1/MINIAODSIM'
config.Data.inputDataset =           '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA4500Ult-v1/MINIAODSIM'
# #don't run##config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000-v1/MINIAODSIM'

#pu200 age300 ext
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA300-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA300-v1/GEN-SIM-DIGI-RAW'
#pu200 age 1000 ext
#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA1000-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA1000-v1/GEN-SIM-DIGI-RAW'
#config.Data.inputDataset =          '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v3/MINIAODSIM'
#config.Data.secondaryInputDataset = '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v3/GEN-SIM-RECO'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA3000Ultimate-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200EA3000Ultimate-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000Ult-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000Ult-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA4500Ult-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA4500Ult-v1/GEN-SIM-DIGI-RAW'


#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPU_ext1-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPU_ext1-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA300-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA300-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA1000-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA1000-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA3000Ultimate-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA3000Ultimate-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA4500Ultimate-v1/MINIAODSIM'
#config.Data.secondaryInputDataset = '/RelValH125GGgluonfusion_14/CMSSW_9_1_1_patch3-91X_upgrade2023_realistic_v3_D17noPUExtEA4500Ultimate-v1/GEN-SIM-DIGI-RAW'

#config.Data.inputDataset =          '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRSpring17MiniAOD-PU140_91X_upgrade2023_realistic_v3-v3/MINIAODSIM'
#config.Data.secondaryInputDataset = '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRSpring17DR-PU140_91X_upgrade2023_realistic_v3-v3/GEN-SIM-RECO'


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outputDatasetTag = 'CRAB3_test'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/mplaner/B4500_200' 
config.Data.publication = True

config.Site.storageSite = 'T3_US_NotreDame'
