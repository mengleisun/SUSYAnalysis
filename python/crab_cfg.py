from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PhotonJet'
config.General.workArea = 'PhotonJet'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'PhotonJet.py'

config.section_("Data")
config.Data.outputPrimaryDataset = 'PhotonJet'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 20  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = 'PhotonJet'
#Use your own username instead of the "lhx". Keep branch tag in the directory name, e.g., PHYS14_720_Dec23_2014.
config.Data.outLFNDirBase = '/store/user/msun/'

config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
