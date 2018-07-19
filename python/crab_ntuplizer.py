from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'GMSB_ggtree'
config.General.workArea = 'GMSB_ggtree'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_mc_74X.py' 

config.section_("Data")
config.Data.outputPrimaryDataset = 'GMSB_MC'
config.Data.userInputFiles = ['/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_1.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_2.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_3.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_4.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_5.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_6.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_7.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_8.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_9.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_10.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_11.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_12.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_13.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_14.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_15.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_16.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_17.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_18.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_19.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_20.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_21.root',
		'/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_22.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_23.root',
		'/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_24.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_25.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_26.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_27.root',
		'/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_28.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_29.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_30.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_31.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_32.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_33.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_34.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_35.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_36.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_37.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_38.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_39.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_40.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_41.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_42.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_43.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_44.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_45.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_46.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_47.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_48.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_49.root',
		 '/store/user/msun/GMSB_MC/GMSBMC/151104_191538/0000/slhafragment_cff_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_50.root']
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
NJOBS = 5  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = 'GMSB_ggtree'
#Use your own username instead of the "lhx". Keep branch tag in the directory name, e.g., PHYS14_720_Dec23_2014.
config.Data.outLFNDirBase = '/store/user/msun/'

config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
