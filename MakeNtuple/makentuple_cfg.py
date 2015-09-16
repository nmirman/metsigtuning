import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      'MCRUN2_74',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "CMS Global Tag"
      )

options.parseArguments()

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("METSigTuning.MakeNtuple.makentuple_cfi")
process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = ( options.globalTag+'::All' )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/RECO/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/0002583F-CCDB-E411-91BB-003048FFCB6A.root'
       #'root://xrootd.unl.edu//store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/RECO/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00001/88B29EF6-9FDB-E411-8879-0025905A60BE.root'
    )
)

#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

#from METSigTuning.MakeNtuple.JetCorrection_cff import *
#process.load('METSigTuning.MakeNtuple.JetCorrection_cff')

## met corrections setup copied from
## https://github.com/TaiSakuma/WorkBookMet/blob/master/corrMet_cfg.py
##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")

process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
# process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")

process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
# process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

##____________________________________________________________________________||
process.metcorr = cms.Sequence(
      process.correctionTermsPfMetType1Type2 +
      process.correctionTermsPfMetType0RecoTrack +
      process.correctionTermsPfMetType0PFCandidate +
      process.correctionTermsPfMetShiftXY +
      process.correctionTermsCaloMet +
      process.caloMetT1 + 
      process.caloMetT1T2 + 
      process.pfMetT0rt +
      process.pfMetT0rtT1 +
      process.pfMetT0pc +
      process.pfMetT0pcT1 +
      process.pfMetT0rtTxy +
      process.pfMetT0rtT1Txy +
      process.pfMetT0pcTxy +
      process.pfMetT0pcT1Txy +
      process.pfMetT1 +
      process.pfMetT1Txy
      )

process.test = cms.EDAnalyzer('MakeNtuple',
      src = cms.InputTag("particleFlow"),
      jets = cms.InputTag("ak4PFJets"),
      #jets = cms.InputTag("ak4PFchsJetsL1FastL2L3"),
      leptons = cms.VInputTag("gedGsfElectrons", "muons", "photons"),
      met = cms.InputTag("pfMet"),
      muons = cms.InputTag("muons"),
      vertices = cms.InputTag("offlinePrimaryVertices"),
      pfjetCorrectorL1 = cms.untracked.string("ak4PFchsL1Fastjet"),
      pfjetCorrectorL123 = cms.untracked.string("ak4PFchsL1FastL2L3")
)
#if not options.runOnMC:
#      process.demo.pfjetCorrectorL123 = 'ak5PFL1FastL2L3Residual'


process.p = cms.Path(
      #process.metcorr *
      #process.ak4PFchsJetsL1FastL2L3 *
      process.test
      )
