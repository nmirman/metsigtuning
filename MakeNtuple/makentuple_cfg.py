import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      #'74X_dataRun2_v2',
      'MCRUN2_74_V9',
      #'GR_P_V56',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "CMS Global Tag"
      )

options.register( 'runOnMC',
      True,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "mc or data"
      )

options.parseArguments()

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("METSigTuning.MakeNtuple.makentuple_cfi")
process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = ( options.globalTag )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'root://xrootd.unl.edu//store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/162/00000/F6A6BB6F-4227-E511-BAAF-02163E014343.root'
       #'/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/162/00000/F6A6BB6F-4227-E511-BAAF-02163E014343.root'
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/10000/00BA30CE-9001-E511-AA08-0025905A60D0.root'
    )
)

# correct jets
#process.ak4PFchsCorr = cms.ESSource(
#      'LXXXCorrectionService',
#      #era = cms.string('Jec11V12'),
#      #section   = cms.string(''),
#      #level     = cms.string('L2Relative'),
#      # the above 3 elements are needed only when the service is initialized from local txt files
#      algorithm = cms.string('AK4PFchs'),
#      # the 'algorithm' tag is also the name of the DB payload
#      useCondDB = cms.untracked.bool(True)
#      )

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *

process.ak4PFCHSJetsCorrected = cms.EDProducer('PFJetCorrectionProducer',
      #src         = cms.InputTag('ak4PFCHSJets'),
      #correctors  = cms.vstring('AK4PF')
      src = cms.InputTag('ak4PFJets'),
      correctors = cms.vstring("ak4PFL1FastL2L3") 
      )

#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

#from METSigTuning.MakeNtuple.JetCorrection_cff import *
#process.load('METSigTuning.MakeNtuple.JetCorrection_cff')

## met corrections setup copied from
## https://github.com/TaiSakuma/WorkBookMet/blob/master/corrMet_cfg.py
##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")

process.corrPfMetType1.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3Corrector")
if not options.runOnMC:
  process.corrPfMetType1.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3ResidualCorrector")

##____________________________________________________________________________||
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

##____________________________________________________________________________||
process.metcorr = cms.Sequence(
      process.correctionTermsPfMetType1Type2 +
      process.pfMetT1
      )

process.test = cms.EDAnalyzer('MakeNtuple',
      src = cms.InputTag("particleFlow"),
      #jets = cms.InputTag("ak4PFJetsCHS"),
      jets = cms.InputTag("ak4PFCHSJetsCorrected"),
      leptons = cms.VInputTag("gedGsfElectrons", "muons", "photons"),
      met = cms.InputTag("pfMetT1"),
      muons = cms.InputTag("muons"),
      vertices = cms.InputTag("offlinePrimaryVertices"),
      #pfjetCorrector = cms.untracked.string("ak4PFchsL1FastL2L3")
)
if not options.runOnMC:
      process.demo.pfjetCorrector = 'ak4PFchsL1FastL2L3Residual'


process.p = cms.Path(
      process.metcorr *
      #process.ak4PFchsJetsL1FastL2L3 *
      process.ak4PFCHSJetsCorrected * 
      process.test
      )

#process.out = cms.OutputModule( "PoolOutputModule"
#      , fileName = cms.untracked.string( "test.root" ),
#      SelectEvents = cms.untracked.PSet(
#         SelectEvents = cms.vstring('p')
#         )
#      )
#process.outpath = cms.EndPath( process.out )
