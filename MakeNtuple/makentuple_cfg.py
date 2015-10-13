import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      '74X_dataRun2_v2',
      #'74X_dataRun2_Prompt_v1',
      #'MCRUN2_74_V9',
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

options.register( 'mfTag',
      'PAT',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "for MET filters"
      )

options.parseArguments()

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("METSigTuning.MakeNtuple.makentuple_cfi")
process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.GlobalTag.globaltag = ( options.globalTag )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
#process.options = cms.untracked.PSet(
#      allowUnscheduled = cms.untracked.bool(True),
#      wantSummary = cms.untracked.bool(False)
#      )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'root://xrootd.unl.edu//store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/162/00000/F6A6BB6F-4227-E511-BAAF-02163E014343.root'
       #'/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/162/00000/F6A6BB6F-4227-E511-BAAF-02163E014343.root'
       #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/10000/00BA30CE-9001-E511-AA08-0025905A60D0.root'
       #'/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/12284DB9-4227-E511-A438-02163E013674.root'
       '/store/data/Run2015B/DoubleMuon/MINIAOD/17Jul2015-v1/30000/3226AF16-C22E-E511-B9D7-0025905A613C.root'
       #'root://eoscms.cern.ch//store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/252/00000/263D331F-AF27-E511-969B-02163E012627.root'
    )
)

# re-apply JECs
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
      src = cms.InputTag("slimmedJets"),
      levels = ['L1FastJet', 
         'L2Relative', 
         'L3Absolute'],
      payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
      jetSource = cms.InputTag("slimmedJets"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
      )

# re-apply MET corrections
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#uncertainty file
jecUncertaintyFile="PhysicsTools/PatUtils/data/Summer15_50nsV4_DATA_UncertaintySources_AK4PFchs.txt"

#-----------------  switch off genMatching
#runMetCorAndUncFromMiniAOD(process,
#      isData = not options.runOnMCo,
#      jecUncFile=jecUncertaintyFile,
#      )

# produce MET
from RecoMET.METProducers.PFMET_cfi import pfMet
metType = "PF"
postfix = "RERUN"
setattr(process, "pfMet"+postfix, pfMet.clone() )
getattr(process, "pfMet"+postfix).src = cms.InputTag('packedPFCandidates')
getattr(process, "pfMet"+postfix).calculateSignificance = False

import PhysicsTools.PatAlgos.tools.helpers as configtools
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
setattr(process, 'pat'+metType+'Met'+postfix, getattr(process,'patPFMet' ).clone() )
from PhysicsTools.PatUtils.patPFMETCorrections_cff import patPFMet

configtools.cloneProcessingSnippet(process, getattr(process,"producePatPFMETCorrections"), postfix)
configtools.cloneProcessingSnippet(process, getattr(process,"patPFMetTxyCorrSequence"), postfix)

getattr(process, "patPFMet"+postfix).metSource = cms.InputTag("pfMet"+postfix)
getattr(process, "patPFMet"+postfix).addGenMET  = False

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
setattr(process, "selectedPatJets"+postfix, selectedPatJets.clone() )

getattr(process,"selectedPatJets"+postfix).src = cms.InputTag("patJets"+postfix)
getattr(process,"selectedPatJets"+postfix).cut = cms.string("pt > 10")

from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
setattr(process, "slimmedMETs"+postfix, slimmedMETs.clone() )
getattr(process,"slimmedMETs"+postfix).src = cms.InputTag("patPFMetT1"+postfix)

# produce MET corrections
patMetCorrectionSequence = cms.Sequence()
patMetModuleSequence = cms.Sequence()
metModName = "pat"+metType+"Met"+postfix
correctionLevel=["T1","Txy"]

corNames = { #not really needed but in case we have changes in the future....
      "T0":"T0pc",
      "T1":"T1",
      "T2":"T2",
      "Txy":"Txy",
      "Smear":"Smear",
      }

corModNames = {
      "T0": "patPFMetT0CorrSequence"+postfix,
      "T1": "patPFMetT1T2CorrSequence"+postfix,
      "T2": "patPFMetT2CorrSequence"+postfix,
      "Txy": "patPFMetTxyCorrSequence"+postfix,
      "Smear": "patPFMetSmearCorrSequence"+postfix,
      "T2Smear": "patPFMetT2SmearCorrSequence"+postfix
      }

configtools.cloneProcessingSnippet(process, process.patPFMetT0CorrSequence, postfix)
configtools.cloneProcessingSnippet(process, process.patPFMetT1T2CorrSequence, postfix)
configtools.cloneProcessingSnippet(process, process.patPFMetT2CorrSequence, postfix)
configtools.cloneProcessingSnippet(process, process.patPFMetSmearCorrSequence, postfix)
configtools.cloneProcessingSnippet(process, process.patPFMetT2SmearCorrSequence, postfix)

corModules = {}
for mod in corModNames.keys():
   corModules[mod] = getattr(process, corModNames[mod] )

corTags = {
   "T0":cms.InputTag('patPFMetT0Corr'+postfix),
   "T1":cms.InputTag('patPFMetT1T2Corr'+postfix, 'type1'),
   "T2":cms.InputTag('patPFMetT2Corr'+postfix,   'type2'),
   "Txy": cms.InputTag('patPFMetTxyCorr'+postfix),
   "Smear":cms.InputTag('patPFMetSmearCorr'+postfix, 'type1'),
   "Smear":cms.InputTag('patPFMetT1T2SmearCorr'+postfix, 'type1'),
   "T2Smear":cms.InputTag('patPFMetT2SmearCorr'+postfix, 'type2') 
   }

corScheme=""
corrections = []
correctionSequence = []
for cor in correctionLevel:
   corScheme += corNames[cor]
   corrections.append(corTags[cor])
   correctionSequence.append(corModules[cor])

#Txy parameter tuning
import PhysicsTools.PatUtils.patPFMETCorrections_cff as metCors
if "Txy" in correctionLevel:
   getattr(process, "patPFMetTxyCorr"+postfix).parameters = metCors.patMultPhiCorrParams_T1Txy_50ns

#create the main MET producer
metModName = "pat"+metType+"Met"+corScheme+postfix
corMetProducer = cms.EDProducer("CorrectedPATMETProducer",
      src = cms.InputTag('pat'+metType+'Met' + postfix),
      srcCorrections = cms.VInputTag(corrections)
      )
sequenceName="patMetCorrectionSequence"
for corModule in correctionSequence:
   patMetCorrectionSequence += corModule
   setattr(process, sequenceName+postfix, patMetCorrectionSequence)

setattr(process,metModName, corMetProducer)

setattr(process, "patMetCorrectionSequence"+postfix, patMetCorrectionSequence)
setattr(process, "patMetModuleSequence"+postfix, patMetModuleSequence)

#metModuleSequence += getattr(process, metModName)
jetCollectionUnskimmed = cms.InputTag('patJets')
getattr(process,"patPFMetT1T2Corr").src = jetCollectionUnskimmed
getattr(process,"patPFMetT2Corr").src = jetCollectionUnskimmed

getattr(process,'selectedPatJetsForMetT1T2Corr'+postfix).src = cms.InputTag('slimmedJets')
getattr(process,'patPFMetTxyCorr'+postfix).vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
getattr(process,'patPFMetTxyCorr'+postfix).srcPFlow = cms.InputTag('packedPFCandidates','')

process.test = cms.EDAnalyzer('MakeNtuple',
      src = cms.InputTag("packedPFCandidates"),
      #jets = cms.InputTag("slimmedJets"),
      jets = cms.InputTag('patJetsReapplyJEC'),
      leptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      #met = cms.InputTag("slimmedMETs"),
      met = cms.InputTag("patPFMetT1TxyRERUN"),
      muons = cms.InputTag("slimmedMuons"),
      vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      runOnMC = cms.untracked.bool(options.runOnMC),
      addPileupInfo = cms.InputTag("slimmedAddPileupInfo")
)

# trigger filter                
trigger_paths = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*']
trigger_pattern = [path+"*" for path in trigger_paths]
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )

#
## MET filters
#

# CSC and ee badSC filters
process.CSCBeamHaloFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_CSCTightHaloFilter'),
      throw=False
      )
process.eeBadScFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_eeBadScFilter'),
      throw=False
      )

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
      inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
      reverseDecision = cms.bool(False)
      )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
      vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
      minimumNDOF = cms.uint32(4),
      maxAbsZ = cms.double(24),
      maxd0 = cms.double(2)
      )

process.p = cms.Path(
      process.triggerSelection *
      process.CSCBeamHaloFilter *
      process.eeBadScFilter *
      process.HBHENoiseFilterResultProducer * #produces HBHE bools
      process.ApplyBaselineHBHENoiseFilter *  #reject events based 
      process.primaryVertexFilter *
      process.patJetCorrFactorsReapplyJEC *
      process.patJetsReapplyJEC *
      getattr(process, "pfMet"+postfix) *
      getattr(process, 'pat'+metType+'Met'+postfix) *
      getattr(process, 'patMetCorrectionSequence'+postfix) *
      getattr(process, metModName) *
      process.test
      )

#process.out = cms.OutputModule( "PoolOutputModule"
#      , fileName = cms.untracked.string( "test.root" ),
#      SelectEvents = cms.untracked.PSet(
#         SelectEvents = cms.vstring('p')
#         ),
#      )

#process.outpath = cms.EndPath( process.out )

