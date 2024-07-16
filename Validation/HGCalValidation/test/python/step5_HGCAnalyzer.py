import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('HGCANA',Phase2C17I13M9)

process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D110_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi')
process.load('SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

options = VarParsing.VarParsing ('analysis')

# options.inputFiles= 'file:/uscms/home/ahussain/nobackup/YOURWORKINGAREA/CMSSW_12_1_0/src/condor/step3.root'
# options.inputFiles= 'file:/uscms_data/d2/kunori/hgc21/dpg/CMSSW_12_4_0_pre2/src/skclue3d/clue3d_step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.root'
options.inputFiles= 'file:step3.root'
options.outputFile = 'step5_nt.root'
# options.outputFile = 'test.root'
options.maxEvents = -1

options.parseArguments()

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
                                   )

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.recHitTask = cms.Task(process.recHitMapProducer,
                              process.lcAssocByEnergyScoreProducer,
)
process.hgcAnalyzer = cms.EDAnalyzer('HGCAnalyzer',
				    setZside = cms.untracked.int32(1),  # 1 zplus, -1 zminus, 0 both
            trackstersmrg = cms.InputTag('ticlCandidate')
				    )

process.p = cms.Path(process.hgcAnalyzer, process.recHitTask)
#process.p = cms.EndPath(process.p + process.hgcAnalyzer)
