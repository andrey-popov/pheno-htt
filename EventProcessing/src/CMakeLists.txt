add_library(htt SHARED Plugin.cpp AnalysisPlugin.cpp Processor.cpp
    DelphesReaderBase.cpp DelphesReader.cpp DelphesReaderGen.cpp
    LJetsLHEFilter.cpp LJetsSelection.cpp
    NuReco.cpp TTReco.cpp TTRecoInputs.cpp TTRecoPerf.cpp
    SmearMttWriter.cpp SystMttHists.cpp VarWriter.cpp)
target_link_libraries(htt ${DELPHES_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY})
