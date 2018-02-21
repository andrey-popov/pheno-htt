add_library(htt SHARED AnalysisPlugin.cpp DelphesReader.cpp LJetsSelection.cpp NuReco.cpp
    Plugin.cpp Processor.cpp TTReco.cpp TTRecoInputs.cpp TTRecoPerf.cpp VarWriter.cpp)
target_link_libraries(htt ${DELPHES_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY})
