add_executable(htt-tuples htt-tuples.cpp)
target_link_libraries(htt-tuples htt)

add_executable(reco-inputs reco-inputs.cpp)
target_link_libraries(reco-inputs htt)

add_executable(reco-performance reco-performance.cpp)
target_link_libraries(reco-performance htt)

add_executable(htt-tuples-gen htt-tuples-gen.cpp)
target_link_libraries(htt-tuples-gen htt)

add_executable(mtt-hists mtt-hists.cpp)
target_link_libraries(mtt-hists htt ${Boost_PROGRAM_OPTIONS_LIBRARY})
