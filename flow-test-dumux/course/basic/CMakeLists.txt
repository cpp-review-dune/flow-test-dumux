add_input_file_links()

add_executable("2pmain" 2pmain.cc)
target_link_dune_default_libraries("2pmain")

add_executable("2p2cmain" 2p2cmain.cc)
target_link_dune_default_libraries("2p2cmain")

add_executable("2pnimain" 2pnimain.cc)
target_link_dune_default_libraries("2pnimain")