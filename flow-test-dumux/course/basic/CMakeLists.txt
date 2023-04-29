add_input_file_links()

# the immiscible two-phase simulation program
add_executable("2pmain" 2pmain.cc)
target_link_dune_default_libraries("2pmain")

# the compositional two-phase two-component simulation program
add_executable("2p2cmain" 2p2cmain.cc)
target_link_dune_default_libraries("2p2cmain")

# here, add the two-phase non-isothermal simulation program
add_executable("2pnimain" 2pnimain.cc)
target_link_dune_default_libraries("2pnimain")