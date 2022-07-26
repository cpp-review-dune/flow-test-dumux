add_input_file_links()

add_executable("lens1p2cexercise1" lens1p2cexercise1.cc)
target_link_dune_default_libraries("lens1p2cexercise1")

add_executable("lens2pexercise1" lens2pexercise1.cc)
target_link_dune_default_libraries("lens2pexercise1")

# exercise1.input -Problem.Name lens-1p2c --zeroThreshold {"x^N2_liq":1e-18,"X^N2_liq":1e-18}
# exercise1.input -Problem.Name lens-2p