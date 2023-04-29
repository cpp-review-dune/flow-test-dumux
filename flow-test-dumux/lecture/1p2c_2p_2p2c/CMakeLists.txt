add_input_file_links()
add_executable("lens2pexercise3" lens2pexercise3.cc)
target_link_dune_default_libraries("lens2pexercise3")

# exercise3.input -Problem.Name lens-2p
# exercise3.input -Problem.Name lens-2p2c -SpatialParams.Temperature 313.15
# exercise3.input -Problem.Name lens-1p2c