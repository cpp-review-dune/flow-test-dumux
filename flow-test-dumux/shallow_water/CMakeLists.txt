add_input_file_links()
add_executable("shallow_water" shallow_water.cc)
target_link_dune_default_libraries("shallow_water")