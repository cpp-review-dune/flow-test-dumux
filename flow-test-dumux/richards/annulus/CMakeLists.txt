add_input_file_links()
add_executable("main_stationary" main_stationary.cc)
target_link_dune_default_libraries("main_stationary")

#add_executable("main_instationary" main_instationary.cc)
#target_link_dune_default_libraries("main_instationary")