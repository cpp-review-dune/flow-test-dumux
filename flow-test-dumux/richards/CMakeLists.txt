add_input_file_links()
add_executable("richards" richards.cc)
target_link_dune_default_libraries("richards")