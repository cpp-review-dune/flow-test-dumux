add_input_file_links()

add_executable("properties" properties.cc)
target_link_dune_default_libraries("properties")