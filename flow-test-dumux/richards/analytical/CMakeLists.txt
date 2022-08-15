add_input_file_links()
add_executable("analytical" analytical.cc)
target_link_dune_default_libraries("analytical")