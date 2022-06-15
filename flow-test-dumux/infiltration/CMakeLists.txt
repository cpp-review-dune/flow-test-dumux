add_input_file_links()
dune_symlink_to_source_files(FILES "initialsolutioncc.txt")
add_executable("infiltration" infiltration.cc)
target_link_dune_default_libraries("infiltration")