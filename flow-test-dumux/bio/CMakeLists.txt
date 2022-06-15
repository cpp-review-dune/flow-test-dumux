add_input_file_links()
dune_symlink_to_source_files(FILES injection_checkpoints.dat injection_type.dat)
add_executable("bio" bio.cc)
target_link_dune_default_libraries("bio")