add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable("henry2p" henry2p.cc)
target_link_dune_default_libraries("henry2p")

# henry2p -ParameterFile henry2p.input -Grid.File grids/henry.dgf -TimeLoop.TEnd 1e9