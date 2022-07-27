add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable("henry1p2c" henry1p2c.cc)
target_link_dune_default_libraries("henry1p2c")

# henry1p2c -ParameterFile henry1p2c.input -Grid.File grids/henry.dgf -TimeLoop.TEnd 1e10