add_input_file_links()
add_executable("lens" lens.cc)
target_compile_definitions("lens" PRIVATE TYPETAG=RichardsLensCC PUBLIC DIFFMETHOD=DiffMethod::numeric)
target_link_dune_default_libraries("lens")