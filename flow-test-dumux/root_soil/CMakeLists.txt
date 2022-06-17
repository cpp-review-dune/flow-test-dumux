add_input_file_links()
dune_symlink_to_source_files(FILES "tinyroot.dgf" "tinyroot_3d.msh")

add_executable("root_soil" root_soil.cc)
target_compile_definitions("root_soil" PRIVATE SOILTYPETAG=SoilCC)
target_link_dune_default_libraries("root_soil")

add_executable("root_soil_proj" root_soil.cc)
target_compile_definitions("root_soil_proj" PRIVATE SOILTYPETAG=SoilBox COUPLINGMODE=Projection SOILGRID=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)
target_link_dune_default_libraries("root_soil_proj")