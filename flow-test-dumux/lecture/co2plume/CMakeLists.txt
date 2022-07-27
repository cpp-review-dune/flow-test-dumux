add_input_file_links()

add_executable("co2plumeshapeexercise" co2plumeshapeexercise.cc)
target_compile_definitions("co2plumeshapeexercise" PRIVATE TYPETAG=PlumeShapeBoxTypeTag)
target_link_dune_default_libraries("co2plumeshapeexercise")
# co2plumeshapeexercise -TimeLoop.TEnd 7e6 -TimeLoop.EpisodeLength 1e5  -BoundaryConditions.InjectionRate 0.01
