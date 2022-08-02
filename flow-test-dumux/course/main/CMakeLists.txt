add_input_file_links()

add_executable("exercise1pamain" exercise1pamain.cc)
target_compile_definitions("exercise1pamain" PRIVATE TYPETAG=OnePIncompressible)
target_link_dune_default_libraries("exercise1pamain")

add_executable("exercise1pbmain" exercise1pbmain.cc)
target_link_dune_default_libraries("exercise1pbmain")

add_executable("exercise1pcmain" exercise1pcmain.cc)
target_link_dune_default_libraries("exercise1pcmain")