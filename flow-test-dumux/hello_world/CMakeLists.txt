add_input_file_links()
# https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/test/porousmediumflow/1p/isothermal/CMakeLists.txt
add_executable("test_1p_tpfa" test_1p_tpfa.cc)
target_compile_definitions("test_1p_tpfa" PRIVATE TYPETAG=OnePTestCCTpfa)
target_link_dune_default_libraries("test_1p_tpfa")