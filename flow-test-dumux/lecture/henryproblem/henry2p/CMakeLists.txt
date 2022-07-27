add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME henry2p
              TIMEOUT 1800
              SOURCES henry2p.cc
              COMMAND ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/lecture/references/Henry2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/Henry2p-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/henry2p -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/henry2p.input -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/henry.dgf -TimeLoop.TEnd 1e9")
