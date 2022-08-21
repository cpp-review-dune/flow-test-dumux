add_input_file_links()
dune_symlink_to_source_files(FILES run_and_plot_m21.py run_and_plot_m22.py)
add_executable("benchmarks" benchmarks.cc)
target_link_dune_default_libraries("benchmarks")