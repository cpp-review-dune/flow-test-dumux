add_input_file_links()

add_executable("naplinfiltration3p" naplinfiltration3p.cc)
target_link_dune_default_libraries("naplinfiltration3p")
# naplinfiltration3p -TimeLoop.TEnd 604800 -TimeLoop.EpisodeLength 302400