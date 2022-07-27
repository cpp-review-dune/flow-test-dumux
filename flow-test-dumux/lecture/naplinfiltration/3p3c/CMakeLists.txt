add_input_file_links()

add_executable("naplinfiltration3p3c" naplinfiltration3p3c.cc)
target_link_dune_default_libraries("naplinfiltration3p3c")

# --zeroThreshold {"Sn":5e-7,"x^mesitylene_g":1e-8,"x^mesitylene_n":1e-5,"x^mesitylene_w":5e-11}
# naplinfiltration3p3c -TimeLoop.TEnd 604800 -TimeLoop.EpisodeLength 302400
