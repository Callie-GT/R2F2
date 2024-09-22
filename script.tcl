open_project project_1
set_top flexible_float_multiply
#set_top templated_float_multiply

add_files flexible_floating_point.cpp
add_files templated_floating_point.cpp


open_solution "solution_1"
set_part {xc7z020clg484-1}
create_clock -period 10 -name default

config_compile -pipeline_loops 1

csynth_design

exit
