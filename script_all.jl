using Logging

include("bpc.jl")
include("data.jl")

# first run to compile
J, w, E, W = read_file("BPWC_0_6_8.txt") 
solution, z = solve_bpc(J, E, w, W, verbose=0, run_ffd=true, max_iter=10000)



# General log
general_log_io = open("log.txt", "w+")
logger = SimpleLogger(general_log_io)
global_logger(logger) # set as default logger

# csv with time and details file
csv_table = open("time.csv", "w+")

# Files
# file_name = "toy"
# file_name = "toy2"
file_name = "BPWC_0_6_8"

# instance log
global LOG_IO = open("logs/$(file_name).txt", "w+")

# read instance
J, w, E, W = read_file("$(file_name).txt")


# solve and get time elapsed
println(LOG_IO, "J: $(J)\nw: $(w)\nW: $(W)\nE: $(E)")
passed_time = @elapsed solution, z = solve_bpc(J, E, w, W, verbose=0, run_ffd=true, max_iter=10000)
println(LOG_IO, "solution: $(solution)\nz: $(z)")

# register at general log
@info z solution "$(passed_time) seconds"

# save csv
println(csv_table, "$(file_name), $(passed_time)")



# Close the file
close(general_log_io)
