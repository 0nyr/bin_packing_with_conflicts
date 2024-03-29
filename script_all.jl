using Logging

include("bpc.jl")
include("data.jl")

# first run to compile
J, w, E, W = read_file("test/BPWC_0_6_8.txt") 
solution, z = solve_bpc(J, E, w, W, verbose=0, run_ffd=true, max_iter=10000)

# General log
general_log_io = open("log.txt", "w+")
logger = SimpleLogger(general_log_io)
global_logger(logger) # set as default logger

# csv with time and details file
csv_table = open("time.csv", "w+")

# Files

folder_path = "test"
instances = [joinpath(folder_path, i) for i in readdir(joinpath("instances", folder_path))]
mkpath(joinpath("logs", folder_path))

for file_path in instances

    # instance log
    global LOG_IO = open(joinpath("logs/", file_path), "w+")
    
    # read instance
    J, w, E, W = read_file(file_path)
    
    # solve and get time elapsed
    println(LOG_IO, "J: $(J)\nw: $(w)\nW: $(W)\nE: $(E)")
    passed_time = @elapsed solution, z = solve_bpc(J, E, w, W, verbose=0, run_ffd=true, max_iter=10000)
    println(LOG_IO, "solution: $(solution)\nz: $(z)")
    
    # register at general log
    @info z solution "$(passed_time) seconds"
    
    # save csv
    instance_name = basename(file_path)[1:end-4]
    println(csv_table, "$(instance_name), $(passed_time)")
end

# Close the file
close(general_log_io)
