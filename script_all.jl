using Logging

include("bpc.jl")
include("data.jl")

println("starting with $(Threads.nthreads()) threads")

do_compile_run = false
do_compile_run = true

if do_compile_run
    # first run to compile
    println("DOING COMPILE RUN\n\n")
    J, w, E, W = read_file("test/BPWC_0_6_8.txt") 
    solution, z = solve_bpc(J, E, w, W, verbose=0, run_ffd=true, max_iter=10000)
end

# General log
general_log_io = open("log.txt", "a")
logger = SimpleLogger(general_log_io)
global_logger(logger) # set as default logger

# csv with time and details file
csv_table = open("time.csv", "a")


# Files
folder_path = "test"
# folder_path = "Elhedhli"
instances = [joinpath(folder_path, i) for i in readdir(joinpath("instances", folder_path))]
mkpath(joinpath("logs", folder_path))


for file_path in instances

    if isfile(joinpath("logs/", file_path))
        println("instance at $(file_path) already done")
        continue
    else
        println("solving instance at $(file_path)")
    end

    # instance log
    global LOG_IO = open(joinpath("logs/", file_path), "a")
    
    # read instance
    J, w, E, W = read_file(file_path)
    
    # solve and get time elapsed
    println(LOG_IO, "J: $(J)\nw: $(w)\nW: $(W)\nE: $(E)")
    passed_time = @elapsed solution, z, is_optimal = solve_bpc(J, E, w, W, time_limit=600, verbose=0, run_ffd=true, max_iter=10000)
    println(LOG_IO, "solution: $(solution)\nz: $(z)")
    
    # register at general log
    @info z solution "$(passed_time) seconds"
    # println("$(passed_time) seconds")
    
    # save csv
    instance_name = basename(file_path)[1:end-4]
    println(csv_table, "$(instance_name), $(Int(is_optimal)), $(passed_time)")

    flush(LOG_IO)
    flush(csv_table)
    close(LOG_IO)
end

# Close the file
close(csv_table)
close(general_log_io)