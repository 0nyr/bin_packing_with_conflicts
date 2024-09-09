using Logging

include("bpc.jl")
include("data.jl")

println("starting with $(Threads.nthreads()) threads")

do_compile_run = false
# do_compile_run = true

verbose=0
verbose=1

if do_compile_run
    # first run to compile
    println("DOING COMPILE RUN\n\n")
    J, w, E, W = read_file("test/BPWC_0_6_8.txt") 
    # J, w, E, W = read_file("test/BPWC_2_7_2.txt") 
    # J, w, E, W = read_file("biginterval/BPWC_1_1_3.txt") 
    # J, w, E, W = read_file("Elhedhli/BPWC_1_5_5.txt") 
    solution, z = solve_bpc(J, E, w, W, verbose=verbose, run_ffd=true, max_iter=10000)
end

# General log

# instance log
# global LOG_IO = open(joinpath("logs/", file_path), "a")

instance = "test/BPWC_0_6_8.txt"
println("Running instance $(instance)")

J, w, E, W = read_file(instance) 

# solve and get time elapsed
println("J: $(J)\nw: $(w)\nW: $(W)\nE: $(E)")
passed_time = @elapsed solution, z, is_optimal = solve_bpc(J, E, w, W, time_limit=600, verbose=verbose, run_ffd=true, max_iter=10000)
println("solution: $(solution)\nz: $(z)\noptimal: $(is_optimal)\ntime: $(passed_time)")

