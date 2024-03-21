include("bpc.jl")
include("data.jl")

# J, w, E, W = read_file("toy.txt")
J, w, E, W = read_file("toy2.txt")
# J, w, E, W = read_file("BPWC_0_6_8.txt") # currently leads to memory issues

solve_bpc(J, E, w, W, verbose=0, run_ffd=true)