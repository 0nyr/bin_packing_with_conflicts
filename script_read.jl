include("bpc.jl")
include("data.jl")

# J, w, E, W = read_file("toy.txt")
# J, w, E, W = read_file("toy2.txt")
J, w, E, W = read_file("BPWC_0_6_8.txt") 

println("J: $(J)\nw: $(w)\nW: $(W)\nE: $(E)")

solution, z = solve_bpc(J, E, w, W, verbose=1, run_ffd=true)

println("solution: $(solution)\nz: $(z)")