include("data.jl")



sol = [[4, 6, 13, 30, 33], [1, 5, 7, 32, 34, 51], [16, 18, 21, 43], [42, 46, 50, 54, 56], [3, 23, 36, 47, 57], [17, 19, 28, 52], [12, 15, 35, 41, 48, 60], [44, 45, 53, 55, 58], [8, 10, 14, 20, 25], [2, 22, 24, 27, 39, 59], [29, 31, 37, 40, 49], [9, 11, 26, 38]]
instance = "test/BPWC_0_6_8.txt"

J, w, E, W = read_file(instance) 

println("W: $(W)")
println("w: $(w)")
println("E: $(E)\n\n")

println("sol:")
for i in sol
    println("\t$(i)")
end

check_solution_viability(sol, J, w, E, W)

println("\n\n")
println([w[j] for j in sol[2]])