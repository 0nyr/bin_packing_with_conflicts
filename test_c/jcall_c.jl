a = 15
b = 10
c = @ccall "./pure_c/advmath.so".add(a::Cint, b::Cint)::Cint
println("Julia here, I have sent $(a) and $(b) and have received $(c)")
