include("modelA2D.jl")

seed = parse(Int,ARGS[1])
Random.seed!(seed)
const L = parse(Int,ARGS[2])
max = parse(Int,ARGS[3])

file = jldopen(filename, "r")
ϕ = Array(file["ϕ"])
close(file)

m² = -3.824f0
maxt = max*L^2

open("output_$L.dat","w") do io 
	for i in 0:maxt
		if i % (L^2) == 0
			(M,ϕk) = op(ϕ, L)
			Printf.@printf(io, "%i %f", i, M)
			for kx in 1:L
				Printf.@printf(io, " %f %f", real(ϕk[kx]), imag(ϕk[kx]))
			end 
			Printf.@printf(io, "\n")
		end
		thermalize(m², ϕ, L, 1)
	end
end