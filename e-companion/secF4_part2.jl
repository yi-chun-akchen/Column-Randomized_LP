using JuMP, Gurobi

W = 200
w = [3,5,7,10,17,22,30,50]
m = 8
b = [1200,1000,1000,400,500,400,600,200]

function random_a(nDemands,w,W)
	a = zeros(nDemands)
	remaining_width = W
	valid_width_indices = find(w .<= remaining_width)

	while (!isempty(valid_width_indices))
		i = rand(valid_width_indices)
		a[i] += 1
		remaining_width -= w[i]
		valid_width_indices = find(w .<= remaining_width)
	end

	return a
end


srand(1061) 
nSim = 20000
K = 100

pattern_arrays = Array{Array{Int64,1},1}[]
x_arrays = []
obj_arrays = []

for s in 1:nSim

	patterns = zeros(m,K)

	for k in 1:K
		patterns[:,k] = random_a(m,w,W) 
	end

	mod = Model(solver = GurobiSolver(Method = 0))

	@variable(mod, x[1:K] >= 0)

	for i in 1:m
		@constraint(mod, sum(patterns[i,j] * x[j] for j in 1:K) >= b[i])
	end 

	setObjective(mod, :Min, sum(x[j] for j in 1:K))

	solve(mod)

	x_val = getValue(x)

	objval = getObjectiveValue(mod)

	nz_patterns = find(x_val .> 0)


	push!( pattern_arrays, [patterns[:,k] for k in nz_patterns])
	push!( x_arrays, x_val[nz_patterns])
	push!( obj_arrays, objval)
end

pattern_arrays_0 = deepcopy(pattern_arrays)
obj_arrays_0 = deepcopy(obj_arrays)
x_arrays_0 = deepcopy(x_arrays)

eps_target_array = [Inf,20.0,10.0,5.0,2.0, 1.0, 1e-8]


outcsvfilepath = string("CS_explo_out.csv");
outcsvhandle = open(outcsvfilepath, "w")
row_string = string("eps_target", ",",
					"numSolutions", ",",
					"numColumns", ",",
					"numUniqueColumns", ",",
					"numUniqueBases", ",",
					"averageIncidence", ",",
					"maximumIncidence", ",",
					"numColumnsWithIncidence1")
row_string = string(row_string, "\n");
print(outcsvhandle, row_string);
flush(outcsvhandle);

for eps_target in eps_target_array

	println("eps = ", eps_target)

	within_eps = find(obj_arrays_0 .<= 324.5 + eps_target )
	pattern_arrays = pattern_arrays_0[within_eps]
	obj_arrays = obj_arrays_0[within_eps]
	x_arrays = x_arrays_0[within_eps]


	# What is the number of overall solutions (including duplicates)?
	numSolutions = length(pattern_arrays)
	println("Number of solutions in pattern_arrays: ", numSolutions)

	# What is the number of patterns (including duplicates)?
	numColumns = sum(map(length, pattern_arrays))
	println("Number of patterns in pattern_arrays: ", numColumns)

	# What is the set of unique columns?
	# temp1 = map(t -> Set(t), pattern_arrays)
	unique_column_set = union( pattern_arrays... )
	numUniqueColumns = length(unique_column_set)
	println("Number of unique columns: ", numUniqueColumns)

	# What is the set of unique bases?
	temp1 = map(t -> Set(t), pattern_arrays)
	temp2 = unique(temp1)
	unique_basis_set = temp2
	numUniqueBases = length(unique_basis_set)
	println("Number of unique bases: ", numUniqueBases)


	# What is the set of unique bases, when we intentionally screw things up.
	# temp0 = [pattern_arrays[1][randperm(length(pattern_arrays[1]))] for i in 1:nSim]
	# temp1 = map(t -> Set(t), temp0)
	# temp2 = Set(temp1)
	# unique_basis_set_2 = temp2
	# println("Number of unique bases (sanity check): ", length(unique_basis_set_2))

	# How many time does each unique column appear in each unique basis? 
	count_by_unique_column = map( j -> sum(map(S -> j in S, unique_basis_set)), unique_column_set)
	averageIncidence = mean(count_by_unique_column)
	maximumIncidence = maximum(count_by_unique_column)

	numColumnsWithIncidence1 = sum( count_by_unique_column .== 1)
	println("Average incidence: ", averageIncidence)
	println("Maximum incidence: ", maximumIncidence)
	println("Columns appearing only once: ", numColumnsWithIncidence1)


	row_string = string(eps_target, ",",
					numSolutions, ",",
					numColumns, ",",
					numUniqueColumns, ",",
					numUniqueBases, ",",
					averageIncidence, ",",
					maximumIncidence, ",",
					numColumnsWithIncidence1)
	row_string = string(row_string, "\n");
	print(outcsvhandle, row_string);
	flush(outcsvhandle);
end

close(outcsvhandle)

num_optimal = sum( obj_arrays_0 .<= 324.5 + 1e-8)
println("Num. optimal solutions: ", num_optimal)

which_optimal = find( obj_arrays_0 .<= 324.5 + 1e-8)

for t in 1:5
	for j in 1:length(pattern_arrays_0[which_optimal[t]])
		println("\\Ab_", j, " = ", pattern_arrays_0[which_optimal[t]][j], " && x_", j, " = ", round(x_arrays_0[which_optimal[t]][j], 4), "\\\\")
	end 
	println()
end 
# @show pattern_arrays_0[which_optimal[2]]
# @show pattern_arrays_0[which_optimal[3]]
# @show pattern_arrays_0[which_optimal[4]]
# @show pattern_arrays_0[which_optimal[5]]

println("Test for overlap:")
for t in 1:4
	for t2 in (t+1):5
		println(t, ",", t2, ": ", intersect(pattern_arrays_0[which_optimal[t]], pattern_arrays_0[which_optimal[t2]]))
	end 
end 

