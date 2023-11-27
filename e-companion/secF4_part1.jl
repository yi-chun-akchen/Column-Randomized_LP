using JuMP, Gurobi

W = 200
w = [3,5,7,10,17,22,30,50]
m = 8
b = [1200,1000,1000,400,500,400,600,200]

#### DO EXACT CG HERE ####
m_dual = Model(solver = GurobiSolver(OutputFlag = 0))
@variable(m_dual, p[1:m] >= 0)
setObjective(m_dual, :Max, sum(b[i] * p[i] for i in 1:m))

pattern_list = Array{Int64,1}[]
first_pattern = ones(Int64,m)
push!(pattern_list, first_pattern)
@constraint(m_dual, sum(first_pattern[i] * p[i] for i in 1:m) <= 1)

solve(m_dual)



println("Solving subproblem!!!")
m_sub = Model(solver = GurobiSolver(OutputFlag = 0))
@variable(m_sub, a[1:m] >= 0, Int)
@constraint(m_sub, sum(w[i] * a[i] for i in 1:m) <= W)
function separation(p_val)
	setObjective(m_sub, :Max, sum(a[i] * p_val[i] for i in 1:m))
	solve(m_sub)

	a_val = getValue(a)
	sub_obj = getObjectiveValue(m_sub)

	return a_val, sub_obj
end

p_val = getValue(p)
a_val, sub_obj = separation(p_val)

dual_obj = getObjectiveValue(m_dual)

CG_iter = 0 
println("Iteration ", CG_iter, ": ", dual_obj)

tol_eps = 1e-5 
isOptimal = (sub_obj <= (1 + tol_eps))
println("sub_obj = ", sub_obj)
@show a_val
@show p_val

while (!isOptimal)
	CG_iter += 1
	push!(pattern_list, a_val)

	@constraint(m_dual, sum(a_val[i] * p[i] for i in 1:m) <= 1)

	solve(m_dual)

	dual_obj = getObjectiveValue(m_dual)
	p_val = getValue(p)

	println("Iteration ", CG_iter, ": ", dual_obj)


	a_val, sub_obj = separation(p_val)
	tol_eps = 1e-5 
	isOptimal = (sub_obj <= (1 + tol_eps))
	println("sub_obj = ", sub_obj)


end 


m_primal = Model(solver = GurobiSolver())

n = length(pattern_list)
@variable(m_primal, x[1:n] >= 0)

for i in 1:m 
	@constraint(m_primal, sum(pattern_list[j][i] * x[j] for j in 1:n) >= b[i])
end 

setObjective(m_primal, :Min, sum(x[j] for j in 1:n))

solve(m_primal)

x_val = getValue(x) 

nzind = x_val .> 0 

pattern_list = pattern_list[ nzind ]
x_val = x_val[nzind]


println()
for j in 1:length(x_val)
	println("\\Ab_", j , " = ", pattern_list[j], " && x_", j, " = ", round(x_val[j], 4), " \\\\")
end 


