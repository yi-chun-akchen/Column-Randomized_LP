using JuMP, Gurobi, Distributions, StatsBase, DataStructures, DataFrames, CSV

# -------------------------------------------------------- #
# Column Generation Approach for the Cutting Stock Problem #
# -------------------------------------------------------- #

function subproblem_cutting_stock(p,w,W)

    # This function is a suproblem solver (equation 13 in the paper)
    # for the column generation approach to the cutting stock problem.
    # p is the dual solution
    # w is the width for each demand
    # W is the width of the large rolls
    # See the definition in the beginning paragraph in Section 5

    M = size(p)[1]

    #subproblem = Model(solver = GurobiSolver(OutputFlag = 0,TimeLimit=120))
    subproblem = Model(solver = GurobiSolver(OutputFlag = 0))

    @variable(subproblem,a[1:M],Int)

    @constraint(subproblem, sum( w[ind_w] * a[ind_w]  for ind_w in 1:M) <= W)
    @constraint(subproblem, a[1:M] .>= 0)

    @objective(subproblem, Max, sum(p[ind_w] * a[ind_w] for ind_w in 1:M))

    solve(subproblem)

    return (getobjectivevalue(subproblem),getValue(a))

end



function dual_of_cutting_stock(b_demands,w_sizes,W,iter)

    # This function solves the cutting stock problem with the column generation approach.
    # We solve the dual of problem (12).
    # b_demands is the demand b_i
    # w_sizes is the width w_i
    # W is the width of the large rolls
    # iter is the upper bound of iterations for the column generation approach.
    # We set iter to be very high so that it doesn't affect the optimality while we can
    # make sure the function for sure terminates.

    t0 = time_ns()
    M = size(b_demands)[1]

    dual_problem = Model(solver = GurobiSolver(OutputFlag = 0))

    @variable(dual_problem,p[1:M])

    A_transpose = zeros(Float64,iter,M)

    # first constraint
    a_first = zeros(Float64,M); a_first[1] = 1
    @constraint(dual_problem, sum(a_first[ind_w] * p[ind_w] for ind_w in 1:M) <= 1)
    A_transpose[1,:] = a_first

    @objective(dual_problem, Max, sum(b_demands[ind_w] * p[ind_w] for ind_w in 1:M))

    solve(dual_problem)
    dual_obj = getobjectivevalue(dual_problem)

    traj = zeros(Float64,iter,2)
    num_add = 1
    for it_time = 1 : iter-1

        p_value = getvalue(p)

        subproblem_outcome = subproblem_cutting_stock(p_value,w_sizes,W)

        if subproblem_outcome[1] > 1.0
            num_add += 1
            a_new = subproblem_outcome[2]
            @constraint(dual_problem, sum(a_new[ind_w] * p[ind_w] for ind_w in 1:M) <= 1.0)
            A_transpose[num_add,:] = a_new

        else
            break
        end

        solve(dual_problem)
        dual_obj = getobjectivevalue(dual_problem)
        t1 = time_ns()
        traj[it_time,:] = [(t1-t0)/1.0e9,dual_obj]
        println(string(it_time,"th iteration with dual obj = ",dual_obj," with time ",(t1-t0)/1.0e9))
    end

    primal_problem = Model(solver = GurobiSolver(OutputFlag = 0))
    A = transpose(A_transpose[1:num_add,:])

    @variable(primal_problem,x[1:num_add] >= 0)
    @constraint(primal_problem, A * x .>= b_demands)
    @objective(primal_problem, Min,sum(x[ind_a] for ind_a in 1:num_add))

    solve(primal_problem)

    t1 = time_ns()
    traj[num_add,:] = [(t1-t0)/1.0e9,getobjectivevalue(primal_problem)]

    return (getobjectivevalue(primal_problem),A,traj[1:num_add,:])
end

# --------------------------- #
# Column Randomization Method #
# --------------------------- #

function sample_one_a(w_sizes,W)

    # We will sample a column according to Algorithm 5 in the paper.

    a_sampled = zeros(Float64,length(w_sizes))
    W_copy = copy(W)

    while W_copy > 0

        cand_set = findall(x -> x <= W_copy,w_sizes)
        if length(cand_set) < 1
            break
        else
            idx_sampled = rand(cand_set,1)[1]
            a_sampled[idx_sampled] += 1
            W_copy -= w_sizes[idx_sampled]
        end
    end

    return a_sampled
end

function sampling(b_demands,w_sizes,W,n_samples)

    # n_samples in the number of columns that you're going to sample,
    # which is K in the paper.

    tt0 = time_ns()
    num_add = 0
    M = size(w_sizes)[1]

    A_constraint = zeros(Float64,M,n_samples)

    for ind_sample = 1 : n_samples
        new_a = sample_one_a(w_sizes,W)
        A_constraint[:,ind_sample] = new_a
    end

    return A_constraint
end

function CRLP_with_sampled_columns(b_demands,w_sizes,W,sampled_A)

    # In this function, after we sample the constraint matrix as sampled_A,
    # we solve the corresponding column-randomized LP.

    CRLP = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(CRLP,x_crlp[1:size(sampled_A)[2]] >= 0)
    @constraint(CRLP, sampled_A * x_crlp .>= b_demands)
    @objective(CRLP, Min, sum(x_crlp[ind_a] for ind_a in 1:size(sampled_A)[2]))
    solve(CRLP)

    return (getobjectivevalue(CRLP),sampled_A)
end

function CRLP_primal(b_demands,w_sizes,W,n_samples)

    A_constraint = sampling(b_demands,w_sizes,W,n_samples)
    CRLP = CRLP_with_sampled_columns(b_demands,w_sizes,W,A_constraint)
    return CRLP
end



# --------------- #
# Reject Sampling #
# --------------- #

function reject_sampling_cutting_stock(W,w_sizes)

    # We sample a column by reject sampling
    # We return a column only if it is feasible to the constraint a^T w \leq W

    neg_feasibility = true
    while neg_feasibility
        a_column = zeros(Int64,length(w_sizes))
        for ind_type = 1 : length(w_sizes)
            max_num = floor(W / w_sizes[ind_type])
            a_column[ind_type] = rand(0:max_num)
        end
        if dot(a_column,w_sizes) <= W
            neg_feasibility = false
            return a_column
        end
    end
end


function reject_sampling_contraint_matrix(W,w_sizes,num_columns)
    # We sample a constraint matrix where each column is sampled by reject sampling
    sampled_matrix = zeros(Int64,(length(w_sizes),num_columns))
    for ind_col = 1 : num_columns
        sampled_matrix[:,ind_col] = reject_sampling_cutting_stock(W,w_sizes)
    end
    return sampled_matrix

end

# ----------------------------#
# Biased Incremental Sampling #
# ----------------------------#

function sample_one_a_BI(w_sizes,W,b_vec)

    # We sample according to Algorithm 10

    a_sampled = zeros(Float64,length(w_sizes))
    W_copy = copy(W)

    while W_copy > 0
        cand_set = findall(x -> x <= W_copy,w_sizes)
        if length(cand_set) < 1
            break
        else
            b_components = sqrt.(b_vec[cand_set])
            normalized_b_components = b_components / sum(b_components)
            idx_sampled = sample(cand_set,Weights(normalized_b_components))
            a_sampled[idx_sampled] += 1
            W_copy -= w_sizes[idx_sampled]
        end
    end
    return a_sampled
end

function sampling_BI(b_demands,w_sizes,W,n_samples)

    tt0 = time_ns()
    num_add = 0
    M = size(w_sizes)[1]
    A_constraint = zeros(Float64,M,n_samples)

    for ind_sample = 1 : n_samples
        new_a = sample_one_a_BI(w_sizes,W,b_demands)
        A_constraint[:,ind_sample] = new_a
    end

    return A_constraint

end
#test_CRLP = CRLP_primal(b_demands,w_sizes,W,40000)

# for m_width in [1000,2000,4000]

#     w_sizes = sort(sample(W/10:W/4,m_width,replace=false))
#     b_demands = sample(1:100,width_types)

#     test = dual_of_cutting_stock(b_demands,w_sizes,W,upper_limit_iter)
# end