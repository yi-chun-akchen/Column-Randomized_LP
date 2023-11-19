using JuMP, Gurobi, Distributions, StatsBase, Distances, DataFrames, CSV

# ------------------------------------------------------------------------ #
# Column Generation Approach for the Nonparametric Choice Model Estimation #
# ------------------------------------------------------------------------ #

function solve_subproblem_for_ranking(u, z)
    # This function solves the subproblem of the ranking-based model estimation.
    # The first inuput u is the dual variable of constraint (14b).The second input
    # z is a N by M binary matrix, which represents a collection of assortments.
    # Please refer to van Ryzin and Vulcano (2015) and Misic (2016) for the formation
    # of the subproblem.

    # M is the number of assortments
    M = size(z)[2]
    # N is the number of products
    N = size(z)[1]

    sub_rank = Model(solver = GurobiSolver((OutputFlag = 0)))
    @variable(sub_rank, x[1:N,1:N], Bin)
    @variable(sub_rank, a[1:N,1:M], Bin)


    for ind_s = 1 : M
        @constraint(sub_rank, sum(a[ind_p,ind_s] for ind_p in 1:N) == 1)
        for ind_p1 = 1 : N
            @constraint(sub_rank,a[ind_p1,ind_s] <= z[ind_p1,ind_s])
            for ind_p2 = 1 : N
                if ind_p1 != ind_p2 && z[ind_p1,ind_s] == 1 && z[ind_p2,ind_s] == 1
    	            @constraint(sub_rank, a[ind_p1,ind_s] <= x[ind_p1,ind_p2])
                end
           end
        end
    end

    for ind_p1 = 1 : N
        for ind_p2 = 1 : N
            if ind_p1 != ind_p2
                @constraint(sub_rank, x[ind_p1,ind_p2] + x[ind_p2,ind_p1] == 1)
            end
        end
    end


    for ind_p1 = 1 : N
        for ind_p2 = 1 : N
            for ind_p3 = 1 : N
                if ind_p1 != ind_p2 && ind_p2 != ind_p3 && ind_p1 != ind_p3
                    @constraint(sub_rank, x[ind_p1,ind_p2] + x[ind_p2,ind_p3] <= 1 + x[ind_p1,ind_p3])
                end
            end
        end
    end

    @objective(sub_rank, Max,  sum(sum(u[ind_p,ind_s] * a[ind_p,ind_s] for ind_p in 1:N) for ind_s in 1:M))
    solve(sub_rank)


    obj_value = getobjectivevalue(sub_rank)

    rank_rep = round.(Int64,getvalue(x))
    rank_sum = zeros(Int64,N)
    for ind_p= 1 : N
        rank_sum[ind_p] = sum(rank_rep[ind_p,:])
    end

    rank = sortperm(rank_sum,rev = true)

    new_A_row = reshape(getvalue(a),(1,N*M))

    return [obj_value, new_A_row, rank]

end

function purchase_decision_and_rank_mapping(assorment,rank)
    # this function returns the choice of a ranking under a given assortment
    for ind_r = 1 : length(rank)
        if assorment[rank[ind_r]] == 1
            return rank[ind_r]
        end
    end
end



function rankings_select_o_in_S_m(z,R)
    # z is the collection of assortments, a N-by-M binary matrix
    # R is the collection of rankings
    M = length(z[1,:]); N = length(z[:,1]);

    map_om_to_ranks = Dict()


    for ind_m = 1 : M
        for ind_o = 1 : N
            map_om_to_ranks[(ind_o,ind_m)] = []
        end
    end

    for ind_m = 1 : M
        for ind_r = 1 : length(R)
            o_selected = purchase_decision_and_rank_mapping(z[:,ind_m],R[ind_r])
            map_om_to_ranks[(o_selected,ind_m)] = [map_om_to_ranks[(o_selected,ind_m)]; ind_r]
        end
    end

    return map_om_to_ranks
end


function predict_distri_by_ranks(R,weights,z)

    # R is a collection of rankings. weights is a distribution over R
    # (R,weights) together is a ranking-based model.
    # This function returns the resulting choice probability on assortments
    # z[:,m], for m = 1,...,M, by the ranking-based model.

    M = length(z[1,:]); N = length(z[:,1]);
    pred_prob = zeros(Float64,N,M)

    for ind_s = 1 : M
        for ind_r = 1 : length(R)
            p = purchase_decision_and_rank_mapping(z[:,ind_s],R[ind_r])
            pred_prob[p,ind_s] += weights[ind_r]
        end
    end

    return pred_prob
end


function rank_estimation_L1(z,v,iter)
    # z is a N by M binary matrix that represents a collection of assortments.
    # Here each z[:,m], m = 1,...,M, represents an assortment. z[i,m] shows
    # whether product i is in assortment m, for i=1,..,N-1.
    # z[N,m] is always 1, as it represents the no-purchase option.
    # v is the corresponding choice probability; see Section 6.
    # iter is the maximum number of rankings to be added; set to a high number.
    # It is used to ensure the termination of the iterations.

    tt0 = time_ns()

    M = length(z[1,:]); N = length(z[:,1])
    traj = zeros(Float64,iter,2)

    ##### Solve the master problem by constraint generation

    master_problem = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(master_problem, -1 <= beta[1:M*N] <= 1)
    @variable(master_problem, mu)
    @objective(master_problem, Max, sum(v[ind_i] * beta[ind_i] for ind_i in 1 : M*N) + mu)
    solve(master_problem)

    rank_added = Array{Any}(iter)
    n_rank_added = 0
    A_transpose = []

    dual_obj = 0

    for iteration = 1 : iter

        beta_copy = getvalue(beta)
        mu_copy = getvalue(mu)

        alpha = reshape(beta_copy,(N,M))
        subproblem_result = solve_subproblem_for_ranking(alpha, z)
        #println(" new rank = ",subproblem_result[3], subproblem_result[2])

        if subproblem_result[1] + mu_copy > 0 #- eps(Float32)

            n_rank_added = n_rank_added + 1

            new_row = subproblem_result[2]
            @constraint(master_problem,
                        sum(new_row[ind_i] * beta[ind_i] for ind_i in 1 : M*N) + mu <= 0)

            rank_added[iteration] = subproblem_result[3]
            A_transpose = append!(A_transpose,new_row)
            #print("Add rank #",n_rank_added)
        else

            break

        end

        solve(master_problem)
        dual_obj = getobjectivevalue(master_problem)
        tt1 = time_ns()
        traj[n_rank_added,:] = [(tt1-tt0)/1.0e9,dual_obj]
        println("iteration = ",n_rank_added," obj = ",dual_obj, " with time = ", (tt1-tt0)/1.0e9)

        if dual_obj <= 0
            break
        end

    end

    ###### Solve the primal problem

    A = reshape(transpose(A_transpose),(M*N,n_rank_added))
    primal_problem = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(primal_problem,ep[1:M*N] >= 0)
    @variable(primal_problem,en[1:M*N] >= 0)
    @variable(primal_problem,lambda[1:n_rank_added] >= 0)

    @constraint(primal_problem,sum(lambda[ind_rank] for ind_rank in 1:n_rank_added) == 1)

    for ind_i = 1 : M*N
        @constraint(primal_problem, sum(A[ind_i,ind_rank] * lambda[ind_rank] for ind_rank in 1:n_rank_added)
                                        + ep[ind_i] - en[ind_i] == v[ind_i])
    end

    @objective(primal_problem, Min, sum(ep[ind_i] + en[ind_i] for ind_i in 1:M*N))

    solve(primal_problem)

    tt1 = time_ns()
    traj[n_rank_added,:] = [(tt1-tt0)/1.0e9,getobjectivevalue(primal_problem)]

    return getobjectivevalue(primal_problem), getValue(lambda), rank_added, getValue(ep), getValue(en), A, traj[1:n_rank_added,:]

end

# ------------------------------------------------------------------------- #
# Column Randomization Method for the Nonparametric Choice Model Estimation #
# ------------------------------------------------------------------------- #


function sample_rankings(N,N_R)
    # N is number of options (including no-purchase option)
    # N_R is number of rankings to sample
    R = Array{Any}(N_R)
    for r = 1 : N_R
        R[r] = sample(1:N, N, replace = false)
    end
    return R
end

function constraint_matrix_assortment_and_ranking_set(z,R)
    # We construct the constraint matrix in Equation 14(b) given a set of rankings R.
    # Here z is a collection of assortments, where z[:,m] is an assortment.
    M = size(z)[2]; N = size(z)[1]

    constraint_matrix = zeros(Float64,M*N,length(R))
    for ind_r = 1 : length(R)
        for ind_s = 1 : M
            assort = z[:,ind_s]
            ranking = R[ind_r]
            choice = purchase_decision_and_rank_mapping(assort,ranking)
            constraint_matrix[N*(ind_s-1)+choice,ind_r] = 1
        end
    end
    return constraint_matrix
end


function random_rankings_estimation(z,v,N_R)

    # z and v follow the definition in function rank_estimation_L1(z,v,iter)
    # N_R is the number of rankings (i.e., columns) we are going to sample
    (N,M) = size(z)
    rankings_set = sample_rankings(N,N_R)
    sampled_A_matrix = constraint_matrix_assortment_and_ranking_set(z,rankings_set)

    RTS = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(RTS,ep[1:M*N] >= 0)
    @variable(RTS,en[1:M*N] >= 0)
    @variable(RTS,lambda[1:N_R] >= 0)

    @constraint(RTS,sum(lambda[ind_rank] for ind_rank in 1:N_R) == 1)

    for ind_i = 1 : M*N
        @constraint(RTS, sum(sampled_A_matrix[ind_i,ind_rank] * lambda[ind_rank] for ind_rank in 1:N_R)
                                        + ep[ind_i] - en[ind_i] == v[ind_i])
    end

    @objective(RTS, Min, sum(ep[ind_i] + en[ind_i] for ind_i in 1:M*N))

    solve(RTS)

    lambda_copy = getValue(lambda)
    n_nz = 0
    A_nonzero = zeros(Float64,M*N,N_R)

    for n_sam = 1 : N_R
        if lambda_copy[n_sam] > 0
            n_nz += 1
            A_nonzero[:,n_nz] = sampled_A_matrix[:,n_sam]
        end
    end
    rankings_nz = Array{Any}(n_nz)
    n_nzz = 0
    for n_sam = 1 : N_R
        if lambda_copy[n_sam] > 0
            n_nzz += 1
            rankings_nz[n_nzz] = rankings_set[n_sam]
        end
    end
    
    return (getobjectivevalue(RTS),A_nonzero[:,1:n_nz],rankings_nz)
    # We return the objective value, the columns with nonzero value in the solution
    # , and the corresponding columns.
end


function generate_different_assortments(N,M)

    assort_index = sample(1:2^(N-1),M,replace=false)
    z = zeros(Int64,N,M); z[N,:] = 1

    for ind_ass = 1 : M
        assort_binary = split(string(assort_index[ind_ass]-1,base=2),"")
        len_str = length(assort_binary)
        if assort_index[ind_ass] != 1
            for ind_digit = 1:length(assort_binary)
                location = N - 1 - len_str + ind_digit
                z[location,ind_ass] = parse(assort_binary[ind_digit])
            end
        end
    end
    return z
end


function MNL_predict(z_out,u)

    M = size(z_out)[2]; N = size(z_out)[1]; v = zeros(Float64,size(z_out))

    for ind_s = 1 : M
        v[N,ind_s] = 1.0
        for ind_p = 1 : N-1
            v[ind_p,ind_s] = exp(u[ind_p]) * z_out[ind_p,ind_s]
        end
        v[:,ind_s] = v[:,ind_s] / sum(v[:,ind_s])
    end
    return v
end