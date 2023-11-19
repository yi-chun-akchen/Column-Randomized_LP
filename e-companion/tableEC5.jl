include("../main/choice_estimation_main.jl")

function rank_estimation_L1_warmstart(
    z, # z is a N by M binary matrix that shows assortments. N is number of options. M is number of historical assortments
    v, # v is a N by M matrix that shows choice distributions. v[:,1] is choice distr. of assortment 1
    iter, # maximum number of rankings to be learned
    warmstart_rankings # the rankings obtained by sampling
    )
	
	tt0 = time_ns()
    M = length(z[1,:]);N = length(z[:,1])
	time_obj_table = zeros(iter,2)

    ##### Solve the master problem by constraint generation

    master_problem = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(master_problem, -1 <= beta[1:M*N] <= 1)
    @variable(master_problem, mu)

    warmup_A_transpose = transpose(constraint_matrix_assortment_and_ranking_set(z,warmstart_rankings))
    for ind_order = 1 : length(warmstart_rankings)
        @constraint(master_problem, sum(warmup_A_transpose[ind_order,ind_p] * beta[ind_p] for ind_p in 1:M*N) + mu <= 0)
    end

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
        else

            break

        end

        solve(master_problem)
        dual_obj = getobjectivevalue(master_problem)
        tt1 = time_ns()
        println("iteration = ",n_rank_added," obj = ",dual_obj, " with time = ", (tt1-tt0)/1.0e9)
        time_obj_table[n_rank_added,:] = [dual_obj,(tt1-tt0)/1.0e9]

        if dual_obj <= 0
            break
        end

    end

    ###### Solve the primal problem
    total_num_rank = length(warmstart_rankings) + n_rank_added
    A = [transpose(warmup_A_transpose) reshape(transpose(A_transpose),(M*N,n_rank_added))]

    primal_problem = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(primal_problem,ep[1:M*N] >= 0)
    @variable(primal_problem,en[1:M*N] >= 0)
    @variable(primal_problem,lambda[1:total_num_rank] >= 0)

    @constraint(primal_problem,sum(lambda[ind_rank] for ind_rank in 1:total_num_rank) == 1)

    for ind_i = 1 : M*N
        @constraint(primal_problem, sum(A[ind_i,ind_rank] * lambda[ind_rank] for ind_rank in 1:total_num_rank)
                                        + ep[ind_i] - en[ind_i] == v[ind_i])
    end

    @objective(primal_problem, Min, sum(ep[ind_i] + en[ind_i] for ind_i in 1:M*N))

    solve(primal_problem)
	
	return getobjectivevalue(primal_problem), time_obj_table

end



K_list = Dict()
K_list[(7,50)] = [500,1000]
K_list[(9,50)] = [500,1000]
K_list[(9,100)] = [500,1000,1500]
K_list[(11,50)] = [500,1000]
K_list[(11,100)] = [500,1000,1500,2000]
K_list[(11,150)] = [500,1000,1500,2000,2500]
table_header = ["N","M","K","Obj","CR_time","CR-then-CG_time","CG-only_time"]

num_iter = 10000

n_instance = 1

Table_EC5 = zeros(Float64, n_instance*18,7)
row = 0

for instance = 1 : n_instance
    for (N,M) in [(7,50),(9,50),(9,100),(11,50),(11,100),(11,150)]
    	# note that the number of products in the table does not include the
    	# no-purchase option. Thus, the correponding value here is one larger.
		MNL_parameter = rand(N-1)
		asst_mat = generate_different_assortments(N,M)
		choice_vec = MNL_predict(asst_mat,MNL_parameter)

		CG_time = @elapsed CG = rank_estimation_L1(asst_mat,choice_vec,num_iter)

		for K in K_list[(N,M)]

			CR_time = @elapsed CR = random_rankings_estimation(asst_mat,choice_vec,K)
			CR_obj = CR[1]
			sampled_rankings = CR[end]
			CRCG_time = @elapsed CRCG =  rank_estimation_L1_warmstart(asst_mat,choice_vec,num_iter,sampled_rankings)

			row += 1
            Table_EC5[row,:] = [N-1,M,K,CR_obj,CR_time,CR_time+CRCG_time,CG_time]

            output_df = convert(DataFrame,Table_EC5[1:row,:])
            CSV.write(string("Table_EC5_instances.csv"),output_df,header = table_header)
        end
    end
end

