using JuMP, Ipopt, NLopt, Distributions

include("../main/choice_estimation_main.jl")



function solve_MNL(z, count1)
	# This function estimates MNL model by maximum likelihood estimation
	# z is the binary matrix that represents a collection of assortments.
	# count1 is the integer matrix, where count1[i,m] is the number of times
	# that product i is chosen given that assortment m is offered.
    
    M = size(z)[2]; N = size(z)[1]

    cc = sum(count1,1)
    MNL = Model(solver=NLoptSolver(algorithm=:LD_MMA))
    @variable(MNL, -20 <= u[1:N-1] <= 20)

    @NLobjective(MNL, Max,
    sum(sum(count1[ind_p,ind_s] * u[ind_p] for ind_p in 1:N-1) for ind_s in 1:M) -
    sum(cc[ind_s]*log(1+sum(z[ind_p,ind_s]* exp(u[ind_p]) for ind_p in 1:N-1)) for ind_s in 1:M)
    )
    solve(MNL)

    uu = getvalue(u)

    return [uu[:];0]

end



function sampling_column_from_MNL(z,v,num_col)

	# This function follows Algorithm 11 in Section G.1

    pseudo_count = floor.(Int64,v * 1000)
    estimated_MNL = solve_MNL(z, pseudo_count)

    G = Gumbel(0,1)
    sampled_columns = Vector(undef,num_col)
    for ind_col = 1 : num_col
        order_table = zeros(size(z)[1],2)
        order_table[:,2] = 1:size(z)[1]
        order_table[:,1] = rand(G,size(z)[1]) + estimated_MNL
        sampled_rank = Int.(sortrows(order_table)[end:-1:1,2])
        sampled_columns[ind_col] = sampled_rank
    end
    return sampled_columns
end



function random_rankings_estimation_via_MNL(z,v,N_R)

    (N,M) = size(z)

    rankings_set = sampling_column_from_MNL(z,v,N_R)
    A_matrix = constraint_matrix_assortment_and_ranking_set(z,rankings_set)


    RTS = Model(solver = GurobiSolver((OutputFlag = 0)))

    @variable(RTS,ep[1:M*N] >= 0)
    @variable(RTS,en[1:M*N] >= 0)
    @variable(RTS,lambda[1:N_R] >= 0)

    @constraint(RTS,sum(lambda[ind_rank] for ind_rank in 1:N_R) == 1)

    for ind_i = 1 : M*N
        @constraint(RTS, sum(A_matrix[ind_i,ind_rank] * lambda[ind_rank] for ind_rank in 1:N_R)
                                        + ep[ind_i] - en[ind_i] == v[ind_i])
    end

    @objective(RTS, Min, sum(ep[ind_i] + en[ind_i] for ind_i in 1:M*N))

    solve(RTS)

    return getobjectivevalue(RTS)

end



n_instance = 1

K_list = Dict()
K_list[(7,50)] = [500,1000]
K_list[(9,50)] = [500,1000]
K_list[(9,100)] = [500,1000,1500]
K_list[(11,50)] = [500,1000]
K_list[(11,100)] = [500,1000,1500,2000]

table_header = ["N","M","K","Obj_uni","Obj_MNL"]

# ------- #
# Setup 1 #
# ------- #
Table_EC4_setup1 = zeros(Float64, n_instance*13,5)

row = 0
for ind_instance = 1 : n_instance
	for (N,M) in [(7,50),(9,50),(9,100),(11,50),(11,100)]
    	# note that the number of products in the table does not include the
    	# no-purchase option. Thus, the correponding value here is one larger.

    	MNL_parameter = rand(N-1)
		asst_mat = generate_different_assortments(N,M)
		choice_vec = MNL_predict(asst_mat,MNL_parameter)

    	for K in K_list[(N,M)]

    		global row += 1

        	# uniform randomization
        	rho_u_obj = random_rankings_estimation(asst_mat,choice_vec,K)[1]
		
			# MNL randomization
        	rho_M_obj = random_rankings_estimation_via_MNL(asst_mat,choice_vec,K)

        	# Write file
        	Table_EC4_setup1[row,:] = [N,M,K,rho_u_obj,rho_M_obj]
        	output_df = convert(DataFrame,Table_EC4_setup1[1:row,:])
            CSV.write(string("Table_EC4_setup1_instances.csv"),output_df,header = table_header)
        end
    end
end



# ------- #
# Setup 2 #
# ------- #
Table_EC4_setup2 = zeros(Float64, n_instance*13,5)

row = 0
for ind_instance = 1 : num_instance
	for (N,M) in [(7,50),(9,50),(9,100),(11,50),(11,100)]
    	# note that the number of products in the table does not include the
    	# no-purchase option. Thus, the correponding value here is one larger.

    	MNL_parameter = rand(N-1) * 20
		asst_mat = generate_different_assortments(N,M)
		choice_vec = MNL_predict(asst_mat,MNL_parameter)

    	for K in K_list[(N,M)]

    		global row += 1

        	# uniform randomization
        	rho_u_obj = random_rankings_estimation(asst_mat,choice_vec,K)[1]
		
			# MNL randomization
        	rho_M_obj = random_rankings_estimation_via_MNL(asst_mat,choice_vec,K)

        	# Write file
        	Table_EC4_setup2[row,:] = [N,M,K,rho_u_obj,rho_M_obj]
        	output_df = convert(DataFrame,Table_EC4_setup2[1:row,:])
            CSV.write(string("Table_EC4_setup2_instances.csv"),output_df,header = table_header)
        end
    end
end
