include("choice_estimation_main.jl")

table_header = ["N","M","K","Obj","Runtime","CG_Runtime"]

K_list = Dict()
K_list[(7,50)] = [500,1000]
K_list[(9,50)] = [500,1000]
K_list[(9,100)] = [500,1000,1500]
K_list[(11,50)] = [500,1000]
K_list[(11,100)] = [500,1000,1500,2000]
K_list[(11,150)] = [500,1000,1500,2000,2500]

num_iter = 10000


n_instance = 1

Table_2 = zeros(Float64, n_instance*18,6)
row = 0

for instance = 1 : n_instance
    for (N,M) in [(7,50),(9,50),(9,100),(11,50),(11,100),(11,150)]
    	# note that the number of products in the table does not include the
    	# no-purchase option. Thus, the correponding value here is one larger.
		MNL_parameter = rand(N-1)
		asst_mat = generate_different_assortments(N,M)
		choice_vec = MNL_predict(asst_mat,MNL_parameter)

		CG_time = @elapsed CG = rank_estimation_L1(asst_mat,choice_vec,num_iter)
		CG_traj = CG[end]

		for K in K_list[(N,M)]

			CR_time = @elapsed CR = random_rankings_estimation(asst_mat,choice_vec,K)
			CR_obj = CR[1]
			CR_time_in_CG = CG_traj[findfirst(x -> x <= CR_obj, CG_traj[:,2]),1]

			#println((N-1,M,K,CR_obj,CR_time,CR_time_in_CG))
			
			row += 1
            Table_2[row,:] = [N-1,M,K,CR_obj,CR_time,CR_time_in_CG]

            output_df = convert(DataFrame,Table_2[1:row,:])
            CSV.write(string("Table_2_instances.csv"),output_df,header = table_header)
        end
    end
end