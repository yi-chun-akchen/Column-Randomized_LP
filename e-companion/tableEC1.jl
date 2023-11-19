include("../main/cutting_stock_main.jl")


# The following code generates problem instances
# for the results in Table EC1 and saves the results
# as a CSV file

W = 100000
m_width = 5
n_instance = 1
table_header = ["K","Delta_U","Delta_I"]

# setup 1
Table_EC1_setup1 = zeros(Float64,4*n_instance,3)
row = 0
for i_instance = 1 : n_instance
    b_demands = sample(1:100,m_width)
    w_sizes =Int.(sort(sample(ceil(W/1000):floor(W/2),m_width,replace=false),rev=true))

    CG = dual_of_cutting_stock(b_demands,w_sizes,W,1000); opt_obj = CG[1]

    for K in [50,100,200,400]
        global row += 1
        obj_I = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling(b_demands,w_sizes,W,K))[1]
        obj_U = CRLP_with_sampled_columns(b_demands,w_sizes,W,reject_sampling_contraint_matrix(W,w_sizes,K))[1]
        Table_EC1_setup1[row,:] = [K,(obj_U - opt_obj)/opt_obj, (obj_I - opt_obj)/opt_obj]

        output_df = convert(DataFrame,Table_EC1_setup1[1:row,:])
        CSV.write(string("Table_EC1_setup1_instances.csv"),output_df,header = table_header)
    end
end

# setup 2
Table_EC1_setup2 = zeros(Float64,4*n_instance,3)
row = 0
for i_instance = 1 : n_instance
    b_demands = sample(1:100,m_width)
    w_sizes =Int.(sort(sample(ceil(W/10):floor(W/4),m_width,replace=false),rev=true))

    CG = dual_of_cutting_stock(b_demands,w_sizes,W,1000); opt_obj = CG[1]

    for K in [50,100,200,400]
        global row += 1
        obj_I = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling(b_demands,w_sizes,W,K))[1]
        obj_U = CRLP_with_sampled_columns(b_demands,w_sizes,W,reject_sampling_contraint_matrix(W,w_sizes,K))[1]
        Table_EC1_setup2[row,:] = [K,(obj_U - opt_obj)/opt_obj, (obj_I - opt_obj)/opt_obj]

        output_df = convert(DataFrame,Table_EC1_setup2[1:row,:])
        CSV.write(string("Table_EC1_setup2_instances.csv"),output_df,header = table_header)
    end
end
