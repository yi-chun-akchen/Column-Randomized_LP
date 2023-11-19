include("cutting_stock_main.jl")

table_header = ["m","K","Optimality Gap"]
iter_upper = 20000
W = 100000
width_list = [1000,2000,4000]
column_size_list = [1,2,3,4,5,6,7,8] .* 10000
n_instance = 1

Figure_1 = zeros(Float64, n_instance*length(width_list)*length(column_size_list),3)
row = 0

for instance = 1 : n_instance
    for m_width in width_list
        w_sizes = sort(sample(W/10:W/4,m_width,replace=false))
        b_demands = sample(1:100,m_width)
        CG_outcome = dual_of_cutting_stock(b_demands,w_sizes,W,iter_upper)

        A_I = sampling(b_demands,w_sizes,W,80000)

        for K in column_size_list
            CRLP_obj = obj_I = CRLP_with_sampled_columns(b_demands,w_sizes,W,A_I[:,1:K])[1]
            global row += 1
            Figure_1[row,1] = m_width
            Figure_1[row,2] = K
            Figure_1[row,3] = (CRLP_obj - CG_outcome[1])/CG_outcome[1]

            output_df = convert(DataFrame,Figure_1[1:row,:])
            CSV.write(string("Figure_1_instances.csv"),output_df,header = table_header)
        end
    end
end