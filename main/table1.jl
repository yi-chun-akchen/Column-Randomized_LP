include("cutting_stock_main.jl")

table_header = ["m","K","Optimality Gap","Runtime","CG Runtime","CG Runtime total"]
iter_upper = 20000
W = 100000
width_list = [1000,2000,4000]
column_size_list = [2,4,6,8] .* 10000
n_trail = 1

Table_1 = zeros(Float64, n_trail*length(width_list)*length(column_size_list),6)
row = 0

for trial = 1 : n_trail
    for m_width in width_list
        w_sizes = sort(sample(W/10:W/4,m_width,replace=false))
        b_demands = sample(1:100,m_width)
        CG_outcome = dual_of_cutting_stock(b_demands,w_sizes,W,iter_upper)
        CG_traj = CG_outcome[3]

        for K in column_size_list
            time_CRLP = @elapsed CRLP_obj = CRLP_primal(b_demands,w_sizes,W,K)[1]
            global row += 1
            Table_1[row,1] = m_width
            Table_1[row,2] = K
            Table_1[row,3] = (CRLP_obj - CG_outcome[1])/CG_outcome[1]
            Table_1[row,4] = time_CRLP

            # the runtime for CG to reach the same optimality level as column randomization.
            Table_1[row,5] = CG_traj[findfirst(x -> x <= CRLP_obj, CG_traj[:,2]),1]
            # the total runtime for CG to find the optimal solution
            Table_1[row,6] = CG_traj[end,1]

            output_df = convert(DataFrame,Table_1[1:row,:])
            CSV.write(string("Table_1_instances.csv"),output_df,header = table_header)
        end
    end
end