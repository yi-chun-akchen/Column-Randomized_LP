include("../main/cutting_stock_main.jl")

# The following code generates problem instances
# for the results in Table EC2 and saves the results
# as a CSV file

table_header = ["m","K","Delta_I","Delta_BI"]

W = 100000
n_instance = 1
K_list = Dict()
K_list[50] = [100,150,200,250,300]
K_list[100] = [200,250,300,350,400]

# setup 2
    #b_demands = sample(50:100,m_width)

# setup 1
Table_EC2_setup1 = zeros(Float64,n_instance*2*5,4)
row = 0
for i_instance = 1 : n_instance
    for m_width in [50,100]
        
        w_sizes =Int.(sort(sample(ceil(W/1000):floor(W/2),m_width,replace=false),rev=true))
        b_demands = sample(25:100,m_width)
    
        CG = dual_of_cutting_stock(b_demands,w_sizes,W,2000); opt_obj = CG[1]
        #A_I = sampling(b_demands,w_sizes,W,400)
        #A_BI = sampling_BI(b_demands,w_sizes,W,400)

        for K in K_list[m_width]
            global row += 1
            obj_I = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling(b_demands,w_sizes,W,K))[1]
            obj_BI = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling_BI(b_demands,w_sizes,W,K))[1]
            Table_EC2_setup1[row,:] = [m_width,K,(obj_I - opt_obj)/opt_obj, (obj_BI - opt_obj)/opt_obj]

            output_df = convert(DataFrame,Table_EC2_setup1[1:row,:])
            CSV.write(string("Table_EC2_setup1_instances.csv"),output_df,header = table_header)
        end
    end
end

# setup 2
Table_EC2_setup2 = zeros(Float64,n_instance*2*5,4)
row = 0
for i_instance = 1 : n_instance
    for m_width in [50,100]
        
        w_sizes =Int.(sort(sample(ceil(W/1000):floor(W/2),m_width,replace=false),rev=true))
        b_demands = sample(25:100,m_width)
    
        CG = dual_of_cutting_stock(b_demands,w_sizes,W,2000); opt_obj = CG[1]
        #A_I = sampling(b_demands,w_sizes,W,400)
        #A_BI = sampling_BI(b_demands,w_sizes,W,400)

        for K in K_list[m_width]
            global row += 1
            obj_I = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling(b_demands,w_sizes,W,K))[1]
            obj_BI = CRLP_with_sampled_columns(b_demands,w_sizes,W,sampling_BI(b_demands,w_sizes,W,K))[1]
            Table_EC2_setup2[row,:] = [m_width,K,(obj_I - opt_obj)/opt_obj, (obj_BI - opt_obj)/opt_obj]

            output_df = convert(DataFrame,Table_EC2_setup2[1:row,:])
            CSV.write(string("Table_EC2_setup2_instances.csv"),output_df,header = table_header)
        end
    end
end