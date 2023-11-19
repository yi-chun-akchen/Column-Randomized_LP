using DataFrames, CSV

N = 100
delta = 0.1

Fig_EC2_data = zeros(Float64,40,3)
table_header = ["K","IID Bound","Permutation Bound"]
row = 0

for K = 1 : 3: 99

    row += 1

    H = ((N-0.5)/(N-K)) * (1 - (1/(2 * maximum([K,N-K]))))

    bound_1 = (1/(K^0.5)) * ( 1+( 2*log(1/delta) )^0.5)
    bound_2 = (1/(K^0.5)) * (( 1-(K-1)/(N-1) )^0.5 + ( (2/H)*log(1/delta) )^0.5)
	
	Fig_EC2_data[row,:] = [K,bound_1,bound_2]

    output_df = convert(DataFrame,Fig_EC2_data[1:row,:])
    CSV.write(string("Fig_EC2.csv"),output_df,header = table_header)
end
