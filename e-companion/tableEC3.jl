include("../main/cutting_stock_main.jl")

# The following code generates problem instances for the results in Table EC3
# and saves the results as a CSV file

# In the following two functions, we first solve a CRLP and obtain the columns
# with nonzero value in the solution. We then serve these columns as the initial
# solution for the column generation.

function CRLP_primal_nonzero_columns(b_demands,w_sizes,W,n_samples)

    M = size(w_sizes)[1]

    A_constraint = sampling(b_demands,w_sizes,W,n_samples)
    RS_problem = Model(solver = GurobiSolver(OutputFlag = 0))

    @variable(RS_problem,x_RS[1:n_samples] >= 0)
    @constraint(RS_problem, A_constraint * x_RS .>= b_demands)
    @objective(RS_problem, Min,sum(x_RS[ind_a] for ind_a in 1:n_samples))
    solve(RS_problem)
    x_copy = getValue(x_RS)
    n_nz = 0
    A_nonzero = zeros(Float64,M,n_samples)
    for n_sam = 1 : n_samples
        if x_copy[n_sam] > 0
            n_nz += 1
            A_nonzero[:,n_nz] = A_constraint[:,n_sam]
        end
    end
    return (getobjectivevalue(RS_problem),A_nonzero[:,1:n_nz])
end

function dual_of_cutting_stock_with_CRLP(b_demands,w_sizes,W,initial_A_transpose,iter)

    t0 = time_ns()
    M = size(b_demands)[1]

    dual_problem = Model(solver = GurobiSolver(OutputFlag = 0))

    @variable(dual_problem,p[1:M] >= 0)

    # constraints from sampling
    @constraint(dual_problem, initial_A_transpose * p .<= 1)

    @objective(dual_problem, Max, sum(b_demands[ind_w] * p[ind_w] for ind_w in 1:M))

    solve(dual_problem)
    dual_obj = getobjectivevalue(dual_problem)
    previous_obj = Inf

    for it_time = 1 : iter-1

        p_value = getvalue(p)
        subproblem_outcome = subproblem_cutting_stock(p_value,w_sizes,W)

        if subproblem_outcome[1] > 1.0
            a_new = subproblem_outcome[2]
            @constraint(dual_problem, sum(a_new[ind_w] * p[ind_w] for ind_w in 1:M) <= 1.0)
        else
            break
        end

        solve(dual_problem)
        dual_obj = getobjectivevalue(dual_problem)
        t1 = time_ns()
    end

    return dual_obj
end


table_header = ["m","K","Delta_I","T_CR","T_CR_then_CG","T_CG_only"]

W = 100000
n_instance = 10
Table_EC3 = zeros(Float64,500,6)
num_iter = 10000
row = 0
for i_instance = 1 : n_instance
    for m_width in [250,500,750,1000,1500]

        global row += 1

        w_sizes = sort(sample(W/10:W/4,m_width,replace=false))
        b_demands = sample(1:100,m_width)

        time_CG = @elapsed opt_obj = dual_of_cutting_stock(b_demands,w_sizes,W,num_iter)[1]
        time_CRLP = @elapsed CRLP_outcome = CRLP_primal_nonzero_columns(b_demands,w_sizes,W,m_width*10)
        time_CG_after_CR = @elapsed CG_after_CR = dual_of_cutting_stock_with_CRLP(b_demands,w_sizes,W,transpose(CRLP_outcome[2]),num_iter)

        Table_EC3[row,1] = m_width
        Table_EC3[row,2] = m_width*10
        Table_EC3[row,3] = (CRLP_outcome[1] - opt_obj)/opt_obj
        Table_EC3[row,4] = time_CRLP
        Table_EC3[row,5] = time_CG_after_CR
        Table_EC3[row,6] = time_CG

        output_df = convert(DataFrame,Table_EC3[1:row,:])
        CSV.write(string("Table_EC3_instances.csv"),output_df)
    end
end