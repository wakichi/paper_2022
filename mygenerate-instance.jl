# import Pkg
# Pkg.add("Concorde")


## パッケージ
using Concorde
using Random
using Plots
using Plots: plot, plot! # to supress VScode linter warning
using Printf
ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
#import Pkg
#Pkg.add("CPLEX")
#Pkg.build("CPLEX")
using JuMP
using CPLEX
#parameta
have_timewindow = true
is_round = true #time_windowを整数にするか否か

## データ初期化
n = 10 # ターゲット数
q_s = [0; 0]   # initial point   
q_f = [50; 0]  # final point 
q_s_dummy = [0; -1000]   # initial point  (to connect q_s and q_f)
q_f_dummy = [50; -1000]  # final point 
q_min = [0;0] # qの最小のx、y
q_max = [60;60]
p_o = [0;0] # startpoint
p_f = [50;0]# end point
v_v = 18 # km/h # default = 90 vehicle speed
rand_seed = 1024 + n
rand_seed = 1029
time_window_size = [6; 10] # time window の幅
time_window_min_max = [0;n] # uの開始時間の幅
Random.seed!(rand_seed)
function main()
    q,u = make_q()
    p1 = plot_q(q,u)
    exist_res = check_result(q,u)
    if exist_res
        save_q(q,u)
    end
end
function make_q()
    q = rand(2,n) * 48 .+ 1 # [1,49]x[1,49]
    u = [1 2 3 4 5;6 7 8 9 10] 
    u = rand(2,n)*(time_window_min_max[2] - time_window_min_max[1]) .+ time_window_min_max[1]
    u_size = rand(n)*(time_window_size[2] - time_window_size[1]) .+ time_window_size[1]
    u2 = u[1,:]+u_size
    u[2,:] = u2
    if is_round 
        u = round.(u)
    end
    return q, u
end
function plot_q(q,u)
    ## ここから Plots で確認する必要あり
    p1 = plot(xlims = (-5,55), ylims = (-5, 55), legend=:none)
    scatter!(p1, q[1,:], q[2,:])
    for i=1:n
        annotate!(p1, q[1,i]+1, q[2,i], text("q$(i)", :blue, 10))
        if have_timewindow
            annotate!(p1, q[1, i] + 1, q[2, i]+3, text(string("[", u[1,i], " ", u[2,i],"]"), :red, 5))
        end
    end
    scatter!(p1, [q_s[1] q_f[1]], [q_s[2] q_f[2]], color=:blue)
    annotate!(p1, q_s[1]+2, q_s[2], text("q_s", :blue, 10))
    annotate!(p1, q_f[1]+2, q_f[2], text("q_f", :blue, 10))
    display(p1)
    return p1
end

function check_result(q,u)
    ## vehicle_speedで解があるか調べる
    w, exist_answer = MISOCP_solverd_first_seq(q, u)
    println("exist_answer!!!  ",exist_answer)
    if exist_answer
        println("exist one of answer")
        w = matrix_to_vector(w)
        println("result:",w)
    else
        println("not exist answer")
    end
    return exist_answer
end

## 結果をファイルに記入
function save_q(q,u)
    filename1 = @sprintf("./tests/points-n-%04d-seed-%04d.csv", n, rand_seed)
    println("Writing file into $(filename1)")
    file1 = open(filename1, "w")
    # q = [ 0  0  1 46 48 46 50 50 20 30;
    #      40 50 49 50 48 46 48 35  5  0]
    println(size(q))
    println(size(u))
    for k=1:n-1
        print(file1, q[1,k], ",")
    end
    println(file1, q[1,n])
    for k=1:n-1
        print(file1, q[2,k], ",")
    end
    println(file1, q[2,n])

    for k=1:n-1
        print(file1, u[1,k], ",")
    end
    println(file1, u[1,n])
    for k=1:n-1
        print(file1, u[2,k], ",")
    end
    println(file1, u[2,n])

    close(file1)
    println("save Complete!!")
end
# targets = [[init[1]; Xpos][opt_tour] [init[2]; Ypos][opt_tour]],
# pos_init = findfirst(opt_tour.==1)[1]
# targets = [targets[1:pos_init-1, :]; targets[pos_init+1:end, :]]
# n = size(targets,1)

function MISOCP_solverd_first_seq(q, u)
    """
    droneなしVRPを解きます。
    とりあえず速度はトラック。
    間に合わないならドローンの速度にする
    """
    n = size(q, 2)
    model = Model()
    set_optimizer(model, CPLEX.Optimizer)
    # 決定変数
    @variable(model, w[1:n, 1:n],Bin)
    @variable(model, q_min[i]<=Q[i = 1:2, 1:n]<=q_max[i])
    @variable(model, T[i=2:n]>=0) # T=1、n＋1は下に分離
    # println("about T", T)#ok!
    @variable(model, T_first>=0) 
    @variable(model, T_last>=0) 
    if have_timewindow
        @variable(model, U[1:2,1:n]>=0)
    end
    # 目的関数
    @objective(model, Min,   (T_first+sum(T)+T_last))
    # 制約
    @constraint(model, c4, [v_v*T_first;Q[:, 1] - p_o[:,1]] in SecondOrderCone())
    @constraint(model, c5[j = 2:n], [v_v*T[j];Q[:,j-1] - Q[:,j]] in SecondOrderCone())
    @constraint(model, c6,[v_v*T_last;p_f - Q[:, n]] in SecondOrderCone())
    @constraint(model, c81[i = 1:n], (Q[1,i] == sum(w[i,j]*q[1,j] for j in 1:n)))
    @constraint(model, c82[i = 1:n], (Q[2,i] == sum(w[i,j]*q[2,j] for j in 1:n)))
    @constraint(model, c9[i = 1:n], 1 == sum(w[i,:])) # for i
    @constraint(model, c10[j = 1:n], 1 == sum(w[:, j])) # for j
    if have_timewindow
        @constraint(model, w1[i = 1:n], (T_first+sum(T[j+1] for j in 1:i-1))<=U[2,i])
        @constraint(model, w2[i = 1:n], U[1,i]<=(T_first+sum(T[j+1] for j in 1:i-1)))
        @constraint(model, w3[i = 1:n], (U[1,i] == sum(w[i,j]*u[1,j] for j in 1:n)))
        @constraint(model, w4[i = 1:n], (U[2,i] == sum(w[i,j]*u[2,j] for j in 1:n)))
    end
    set_time_limit_sec(model, 10800.0)
    elapsed_time = @elapsed optimize!(model)
    # @show termination_status(model)
    # @show objective_value(model)
    # @show elapsed_time # cputime
    exist_answer = true
    println("hello!!!!!!!!!!!!!!!")
    try
        objective_value(model)
    catch
        exist_answer = false
    end
    println("exist_answer!!?!:  ", exist_answer)
    return value.(w),exist_answer
end

function matrix_to_vector(matrix)
    """
    input: matrix, 正し、各列、各行に必ず1が一つ他は0
    output:vector, matrixの1となっている列の数字を入れたベクトル
    """
    n = size(matrix, 1)
    vector = zeros(Int, n)
    for i in 1:n
        for j in 1:n
            if matrix[i,j] == 1
                vector[i] = j
                break
            end
        end
    end
    return vector
end
main()