using Plots
using LinearAlgebra
using Plots: plot, plot! # to supress VScode linter warningd
using CSV
using DataFrames
using Random
using StatsBase # 重複なく2つ抜き出すのに使用
using BenchmarkTools
using Dates
# NOTE: paper_2022上で動かすことを想定！

ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
#import Pkg
#Pkg.add("CPLEX")
#Pkg.build("CPLEX")
using JuMP
using CPLEX

# パラメータの設定(global)
is_includecsv = true
have_timewindow  = true
is_plot = true
is_save_pic = false
n_iteration =50
is_movie_mode = false

q_min = [0;0] # qの最小のx、y
q_max = [60;60]
v_c = 18 # km/h # default = 18 carrior speed
v_v = 60 # km/h # default = 90 vehicle speed
a = 21/60 # vehicle operation range
p_o = [0;0] # startpoint
p_f = [50;0]# end point
infeasible_score = 500 # infeasibleの場合に返す値
max_runtime = 600.0

test_file = "/Users/wakitakouhei/Lab/paper_2022/src/resources/p0040/points-n-0040-seed-1064.csv"
save_path = "/Users/wakitakouhei/Lab/paper_2022/src/Experiments.csv"
function main()
    # 入力
    q, u = input()
    # ヒューリスティックの実行部分
    runtime = @elapsed begin
    score, seq, p_to, p_land, status, cnt = improved_SA(q, u)
    end
    println("count_all_iter:", cnt)
    println("runtime:",runtime)
    println("score:",score)
    is_infeasible = (score == infeasible_score)
    #status　は以下の二つ
    # NORMAL_TERMINATE
    # TIME_LIMIT
    write_to_csv(save_path, n, runtime, score, n_iteration, true, is_infeasible, status, test_file)
    # 出力
    # 7.00729450915383
    if !is_infeasible
        # solutionが少なくとも一つ
        p2 = output( p_to, p_land, seq, q)
        if is_save_pic
            save_figure("/Users/wakitakouhei/Lab/paper_2022/temp/pictures/", "aa", p2)
        end
    end
end
function input()
    if is_includecsv
        println("loading the data in $(test_file)")
        df = CSV.read(test_file, DataFrame; header=0)    # test data by "generate-instance.jl"
        q = [df[i, j] for i = 1:2, j = 1:size(df, 2)]
        # TODO: 3,4行目が存在しない場合、printした上でhave_timewindowをfalseにする。
        if have_timewindow
            u = [df[i, j] for i = 3:4, j = 1:size(df, 2)]
        else
            u = [;]
        end
    else
        q = [45 10 30 
        45 40 20 ]
        u = [0 0 200 
           200 200 300 ]
    end
    global n = size(q, 2)
    if is_plot
        global p1 = plot(xlims=(-5, 55), ylims=(-5, 55), legend=:none, aspect_ratio=:equal)
        scatter!(p1, q[1, :], q[2, :])
        for i = 1:n
            if have_timewindow
                annotate!(p1, q[1, i] + 1, q[2, i]+3, text(string("[", u[1,i], " ", u[2,i],"]"), :red, 5))
            end
            annotate!(p1, q[1, i] + 1, q[2, i]+1, text("q$(i)", :blue, 10))
        end
        scatter!(p1, [p_o[1] p_f[1]], [p_o[2] p_f[2]], color=:blue)
        annotate!(p1,p_o[1] + 2, p_o[2], text("p_s", :blue, 10))
        annotate!(p1, p_f[1] + 2, p_f[2], text("p_f", :blue, 10))
        display(p1)
    end
    return (q,u)
end

function climing(q,u)
    # 初期解の生成(どういう形？)(1 0 0 0;0 0 1 0)的な形
    seq = make_first_seq(q,u)

    score,p_to,p_land = calc_score(seq, q, u)
    for i = 1:n_iteration
        # 近傍をとる。
        n_seq = make_new_seq(seq)
        # スコアの算出(cplexに投げる)
        n_score,n_p_to, n_p_land = calc_score(n_seq, q, u)
        # 順列を更新するか判断する。
        println("i:", i)
        output_easy(n_seq, score)
        if score>=n_score # 最小化問題なので
            if (is_movie_mode && score!= n_score)
                output(p_to, p_land, n_seq, q)
            end
            score = n_score
            seq = n_seq
            p_land = n_p_land
            p_to = n_p_to
        end
    end
    # 結果を可視化するならここで。
    return (score, seq, p_to, p_land)
end

function annealing(q,u)
    seq = make_first_seq(q,u)
    first_temp = 1
    end_temp = 0.001
    output_easy(seq, 0)
    score,p_to,p_land = calc_score(seq, q, u)
    first_score = score
    best_score, best_p_to, best_p_land,best_seq= score, p_to, p_land, seq
    for i = 1:n_iteration
        temp = first_temp+(end_temp - first_temp)*(i/n_iteration)
        # 近傍をとる。
        n_seq = make_new_seq(seq)
        # スコアの算出(cplexに投げる)
        n_score,n_p_to, n_p_land = calc_score(n_seq, q, u)
        # 順列を更新するか判断する。
        output_easy(n_seq, score)
        prob = exp((score - n_score)/temp)
        # best case 
        if n_score<best_score
            best_score = n_score
            best_p_to = n_p_to
            best_p_land = n_p_land
            best_seq = n_seq
        end
        if prob> rand()
            if (is_movie_mode && score!= n_score)
                output(p_to, p_land, n_seq, q)
            end
            # 更新
            score = n_score
            seq = n_seq
            p_land = n_p_land
            p_to = n_p_to
        end
    end
    println("firstscore:", first_score)
    output_easy(best_seq, best_score)
    return (best_score, best_seq, best_p_to, best_p_land)
end

function improved_SA(q,u)
    """
    終了条件: すくなくとも一つ実行可能解を見つけた上で、指定回数回解が改善されなかったら終了
    もしくは開始から600秒が過ぎた段階で強制終了
    """
    # 初期設定
    have_feasible_solution = false
    update_solution = false # 一回のiterationの間で解が更新されたかどうか
    status = "NORMAL_TERMINATE"
    cnt_iter = 0
    first_temp = 10
    end_temp = 0.001
    start_time = now()

    seq = make_first_seq(q,u)
    score,p_to,p_land = calc_score(seq, q, u)
    best_score, best_p_to, best_p_land,best_seq= score, p_to, p_land, seq

    while(!(have_feasible_solution && !update_solution))
        println("uifhaiughi;oahg;a")
        update_solution = false
        cnt_iter+=1
        #TODO: 時間超過を判断
        runtime = parse(Float64,string(now()-start_time)[1:end-12]) #millisecond
        if (runtime>max_runtime*1000)
            status = "TIME_LIMIT"
            break
        end
        if !have_feasible_solution
            # 解が存在する場合はそのまま使用。それ以外は新しく生成
            seq = make_first_seq(q,u)
            score,p_to,p_land = calc_score(seq, q, u)
            best_score, best_p_to, best_p_land,best_seq= score, p_to, p_land, seq
        end
        #TODO: iterかいSA
        for i = 1:n_iteration
            temp = first_temp+(end_temp - first_temp)*(i/n_iteration)
            # 近傍をとる。
            n_seq = make_new_seq(seq)
            # スコアの算出(cplexに投げる)
            n_score,n_p_to, n_p_land = calc_score(n_seq, q, u)
            # 順列を更新するか判断する。
            output_easy(n_seq, score)
            prob = exp((score - n_score)/temp)
            # best case 
            if n_score<best_score
                best_score = n_score
                best_p_to = n_p_to
                best_p_land = n_p_land
                best_seq = n_seq
                update_solution = true
            end
            if prob> rand()
                # 更新
                score = n_score
                seq = n_seq
                p_land = n_p_land
                p_to = n_p_to
            end
        end
        #TODO: 解が更新されたかや最適解があるかどうか検知
        if (best_score != infeasible_score) 
            have_feasible_solution = true
        end
    end
    return (best_score, best_seq, best_p_to, best_p_land, status, cnt_iter)
end

function make_first_seq(q,u)
    # seq = MISOCP_solverd_first_seq(q, u)
    seq = shuffle_first_seq(q,u)
    # seq = simple_first_seq(q,u)
    return seq
end

function simple_first_seq(q, u)
    n = size(q, 2) # the number of target points 
    seq = Matrix{Int32}(1I, (n, n)) # 単位行列
    return seq
end

function shuffle_first_seq(q, u)
    n = size(q,2)
    vecter_to_matrix
    seq = vecter_to_matrix(shuffle(1:n))
    return seq
end

function Integer_solverd_first_seq(q, u)
    # TODO: 実装
    distance_matrix = calc_distance(q)
    distance_matrix_first = calc_distance_with(q, p_o)# TODO:q1との距離のベクトル
    distance_matrix_last = calc_distance_with(q, p_f)# qnとの距離のベクトル
    M = 10^10

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
    # TODO:k,l, iの調整をする
    @constraint(model, c1[k = 1:n], distance_matrix_first[k] <= v_c*T_first + M*(1-w[1,k]))
    @constraint(model, c2[k = 1:n, l = 1:n, i=2:n ], distance_matrix[k,l] <= v_c*T[i]+M*(1-w[i,k])+M*(1-w[i-1,l]))
    @constraint(model, c3[k = 1:n], distance_matrix_last[k] <= v_c*T_last + M*(1-w[n, k]))
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
    @show termination_status(model)
    @show objective_value(model)
    @show elapsed_time # cputime
    return value.(w)
end

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
    @constraint(model, c4, [v_c*T_first;Q[:, 1] - p_o[:,1]] in SecondOrderCone())
    @constraint(model, c5[j = 2:n], [v_c*T[j];Q[:,j-1] - Q[:,j]] in SecondOrderCone())
    @constraint(model, c6,[v_c*T_last;p_f - Q[:, n]] in SecondOrderCone())
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
    @show termination_status(model)
    @show objective_value(model)
    @show elapsed_time # cputime
    return value.(w)


end

function calc_score(w, q, u)
    # cplexを用いてsequenceからスコアを計算する。
    # スコアは単純に目的関数値を用いる。
    n = size(q, 2) # the number of target points 

    model = Model()
    # set_optimizer(model, Mosek.Optimizer)
    # set_optimizer(model, Gurobi.Optimizer)
    set_optimizer(model, CPLEX.Optimizer)

    @variable(model, q_min[i]<=Q[i = 1:2, 1:n]<=q_max[i])
    # println("LB of Q",lower_bound.(Q))
    # println("UB of Q",upper_bound.(Q))
    # @variable(model, w[1:n, 1:n],Bin)
    @variable(model, p_to[1:2, 1:n])
    @variable(model, p_land[1:2, 1:n])
    @variable(model, t_1[1:n]>=0) 
    @variable(model, t_2[1:n]>=0) 
    @variable(model, t[1:n]>=0) 
    # todo: aとtの関係はconstraintで、
    @variable(model, T[i=2:n]>=0) # T=1は下に分離
    # println("about T", T)#ok!
    @variable(model, T_1>=0) 
    @variable(model, T_last>=0) 
    if have_timewindow
        @variable(model, U[1:2,1:n]>=0)
    end

    @objective(model, Min,  sum(t) + 1.000001*(T_1+sum(T)+T_last))

    # Qの列どうやってとる？？
    @constraint(model, c0[i = 1:n], t[i]<=a)
    @constraint(model, c1[j=1:n], [v_v*t_1[j];Q[:,j] - p_to[:,j]] in SecondOrderCone())
    @constraint(model, c2[j=1:n], [v_v*t_2[j];Q[:,j] - p_land[:,j]] in SecondOrderCone())
    @constraint(model, c3[j=1:n], [v_c*t[j];p_to[:,j] - p_land[:,j]] in SecondOrderCone())
    @constraint(model, c4, [v_c*T_1;p_o - p_to[:,1]] in SecondOrderCone())
    @constraint(model, c5[j = 2:n], [v_c*T[j];p_land[:,j-1] - p_to[:,j]] in SecondOrderCone())
    @constraint(model, c6,[v_c*T_last;p_f - p_land[:, n]] in SecondOrderCone())
    @constraint(model, c7[i = 1:n], (t_1[i]+t_2[i]<=t[i]))
    @constraint(model, c81[i = 1:n], (Q[1,i] == sum(w[i,j]*q[1,j] for j in 1:n)))
    @constraint(model, c82[i = 1:n], (Q[2,i] == sum(w[i,j]*q[2,j] for j in 1:n)))
    # @constraint(model, c9[i = 1:n], 1 == sum(w[i,:])) # for i
    # @constraint(model, c10[j = 1:n], 1 == sum(w[:, j])) # for j

    if have_timewindow
        @constraint(model, w1[i = 1:n], (T_1+t_1[1]+sum(t_2[j]+T[j+1]+t_1[j+1] for j in 1:i-1))<=U[2,i])
        @constraint(model, w2[i = 1:n], U[1,i]<=(T_1+t_1[1]+sum(t_2[j]+T[j+1]+t_1[j+1] for j in 1:i-1)))
        @constraint(model, w3[i = 1:n], (U[1,i] == sum(w[i,j]*u[1,j] for j in 1:n)))
        @constraint(model, w4[i = 1:n], (U[2,i] == sum(w[i,j]*u[2,j] for j in 1:n)))
    end

    set_time_limit_sec(model, 10800.0)
    elapsed_time = @elapsed optimize!(model)
    # @show termination_status(model)
    # @show objective_value(model)
    # @show elapsed_time # cputime
    score = 1.000
    try
        score = objective_value(model)
    catch
        score =infeasible_score # 10^14でかく#TODO: ここの値の調整必要
    end
    return score, p_to, p_land
end

function make_new_seq(seq)
    # ok
    # 近傍をとる関数
    # 2-optも入れる
    # parameta tuning
    if rand()>0.8
        n_seq = simple_swap(seq)
    else
        # n_seq = two_opt(seq)
        n_seq = simple_swap(seq)
        n_seq = simple_swap(n_seq)
    end
    return n_seq
end

function simple_swap(seq)
    i,j = sample(1:n, 2, replace=false)
    seq_i = seq[i,:]
    seq_j = seq[j,:]
    n_seq = copy(seq)
    n_seq[i,:] = seq_j
    n_seq[j,:] = seq_i
    return n_seq
end

function two_opt(seq)
    # seqを行列から巡る順番の配列に直す。
    vec = matrix_to_vector(seq)
    sz = length(vec)
    if sz<4 return seq end
    # 交換する場所を選択する
    x = rand(1:sz-1)
    y = rand(1:sz-1)
    while x== y
        y = rand(1:sz-1)
    end
    #実際に交換する。
    d = min(x,y)
    b = max(x,y)
    while(d+1<b-1)
        d+=1
        b-=1
        vec[d],vec[b] = vec[b], vec[d]
    end
    # 行列になおす
    matrix = vecter_to_matrix(vec)
    # return
    return matrix
end

function output( p_to, p_land, seq, q)
    p2 = deepcopy(p1)
    pto_value = value.(p_to)
    pland_value = value.(p_land)
    Q_value =make_Q(seq, q)
    tlvector = matrix_to_vector(seq)
    # intial to the first takeoff
    plot!(p2, [p_o[1], pto_value[1, 1]], [p_o[2], pto_value[2, 1]],
    color=:black)
    # last landing to the final 
    plot!(p2, [pland_value[1, n], p_f[1]], [pland_value[2,n], p_f[2]],
    color=:black)
    # plot!(pl, [targets[i,1], targets[i+1,1]], [targets[i,2],targets[i+1,2]], color=:black, linewidth=1, linestyle=:dot, label=\"\")

    for k = 1:n
        point = tlvector[k]
        print(point)
        # drone launch point
        scatter!(p2, [pto_value[1,k]], [pto_value[2, k]],
            color=:red)
        annotate!(p2, [pto_value[1, k] + 2], [pto_value[2, k]],
            text("t$(k)", color=:red, 10))
        # drone landing point
        scatter!(p2, [pland_value[1, k]], [pland_value[2, k]],
            color=:red, marker=:rect)
        annotate!(p2, [pland_value[1, k] + 2], [pland_value[2, k]],
            text("l$(k)", color=:red, 10))
        # ship-only route 
        plot!(p2, [pto_value[1, k], pland_value[1, k]],
            [pto_value[2, k], pland_value[2, k]],
            color=:red, linewidht=1, linestyle=:dot,
            label="")
        # # UAV-only route 
        #行き
        plot!(p2, [pto_value[1, k], Q_value[1, k]],
            [pto_value[2, k], Q_value[2, k]],
            color=:blue, linewidht=1, linestyle=:dot,
            label="")
        # 帰り
        plot!(p2, [Q_value[1, k], pland_value[1, k]],
            [Q_value[2, k], pland_value[2, k]],
            color=:blue, linewidht=1, linestyle=:dot,
            label="")

        # ship & UAV [between landing to next takeoff]
        if k < n
            plot!(p2, [pland_value[1, k], pto_value[1,k+1]],
                [pland_value[2, k], pto_value[2, k+1]],
                color=:black, linewidht=1, # linestyle=:dot,
                label="")
        end
    end
    display(p2)
    return p2
end

function output_simple_TSP( seq, q)
    """
    droneなしのTSPを可視化します。
    """
    p3 = deepcopy(p1)
    Q_value =make_Q(seq, q)
    tlvector = matrix_to_vector(seq)
    # intial to the first takeoff
    plot!(p3, [p_o[1], Q_value[1, 1]], [p_o[2], Q_value[2, 1]],
    color=:black)
    # last landing to the final 
    plot!(p3, [Q_value[1, n], p_f[1]], [Q_value[2,n], p_f[2]],
    color=:black)
    # plot!(pl, [targets[i,1], targets[i+1,1]], [targets[i,2],targets[i+1,2]], color=:black, linewidth=1, linestyle=:dot, label=\"\")

    for k = 1:n-1
        point = tlvector[k]
        # ship-only route 
        plot!(p3, [Q_value[1, k], Q_value[1, k+1]],
            [Q_value[2, k], Q_value[2, k+1]],
            color=:black, )
    end
    display(p3)
    return p3
end

function output_easy(seq, score)
    tlvector = Vector{Int}(undef, n) # takeoff and landing
    for i = 1:n,j = 1:n
        if value(seq[i, j]) > 0.9
            # This should be exactly 1 if there is no numerical error
            tlvector[i] = j
        end
    end
    println("objective_value : ", score)
    print("sequence ==> ")
    print(tlvector[1])
        for k = 2:n
            takeoff = tlvector[k]
            print("-> $takeoff ")
        end
    println()
end


function calctest()
    q,u = input()
    seq = make_first_seq(q,u)
    calc_score(seq, q, u)
end

function vecter_to_matrix(vector)
    """vectorには1~sizeまでの数字がひとつづつ入っていることを要求。"""
    n= length(vector)
    matrix = zeros(Int, n,n)
    for idx in vector
        matrix[idx, vector[idx]] = 1
    end
    return matrix
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
    println("vec", vector)
    return vector
end

function make_Q(seq, q)
    """
    Q:
    """
    println(seq)
    vec = matrix_to_vector(seq)
    println(vec)
    Q = zeros(Float64, 2,length(vec))
    for i in 1:length(vec)
        idx_Q = vec[i]
        Q[1, i] = q[1, idx_Q]
        Q[2,i] = q[2, idx_Q]
    end
    return Q
end

function calc_distance(q)
    """
    qの距離行列を返す.
    """
    n = size(q, 2)
    matrix = zeros(Float32, n,n)
    for i in 1:n
        for j in 1:n
            matrix[i,j] = norm(q[:,i] - q[:,j])
        end
    end
    return matrix
end

function calc_distance_with(q, point)
    n = size(q, 2)
    vector = zeros(Float32, n) 
    for i in 1:n
        vector[i] = norm(q[:,i] - point)
    end
    return vector
end

function save_figure(dir_path::String, pic_name::String, plt)
    """
    Plot instance pltをfile_pathとして保存する関数。
    """
    file_path = string(dir_path, "/", pic_name, ".png")
    println("save picture to: ", file_path)
    display(plt)
    savefig(file_path)
    print("save_complete!")
end

function write_to_csv(csv_path, n, compute_time, obj_value, n_iter, is_SA, is_infeasible,status,  data_path)
    seed = data_path[end-7:end-4] # seed値は必ず4桁を補償
    df  = DataFrame(n = [n], compute_time = [compute_time], obj_value = [obj_value], n_iter = [n_iter], is_SA = [is_SA], is_infeasible = [is_infeasible], status = [status],seed = [seed])
    CSV.write(csv_path, df, append = true, writeheader = false)
end


main()
# q,u = input()
# seq = make_first_seq(q,u)
# p3 = output_simple_TSP(seq, q)
# save_figure("temp/pictures", "easy_no_drone_answer", p3)