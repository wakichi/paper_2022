using Plots
using LinearAlgebra
using Plots: plot, plot! # to supress VScode linter warningd
using CSV
using DataFrames
using Random
using StatsBase # 重複なく2つ抜き出すのに使用


ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
import Pkg
Pkg.add("CPLEX")
Pkg.build("CPLEX")
using JuMP
using CPLEX

# パラメータの設定(global)
is_includecsv = true
have_timewindow = false
is_plot = false

function main()
    # 入力
    q, u = input()
    global n = size(q, 2)
    # ヒューリスティックの実行部分
    score, seq = climing(q, u)
    # 出力
    # output_easy(seq, score)
end

function input()
    if is_includecsv
        test_file = "/Users/wakitakouhei/Lab/paper_2022/tests/points_hard.csv"
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
    if is_plot
        # plotする
    end
    return (q,u)
end

function climing(q,u)
    # 初期解の生成(どういう形？)(1 0 0 0;0 0 1 0)的な形
    seq = make_first_seq(q,u)
    score = calc_score(seq, q, u)
    for i = 1:1000
        # 近傍をとる。
        n_seq = make_new_seq(seq)
        # スコアの算出(cplexに投げる)
        n_score = calc_score(n_seq, q, u)
        # 順列を更新するか判断する。
        println("i:", i)
        output_easy(seq, score)
        output_easy(n_seq, n_score)
        if score>n_score # 最小化問題なので
            score = n_score
            seq = n_seq
        end
    end
    # 結果を可視化するならここで。
    return (score, seq)
end

function make_first_seq(q,u)
    # TODO:tspで初期回を生成する。
    # now:まずはq1, q2,...qnを順番に辿る
    n = size(q, 2) # the number of target points 
    seq = Matrix{Int32}(1I, (n, n)) # 単位行列
    return seq
end

function calc_score(w, q, u)
    # cplexを用いてsequenceからスコアを計算する。
    # スコアは単純に目的関数値を用いる。
    n = size(q, 2) # the number of target points 
    q_min = [0;0] # qの最小のx、y
    q_max = [60;60]
    v_c = 18 # km/h # default = 18 carrior speed
    v_v = 60 # km/h # default = 90 vehicle speed
    a = 21/60 # vehicle operation range
    p_o = [0;0] # startpoint
    p_f = [50;0]# end point

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

    @objective(model, Min,  sum(t) + (T_1+sum(T)+T_last))

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
    return objective_value(model)
end

function make_new_seq(seq)
    # ok
    # 近傍をとる関数
    i,j = sample(1:n, 2, replace=false)
    seq_i = seq[i,:]
    seq_j = seq[j,:]
    n_seq = copy(seq)
    n_seq[i,:] = seq_j
    n_seq[j,:] = seq_i
    return n_seq
end

function output(res)
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
main()