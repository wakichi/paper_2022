## Packages 

ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
# import Pkg
# Pkg.add("CPLEX")
# Pkg.build("CPLEX")

using JuMP
# using MosekTools
# using Gurobi
using CPLEX
using Plots
using LinearAlgebra
using Plots: plot, plot! # to supress VScode linter warningd
using CSV
using DataFrames

gr()
# unicodeplots()

## Input parameters

# target points
# The results of this data is as follows:    
#    objective_value(model) = 6.653748364493911
#    elapsed_time = 2.009810452
#    sequence ==> 1 -> 1, 2 -> 3, 4 -> 6, 7 -> 7, 8 -> 8, 9 -> 9, 10 -> 10,


# include("points-n-0020-seed-1024.jl")    # test data by "generate-instance.jl"
# 実行時パラメータを指定する場所
is_includecsv = true
have_timewindow = true

if is_includecsv
    test_file = "/Users/wakitakouhei/Lab/paper_2022/tests/points-n-0010-seed-1029.csv"
    println("loading the data in $(test_file)")
    df = CSV.read(test_file, DataFrame; header=0)    # test data by "generate-instance.jl"
    q = [df[i, j] for i = 1:2, j = 1:size(df, 2)]
    # TODO: 3,4行目が存在しない場合、printした上でhave_timewindowをfalseにする。
    if have_timewindow
        u = [df[i, j] for i = 3:4, j = 1:size(df, 2)]
    end
end

# constants
# 数値実験欄見て埋める
# qは簡単のため4点程度で行う。
# uの可視化どうする？？？
# 案１. アノテーションのように[u1,u2]と点の近くに書く
# q = [0 0 1 46 48 46 50 50 20 30
# 40 50 49 50 48 46 48 35 5 0]
if !is_includecsv
    q = [45 10 30 
     45 40 20 ]
    u = [0 0 200 
        200 200 300 ]
end

n = size(q, 2) # the number of target points 
q_min = [0;0] # qの最小のx、y
q_max = [60;60]
v_c = 18 # km/h # default = 18 carrior speed
v_v = 60 # km/h # default = 90 vehicle speed
a = 21/60 # vehicle operation range
p_o = [0;0] # startpoint
p_f = [50;0]# end point
#todo q とuのサイズ違う場合をチェック


# q_s = [0; 0]   # initial point   
# q_f = [50; 0]  # final point 
# v_c = 18 # km/h # default = 18
# v_h = 90 # km/h # default = 90
# t_hmax = 21 / 60 # 21 minutes = 21/60 hours
# n = size(q, 2) # the number of target points 

# Plot initial locations

p1 = plot(xlims=(-5, 55), ylims=(-5, 55), legend=:none, aspect_ratio=:equal)
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

## distance matrix
# this part may be slow, but this part is not the main bottleneck,
# so we leave this part as it is.
# d = zeros(n, n)
# for i = 1:n, j = i+1:n
#     d[i, j] = sum([norm(q[:, k] - q[:, k+1]) for k = i:j-1])
# end

## big-M
##Here
# M = maximum(d) + 2 * 50 * sqrt(2) # diagonal length of the area


## make the model
model = Model()
# set_optimizer(model, Mosek.Optimizer)
# set_optimizer(model, Gurobi.Optimizer)
set_optimizer(model, CPLEX.Optimizer)

@variable(model, q_min[i]<=Q[i = 1:2, 1:n]<=q_max[i])
# println("LB of Q",lower_bound.(Q))
# println("UB of Q",upper_bound.(Q))
@variable(model, w[1:n, 1:n],Bin)
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
@constraint(model, c9[i = 1:n], 1 == sum(w[i,:])) # for i
@constraint(model, c10[j = 1:n], 1 == sum(w[:, j])) # for j

if have_timewindow
    @constraint(model, w1[i = 1:n], (T_1+t_1[1]+sum(t_2[j]+T[j+1]+t_1[j+1] for j in 1:i-1))<=U[2,i])
    @constraint(model, w2[i = 1:n], U[1,i]<=(T_1+t_1[1]+sum(t_2[j]+T[j+1]+t_1[j+1] for j in 1:i-1)))
    @constraint(model, w3[i = 1:n], (U[1,i] == sum(w[i,j]*u[1,j] for j in 1:n)))
    @constraint(model, w4[i = 1:n], (U[2,i] == sum(w[i,j]*u[2,j] for j in 1:n)))
end

# @variable(model, alpha[i=1:n, j=i:n], Bin)
# @variable(model, f[i=1:n, j=i:n] >= 0)
# @variable(model, s[i=1:n, j=i:n-1] >= 0)

# @constraint(model, e18[i=1:n, j=i:n],
#     f[i, j] - t_hmax <= M * (1 - alpha[i, j]))

# @variable(model, tau[1:2, 1:n])
# @variable(model, l[1:2, 1:n])
# @constraint(model, e19a[i=1:n, j=i:n],
#     [v_c * f[i, j] + M * (1 - alpha[i, j]); tau[:, i] - l[:, j]]
#     in
#     SecondOrderCone())

# @variable(model, tau_q[1:n] >= 0)
# @variable(model, q_l[1:n] >= 0)
# @constraint(model, tau_q_socp[i=1:n],
#     [tau_q[i]; tau[:, i] - q[:, i]] in SecondOrderCone())
# @constraint(model, q_l_socp[j=1:n],
#     [q_l[j]; q[:, j] - l[:, j]] in SecondOrderCone())
# @constraint(model, e19b[i=1:n, j=i:n],
#     tau_q[i] + d[i, j] + q_l[j] - v_h * f[i, j] <= M * (1 - alpha[i, j]))
# @constraint(model, e19c[i=1:n, j=i:n-1],
#     [v_c * s[i, j] + M * (1 - alpha[i, j]); l[:, j] - tau[:, j+1]] in SecondOrderCone())

# @variable(model, qs_tau >= 0)
# @variable(model, ln_qf >= 0)
# @constraint(model, qs_tau_socp, [qs_tau; q_s - tau[:, 1]] in SecondOrderCone())
# @constraint(model, ln_qf_socp, [ln_qf; l[:, n] - q_f] in SecondOrderCone())

# @objective(model, Min, (1 / v_c) * (qs_tau + ln_qf) + sum(f) + sum(s))

# @constraint(model, e2[k=1:n],
#     sum([alpha[i, j] for i = 1:k, j = k:n]) == 1)

## optimize the model
# set_silent(model) # to supress detailed log of MISOCP solver
# set_time_limit_sec(model, 3600.0)
set_time_limit_sec(model, 10800.0)
elapsed_time = @elapsed optimize!(model)
@show termination_status(model)
@show objective_value(model)
@show elapsed_time # cputime


println(value.(p_to))
println(value.(T_1))

## print flight sequence
# alpha_value = zeros(Int,n,n)
tlvector = Vector{Int}(undef, n) # takeoff and landing

for i = 1:n,j = 1:n
    global tlvector
    if value(w[i, j]) > 0.9
        # This should be exactly 1 if there is no numerical error
        tlvector[i] = j
        
    end
end

print("sequence ==> ")
print(tlvector[1])
for k = 2:n
    takeoff = tlvector[k]
    print("-> $takeoff ")
end

if 1==1
## Plot the result
p2 = deepcopy(p1)
pto_value = value.(p_to)
pland_value = value.(p_land)
Q_value = value.(Q)
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
savefig("/Users/wakitakouhei/Lab/paper_2022/hoge.png")
end
