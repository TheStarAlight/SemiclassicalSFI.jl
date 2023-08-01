using SemiclassicalSFI
using Test
using StaticArrays
using OrdinaryDiffEq
using DiffEqGPU, CUDA

@info "# Testing trajectory simulation ..."

@testset verbose=true "Trajectory simulation" begin
    # l = Lasers.Cos2Laser(2e14,800,2,1.0)
    traj = function (u,p,t)
        F0 = 0.05338027007325633 # Lasers.LaserF0(l)
        ω = 0.05695419065625 # Lasers.AngFreq(l)
        N = 2   # l.cyc_num
        φ = 0.0 # l.cep
        ε = 1.0 # l.ellip
        tFx, tFy, tFz = -1*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
        du1 = u[4]
        du2 = u[5]
        du3 = u[6]
        du4 = tFx- F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 1/N*sin(ω*t/(2N))*cos(ω*t+φ))
        du5 = tFy- F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*cos(ω*t+φ) - 1/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
        du6 = tFz
        @SVector [du1,du2,du3,du4,du5,du6]
    end
    traj_t_final = 120.
    traj_rtol = 1e-6
    init = [10.0  9.0  8.0  7.0;    # x
             5.0  1.0  0.0  0.0;    # y
             0.0  0.0  6.0  1.0;    # z
             0.0  1.0  0.0  1.0;    # px
             1.0  1.0 -1.0  2.0;    # py
             0.0  0.0  0.0  3.0;    # pz
             0.0 -5.0  0.0  0.0;]
    final = [-107.586659754  14.769525931  -111.052945103   10.541838348;
               88.693094313  66.732075229  -135.011048697  216.489086428;
                0.0           0.0             0.580938720  357.985889110;
               -0.972143888  -0.017934642    -0.992390595    0.028702222;
                0.852999481   0.666821600    -0.932762255    1.985858411;
                0.0           0.0            -0.047583468    2.972879185;]
    traj_ODE_prob::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,6)), (0,traj_t_final))
    init_traj = (prob,i,repeat) -> remake(prob; u0=SVector{6}([init[k,i] for k in 1:6]), tspan=(init[7,i],traj_t_final))
    ensemble_prob::EnsembleProblem = EnsembleProblem(traj_ODE_prob, prob_func=init_traj, safetycopy=false)
    solc = nothing
    solg = nothing
    @info "Testing CPU..."
    @testset verbose=true "CPU" begin
        @test begin
            solc = solve(ensemble_prob, OrdinaryDiffEq.Tsit5(), EnsembleThreads(), trajectories=size(init,2), adaptive=true, dt=0.01, reltol=traj_rtol, save_everystep=false);
            mapreduce(function((k,i),) ≈(final[k,i],solc[i][end][k],rtol=1e-2) end, *, [(k,i) for k in 1:6, i in 1:size(init,2)])
        end
    end
    @info "Testing GPU..."
    @testset verbose=true "GPU" begin
        @test begin
            solg = solve(ensemble_prob, DiffEqGPU.GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend(),0.0), trajectories=size(init,2), adaptive=true, dt=0.01, reltol=traj_rtol, save_everystep=false);
            mapreduce(function((k,i),) ≈(final[k,i],solg[i][end][k],rtol=1e-2) end, *, [(k,i) for k in 1:6, i in 1:size(init,2)])
        end
    end
end