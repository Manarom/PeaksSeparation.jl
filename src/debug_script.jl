using Plots,BenchmarkTools,Revise,Peaks,Optimization,CSV,DataFrames,OptimizationOptimJL,ForwardDiff
using Profile
includet("PeaksSeparation.jl")

f(x,μ,σ,a) = exp(- 0.5* ^(x-μ,2)/(σ^2))*a
f_g(μ,σ,a) = x->f(x,μ,σ,a)
f_2p(x,a) = a[1] + a[2]*x + f(x,a[3:5]...) + f(x,a[6:8]...)
f_par(xval) = p->f_2p(xval,p)

par_check = [0.0;0.0;300.0;20.5;1.0;600.0;40.5;1.5]

mp = PeaksSeparation.MultiPeaks([300.0],[10.0],2)
PeaksSeparation.discr(par_check,mp)
mp.y
f1 = f_par(300.0) 
f1(par_check)

v_autodiff = ForwardDiff.gradient(f1,par_check)
v_an = similar(v_autodiff)

v_an = PeaksSeparation.fill_from_tuples!(v_an,PeaksSeparation.∇(300.0,1,mp))
Profile.init()
@profview for _ in 1:10000 
    PeaksSeparation.fill_from_tuples!(v_an,PeaksSeparation.∇(300.0,1,mp))
end
@benchmark PeaksSeparation.fill_from_tuples!(v_an,PeaksSeparation.∇(300.0,1,mp))
    #@profview
Profile.print() 
f_1() = 2+2
@profile f_1()
begin
    
    N=300
    test_data = Matrix{Float64}(undef,N,2)
    test_data[:,1] = collect(range(100,1000,N))

    mu_test = 300
    sigma_test = 300
    a_test = 0.5
    plot(test_data[:,1],f_g(mu_test,sigma_test,a_test).(test_data[:,1]),label="μ=$(mu_test),σ=$(sigma_test)")
    sigmas = (50, 20, 45)
    #sigmas = (50,)
    mus = (300, 600, 800)
    #mus = (300, )
    as = (1, 0.8, 0.5 )
    #as = (1,)
    d = @view test_data[:,2]
    fill!(d,0.0)
    x = @view test_data[:,1]
    for (m,s,a) in zip(mus, sigmas,as)
        @show (m,s,a)
        @. d += f_g(m,s,a).(x)
    end	
end





x = test_data[:,1];
y = test_data[:,2];
pl = plot(x,y)
mp = PeaksSeparation.MultiPeaks(x,y,3)
s_v = PeaksSeparation.fill_starting_vector!(Float64[],mp)
PeaksSeparation.fill_pars!(mp,s_v)
plot(mp)
gv = similar(s_v)

mp = PeaksSeparation.MultiPeaks(x,y,3,PeaksSeparation.LorentzPeak)
PeaksSeparation.fill_pars!(mp,s_v)
plot(mp)

(sol,p) = PeaksSeparation.fit_peaks(x,y,optimizer = ParticleSwarm(),N=3,use_constraints=true)
plot(p)

(sol,p_lorentz) = PeaksSeparation.fit_peaks(x,y,optimizer = LBFGS(),N=3,use_constraints=true,
                PeakType=PeaksSeparation.GaussPeak)
plot(p_lorentz)

(sol,p_lorentz) = PeaksSeparation.fit_peaks!(p_lorentz,optimizer=LBFGS(),N=3,use_constraints=true)
pl = plot(p_lorentz)
plot!(pl,x,d,label="data to fit")

# loadig test data
using Plots,BenchmarkTools,Revise,Peaks,Optimization,CSV,DataFrames,OptimizationOptimJL
include("NetzFileParser.jl")
nf = NetzFileParser.NetzFile(joinpath(@__DIR__,"data","kinetics","q5.txt"))
#data_to_fit = CSV.read(joinpath(@__DIR__,"data","q5.txt"),DataFrame)
p1 = plot(nf.data[:,1],nf.data[:,3])
PeaksSeparation.maxima_indices(nf.data[:,3],3)
s_v = PeaksSeparation.fill_starting_vector!(Float64[],p)

include("PeaksSeparation.jl")
(sol,p) = PeaksSeparation.fit_peaks(nf.data[:,1],nf.data[:,3],N=5,use_constraints=true)
(sol,p) = PeaksSeparation.fit_peaks(nf.data[:,1],nf.data[:,3],N=3,optimizer = NelderMead(),use_constraints=false)
p1 = plot(p.x,p.y0)
plot!(p1,p.x,p.y)
(sol,p) = PeaksSeparation.fit_peaks!(p,optimizer = NelderMead(),use_constraints=false)
#pl = plot(data_to_fit.T,data_to_fit.q20)
plot!(p1,p.x,p.y)

out_mat = PeaksSeparation.split_peaks(p,sort_by = :μ,rev=true)
for c in eachcol(out_mat)
    plot!(p1,p.x,c)
end
p1

plot(p)


