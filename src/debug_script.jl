using Plots,BenchmarkTools,Revise,Peaks,Optimization,CSV,DataFrames,OptimizationOptimJL
includet("PeaksSeparation.jl")

begin
    f_g(μ,σ,a)=x->exp(- ^(x-μ,2)/(2*σ^2))*a
    N=300
    test_data = Matrix{Float64}(undef,N,2)
    test_data[:,1] = collect(range(100,1000,N))

    mu_test = 300
    sigma_test = 3e-2
    a_test = 0.5
    plot(test_data[:,1],f_g(mu_test,sigma_test,a_test).(test_data[:,1]),label="μ=$(mu_test),σ=$(sigma_test)")
    sigmas = (50, 20, 45)
    mus = (300, 600, 800)
    as = (1, 0.8, 0.5 )
    d = @view test_data[:,2]
    fill!(d,0.0)
    x = @view test_data[:,1]
    for (m,s,a) in zip(mus, sigmas,as)
        @show (m,s,a)
        @. d += f_g(m,s,a)(x)
    end	
end



x = test_data[:,1];
y = test_data[:,2];
pl = plot(x,y)
mp = PeaksSeparation.MultiPeaks(x,y,3)
s_v = PeaksSeparation.fill_starting_vector!(Float64[],mp)
PeaksSeparation.fill_pars!(mp,s_v)
PeaksSeparation.residual!(mp)
plot!(pl,mp.x,mp.y)


(sol,p) = PeaksSeparation.fit_peaks(x,y,optimizer = NelderMead(),N=3,use_constraints=false)
pl = plot(p.x,p.y)
plot!(pl,p.x,p.y0,label="data to fit")

(sol,p) = PeaksSeparation.fit_peaks!(p)
pl = plot(p.x,p.y);
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


