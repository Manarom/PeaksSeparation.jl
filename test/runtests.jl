using PeaksSeparation
using Test,ForwardDiff

@testset "PeaksSeparation.jl" begin
    # Write your tests here.
    for (f,grad) in zip((PeaksSeparation.f_gauss,PeaksSeparation.f_lorentz,PeaksSeparation.f_voigt),(PeaksSeparation.∇gauss,PeaksSeparation.∇lorentz,PeaksSeparation.∇voigt))
        f_p = p->f(400.0,p...) 
        par_2_check = (f == PeaksSeparation.f_voigt) ? [300.0,50.0,1.0,0.3] : [300.0,50.0,0.3]
        grad_auto = ForwardDiff.gradient(f_p,par_2_check)
        grad_calc = grad(400.0,par_2_check...)
        for i in eachindex(grad_auto)
            @show "$(f) at $(i)"
            @test grad_auto[i] ≈ grad_calc[i] 
        end
    end
end
