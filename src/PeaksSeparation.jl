module PeaksSeparation
    using Optimization, Peaks, StaticArrays,
        OptimizationOptimJL,ForwardDiff,RecipesBase,
        AllocCheck,Statistics,LinearAlgebra,
        QuadGK

    export MultiPeaks,
            fit_peaks,
            fit_peaks!,
            split_peaks,
            GaussPeak,
            LorentzPeak,
            VoigtPeak,
            ConstantBaseLine,
            LinearBaseLine,
            QuadraticBaseLine,
            CubicBaseLine,
            WeibullPeak

    const default_optimizer = Ref(ParticleSwarm(n_particles=100))
    const fwhm = 2*sqrt(2*log(2.0)) # contant to evaluate the dispersion of Gaussian from peaks's fwhm
    include("Utils.jl")
    sqr(x) = x^2
    f_gauss(x,μ,σ,a) = abs(σ)>1e-16 ? a*exp(- 0.5*^(x-μ,2) * ^(σ,-2)) : a*exp(- ^(x - μ, 2) * 0.5 * 1e32)
     
   
    function ∇gauss(x,μ,σ,a) 
        abs(σ) > 1e-16 || (σ = 1e-16) 
        dfda =  f_gauss(x, μ, σ, 1.0) 
        dfdmu =  a * (x - μ) * dfda /( σ^2 )
        return (dfdmu, dfdmu*(x-μ)/σ ,dfda )
    end

    f_lorentz(x,μ,σ,a)  = abs(σ)>1e-16 ? a / (1.0 + ((x - μ) / σ) ^ 2) : a / (1.0 + 1e32*(x - μ) ^ 2)
    function ∇lorentz(x,μ,σ,a) 
        abs(σ) > 1e-16 || (σ = 1e-16) 
        dfda = f_lorentz(x,μ,σ,1.0)
        dfdmu = 2*a*((x - μ) / σ^2)*dfda^2
        return (dfdmu, dfdmu*(x-μ)/σ ,dfda )
    end

    f_voigt(x,μ,σ,a,b) = f_gauss(x,μ,σ,a)*(1-b) + f_lorentz(x,μ,σ,a)*b
    function ∇voigt(x,μ,σ,a,b) 
        abs(σ) > 1e-16 || (σ = 1e-16) 
        fl1 = f_lorentz(x, μ, σ, 1.0) 
        fg1 = f_gauss(x, μ, σ, 1.0) 
        dfdmu = a * (x - μ) * (fg1 * (1 - b)  + 2 * b * fl1^2) /( σ^2 )
        return (dfdmu, 
                dfdmu*(x-μ)/σ, #dfds
                fg1 * (1 - b) + fl1 * b, #dfda
                a*(fl1 - fg1) # dfdb 
                )
    end

    function f_weibull(x,μ,σ,a,b)
      (x - μ) > 0 || return 0.0
      return  abs( σ ) > 1e-16 ? a * ^((x - μ)/σ, b - 1) * exp(- ^((x - μ)/σ , b) ) : a*^(1e16*(x - μ), b - 1)*exp(- ^( 1e16*(x-μ) ,b) ) 
    end
    function ∇weibull(x,mu,s,a,b)
        #error("Gradient is under construction")
        (x - mu) > 0 || return (0.0,0.0,0.0,0.0)
        x_sn = (x - mu)/s
        x_b = x_sn^(b - 1)
        e_b = exp(- x_sn^b )*x_b
        dfdmu = a * e_b * (b * x_b  - (b - 1)/x_sn)/s
        (
            dfdmu,
            (x - mu) * dfdmu/s,
            e_b,
            a*log(x_sn) * e_b * ( 1 -  x_sn * x_b)
        )
    end
    # N - is the number of parameters of this component
    abstract type AbstractCurveComponent{N} end
  
eval_to!(y,cc::AbstractCurveComponent,x) = map!(cc,y,x)
    """
    add_to!(y,cc::AbstractCurveComponent,x)

 Adds result of curve component (like peak or baseline) evaluation to the vector y
"""
add_to!(y,cc::AbstractCurveComponent,x) = @. y += cc(x)
add_to!(y,cc::AbstractCurveComponent,x,bynary_operator::Function) = @. y = bynary_operator(y,cc(x))
    function add_to!(y,curve_component_iterator,x)
        for cc in curve_component_iterator 
            add_to!(y,cc,x)
        end
        return y
    end
    names(cc::AbstractCurveComponent) = getfield(cc,:names)
    parameters(cc::AbstractCurveComponent) = getfield(cc,:parameters)
    flags(cc::AbstractCurveComponent) = getfield(cc,:flags)
    parnumber(::AbstractCurveComponent{N}) where N = N
    Base.length(cc::AbstractCurveComponent) = parnumber(cc)

    parnumber(::Type{<:AbstractCurveComponent{N}}) where N = N
    parindex(cc::AbstractCurveComponent{N},field::Symbol) where N = begin
        nm = names(cc)
        for i in 1:N
            field != nm[i] || return i
        end
        return 0
    end
    function Base.getproperty(cc::AbstractCurveComponent{N},field::Symbol) where N
        !(field == :parameters || field == :names || field == :flags) || return getfield(cc,field)
        nm = names(cc)
        for i in 1:N
            field != nm[i] || return parameters(cc)[i]
        end
        getfield(cc,field)
    end
    # baseline supertype
    abstract type AbstractBaseLine{N}<:AbstractCurveComponent{N} end
    abstract type PolynomialBaseLine{N}<:AbstractBaseLine{N} end
    abstract type AbstractPeak{N}<:AbstractCurveComponent{N} end
    """
    peaks_distance(p1::AbstractCurveComponent,p2::AbstractCurveComponent)

Evaluates Euclidean distance between two curves (including arbitrary shape peaks)
as ∫(p₁(x)-p₂(x))²dx/(∫p₁²(x)dx + ∫p₂²(x)dx) = I₁₂/(I₁ + I₂)

retuns a tupe of distance and error functions thich is evaluated as:

error = √(e₁₂/(I₁ + I₂))² + (e₁² + e₂²)/(I₁ + I₂)⁴ 

"""
function peaks_distance(p1::T,p2::P) where {T<:AbstractCurveComponent,P<:AbstractCurveComponent}
        f1 = x->^(p1(x),2)
        f2 = x->^(p2(x),2)
        f  = x-> ^(p1(x)-p2(x),2)
        tasks = Vector{Task}(undef, 3)
        out = Vector{NTuple{2,Float64}}(undef,3)
        Threads.@sync begin
            for (i,fi) in enumerate((f,f1,f2))
                tasks[i] = Threads.@spawn QuadGK.quadgk(fi,-Inf,Inf)
            end
        end
        map!(fetch,out,tasks)
        (f1_f2,e1_e2) = out[1] 
        (f1_m,    e1) = out[2]
        (f2_m,    e2) = out[3]
        e1 = e1^2
        e2 = e2^2
        dist = f1_f2/(f1_m + f2_m) # integral Euclidean-like distance between functions
        e = sqrt(sum(t->t^2,
                        (e1_e2/(f1_m + f2_m) ,(e1+e2)/^(f1_m + f2_m,2)) # error estimation 
                    )) # error estimation 
        return (dist, e)   
    end
    """
    fill_pars!(curve_component_iterator,itr)

curve_component_iterator - iterable object of curve_components
itr - statefull iterator!
"""
function fill_pars!(curve_component_iterator,itr)
        for cc in curve_component_iterator
            fill_pars!(cc,itr)
        end
    end
    """
    fill_pars!(::PolynomialBaseLine,x)

function to fill baseline polynomial coefficients from the iterable object x

"""
function fill_pars!(::PolynomialBaseLine,x) end
# generating functions
for i in 1:5 # reserving function for up to five 
        @eval function fill_pars!(bl::PolynomialBaseLine{$i},x) 
                setfield!(bl,:parameters , @iterator_unpack($i,x))
                return bl
        end
        @eval function fill_pars!(p::AbstractPeak{$i},x)
            setfield!(p,:parameters, @iterator_unpack($i,x))
        end
        @eval function ∇(x,::PolynomialBaseLine{$i})
            ntuple(i->^(x,i-1),Val($i))
        end
    end

    (bl::PolynomialBaseLine)(x) = @evalpoly(x,getfield(bl,:parameters)...)
    # generating concrete polynomial baseline type up to 0...3
    _names = Symbol[]
    _pars = Float64[] 
    _flags = Bool[]
    for (i,bl_name) in enumerate((:ConstantBaseLine,:LinearBaseLine,:QuadraticBaseLine,:CubicBaseLine))
        push!(_names,Symbol("b$(i)"))
        push!(_pars,0.0)
        push!(_flags,false)
        args = (Tuple(_pars),Tuple(_names),Tuple(_flags))
        @eval mutable struct $bl_name <: PolynomialBaseLine{$i}
        parameters::NTuple{$i,Float64}
        names::NTuple{$i,Symbol}
        flags::NTuple{$i,Bool}
        $bl_name() = begin 
                new($args...)
            end
        end
    end
    #_names = nothing;_pars = nothing;_flags = nothing
    # concrete peaks types 
    # gauss and lorenz peaks types 
   for peak_name in (:GaussPeak,:LorentzPeak)
        @eval mutable struct $peak_name <:AbstractPeak{3}
            parameters::NTuple{3,Float64}
            names::NTuple{3,Symbol}
            flags::NTuple{3,Bool}
            index::Int
            $peak_name(i::Int=1,mu=0.0,sigma=1.0,a=1.0) = begin 
                new((mu,sigma,a),(:μ,:σ,:a),(false,false,false),i)
            end
        end
   end
   #pseudo-voigt peak
   for peak_name in (:VoigtPeak,:WeibullPeak)
    @eval mutable struct $peak_name <:AbstractPeak{4}
        parameters::NTuple{4,Float64}
        names::NTuple{4,Symbol}
        flags::NTuple{4,Bool}
        index::Int
            $peak_name(i::Int=1,mu=0.0,sigma=1.0,a=1.0,b=0.5) = begin 
                new((mu,sigma,a,b),(:μ,:σ,:a,:b),(false,false,false,false),i)
            end
        end
    end
    # generating functions and gradients
   for (peak_name,f_name,grad_name) in zip((:GaussPeak,:LorentzPeak,:VoigtPeak,WeibullPeak),
                            (:f_gauss, :f_lorentz, :f_voigt, :f_weibull),
                            (:∇gauss,  :∇lorentz , :∇voigt, :∇weibull))
        #parnmb = parnumber(eval(peak_name))
        @eval (p::$peak_name)(x::T) where T= $f_name(x,getfield(p,:parameters)...)::T
        @eval ∇(x,p::$peak_name)= $grad_name(x,getfield(p,:parameters)...)
   end
   is_has_gradient(::Type{P}) where P <: AbstractPeak = P<:GaussPeak || P<:LorentzPeak || P<: VoigtPeak || P<: WeibullPeak
   is_has_gradient(::P) where P<:AbstractPeak = is_has_gradient(P)
   is_has_hessian(::Type{P}) where P<:AbstractPeak = false
   is_has_hessian(::AbstractPeak) = false
 #peak functions   
    mutable struct MultiPeaks{N,B<:AbstractBaseLine,P<:AbstractPeak,T<:AbstractVector}
        x::T
        y0::T
        y::T
        r::T
        peaks::SVector{N,P}
        flags::SVector{N,Bool}
        baseline::B
        """
    MultiPeaks(x::T,y0::T,N::Int,
                    ::Type{P}=GaussPeak,
                    ::Type{B}=LinearBaseline) where {T<: AbstractVector,P<:APeak,B<:BaseLine}

x - coordinate
y0 - experimental curve
P - type of peaks 
B - type of baseline function 
"""
MultiPeaks(x::T,y0::T,N::Int,
                ::Type{P}=GaussPeak,
                ::Type{B}=LinearBaseLine) where {P<:AbstractPeak,B<:AbstractBaseLine,T}= begin
                            @assert N>=1 "Number of peaks must be greater or equal 1"
                            mp = new{N,B,P,typeof(x)}(x,y0,
                                        similar(y0),similar(y0),
                                        SVector(ntuple(P,Val(N))),
                                        SVector(ntuple(_->false,Val(N))),
                                        B())
                            residual!(mp)
                            return mp
        end
    end

function Base.getindex(mp::MultiPeaks{N},i) where N
    (i >N || i < 0) ? error("Index out of bounds") : return i == 0 ? getfield(mp,:baseline) : getfield(mp,:peaks)[i]
end
"""
    is_new_params(p,mp::MultiPeaks{N,B,P}) where {N,B,P}

Checks if parameters in vector p are the same as currently in settled in for 
this object (used to check if the gradient should be recalculated)
"""
function is_new_params(p,mp::MultiPeaks{N,B,P}) where {N,B,P}
    length(p) == parnumber(mp) || return true
    Base_pars_number = parnumber(B)
    Peak_pars_number = parnumber(P)
    _pb = @view p[1:Base_pars_number]
    for (bj,pj) in zip(parameters(mp[0]),_pb)
        bj == pj || return true 
    end
    startind = Base_pars_number + 1
    endind = 0 
    for ith_peak in mp.peaks
        endind = startind + Peak_pars_number -1
        _pb = @view p[startind:endind]
        for (bj,pj) in zip(parameters(ith_peak),_pb)
            bj == pj || return true 
        end
        startind = endind + 1
    end
    return false
end
"""
    parnumber(::MultiPeaks{N,B,P}) where {N,B,P}

Total number of parameters used to fit the curve 
"""
parnumber(::MultiPeaks{N,B,P}) where {N,B,P} = parnumber(B)+ N*parnumber(P) # total parnumber of the curve 
peaknumber(::MultiPeaks{N}) where N = N
#
    """
    add_baseline!(mp::MultiPeaks)

Adds baseline to fitting results
"""
add_baseline!(mp::MultiPeaks) = add_to!(mp.y,mp.baseline,mp.x)

    """
    add_peaks!(mp::MultiPeaks)

Adds peaks to fitting results
"""
add_peaks!(mp::MultiPeaks) = add_to!(mp.y,mp.peaks,mp.x)
    """
    fill_y!(mp::MultiPeaks)

Evaluates curve fitter and fills y field 
"""
function fill_y!(mp::MultiPeaks) 
        fill!(mp.y,0.0)
        add_baseline!(mp)
        add_peaks!(mp)
        return mp
    end
    """
    residual!(mp::MultiPeaks{N}) where N

Evaluates and fills the residual vector
"""
function residual!(mp::MultiPeaks{N}) where N
        fill_y!(mp)
        copyto!(mp.r,mp.y)
        @. mp.r -= mp.y0
        return mp
    end

    """
    fill_pars!(mp::MultiPeaks{N},x) where N

Fills all curve components parameters from iterable object x

"""
function fill_pars!(mp::MultiPeaks{N},x) where N
        itr = Iterators.Stateful(x)
        fill_pars!(getfield(mp,:baseline),itr)
        fill_pars!(getfield(mp,:peaks),itr)
        residual!(mp) # recalculating residuals
        return nothing
    end
    """
    discr(x,mp::MultiPeaks{N}) where N

Evaluates discrepancy value 
"""
function discr(x,mp::MultiPeaks{N}) where N
        fill_pars!(mp,x)
        return sum(sqr,mp.r)/2
    end

    discr_fun(x,p) = x->discr(x,p)
    """
    ∇(xi,ri,mp::MultiPeaks{N}) where N

Gradient of fitter with respect to all fitting variables, returns a Tuple of Tuples 
each tuple contains the gradient value (xi and ri are supposed to be a single number)  
"""
∇(xi,ri,mp::MultiPeaks{N}) where N =  ntuple(i->tuple_mult(∇(xi,mp[i-1]),ri), Val(N+1))

"""
    grad!(g::AbstractVector, x ,mp::MultiPeaks{N,B,P}) where {N,B,P}

Fills in-pace calculated gradient to g vector 
"""
function grad!(g::AbstractVector, x ,mp::MultiPeaks{N,B,P}) where {N,B,P}
        !is_new_params(x,mp) || fill_pars!(mp,x) 
        fill!(g,0.0)
        for (i,xi) in enumerate(mp.x)
            ri = getindex(getfield(mp,:r),i)
            add_from_tuples!(g,∇(xi,ri,mp))    
        end
    end

function hess!(h, x, mp::MultiPeaks{N,B,P}) where {N,B,P}
        !is_new_params(x,mp) || fill_pars!(mp,x) 
        npoints = length(mp.x)
        npars = parnumber(mp)
        if is_has_gradient(P)
            grad_fill! = (g,x) -> fill_from_tuples!(g,∇(x,1.0,mp)) 
        else # if gradient is not provided 
            #   TBW
 #           grad_fill! = (g,x) -> copyto!(g, ForwardDiff.gradient())
        end
        J = Matrix{Float64}(undef,npoints,npars)
        #H = Matrix{Float64}(undef,npars,npars)
        for (i,xi) in enumerate(mp.x)
            ji = @view J[i,:]
            grad_fill!(ji,xi)
        end
        # calculating the approximate hessian J'*J
        mul!(h,transpose(J),J)
        return nothing
end    
 
    """
    fill_from_tuples!(v,  args::Vararg;resizable::Bool=false)

Puts all content of args (tuple of tuple or vector of vectors et.c) into a single vector
"""
    function fill_vector_with_pars!(v,p::MultiPeaks{N};resizable::Bool=false) where N 
        fill_from_tuples!(v,ntuple(i->parameters(p[i-1]), Val(N+1)),resizable=resizable)
    end

    """
    maxima_indices(y,N::Union{Int,Nothing})

Returns indices of local maxima of y vector if N=nothing than the number of peaks 
is obtained using `argmaxima` function from `Peaks` package, if N is fixed returns N 
indeices of peaks in y, if N>L (L is the number of peaks obtained by `argmaxima`) than
repeates indices of peaks with the largest height, if N<L takes N peaks with the largest 
height.
"""
function maxima_indices(y,N::Union{Int,Nothing})
        peaks_inds = Peaks.argmaxima(y) 
        L = length(peaks_inds)
        !isnothing(N) || return peaks_inds
        N != L || return peaks_inds
        if N<L # we need smaller ammount of peaks thn the actua; one
            idxs = partialsortperm(view(y,peaks_inds),1:N,rev=true)#sortperm(peaks.heights)
            return peaks_inds[idxs]
        end
        # we need more peaks !!!
        idxs = similar(peaks_inds,N)
        copyto!(idxs,peaks_inds);
        v = @view idxs[L+1 : N]
        map!(_->rand(peaks_inds),v,1:N-L)
        return idxs
    end
    """
    param_estimate(::Type{<:AbstractPeak},x,y,indices::Vector{Int})

Finds peaks and estimates and estimates the starting values of peaks location and width using Peaks.jl 
package function peakwidths.

Input variables:
x coordinate
y values 
indices  - indices in x vector of peaks maxima location

Returns tuple of tuples with (peak location,peaks width, peak amplitude) the number of 
tuples is the same as the number indices  
"""
function param_estimate(::Type{<:AbstractPeak},x,y,indices::Vector{Int}) end
function param_estimate(::Type{<:Union{GaussPeak,LorentzPeak}},x,y,indices::Vector{Int})
        N = length(indices)
        as = view(y,indices)
        npoints = length(x)
        cor = abs(-(extrema(x)...))/npoints
        (inds,sigmas) = peakwidths(indices,y,0.5*as)
        if any(isnan,sigmas)
            replace_nans!(sigmas,npoints/(N+1))
        end
        sigmas *=cor/fwhm
        mus = @view x[indices]
        return ntuple(i->(mus[i],sigmas[i],as[i]),N)
    end

    """
    param_estimate(T::Type{<:AbstractPeak},x,y,N::Union{Int,Nothing})


The same as [`param_estimate(::Type{<:Union{GaussPeak,LorentzPeak}},x,y,indices::Vector{Int})`](@ref) but 
the number of peaks is specified by N
"""
function param_estimate(T::Type{<:AbstractPeak},x,y,N::Union{Int,Nothing})  throw(DomainError("Undefined for $(T)")) end
param_estimate(T::Type{<:Union{GaussPeak,LorentzPeak}},x,y,N::Union{Int,Nothing}) = param_estimate(T,x,y,maxima_indices(y,N))
param_estimate(::Type{VoigtPeak},x,y,N::Union{Int,Nothing}) = (param_estimate(GaussPeak,x,y,maxima_indices(y,N))...,0.5)
param_estimate(::Type{WeibullPeak},x,y,N::Union{Int,Nothing}) = (param_estimate(GaussPeak,x,y,maxima_indices(y,N))...,1.0)
    
    """
    fill_starting_vector!(starting_vector,mp::MultiPeaks{PeakNum,BaseType,PeakType}) where {PeakNum,BaseType,PeakType}

Fills starting vector for the optimization by 
"""
function fill_starting_vector!(starting_vector,mp::MultiPeaks{PeakNum,BaseType,PeakType}) where {PeakNum,BaseType,PeakType}
        length(starting_vector) == parnumber(mp) || resize!(starting_vector,parnumber(mp))
        base_param_number = parnumber(BaseType)
        miny = minimum(mp.y0)
        return fill_from_tuples!(starting_vector,  
                    ((miny,),ntuple(_->0.0,base_param_number-1),
                    param_estimate(PeakType,mp.x,mp.y0,PeakNum)...)
                    ;resizable=false)
    end
    int_floor = Int ∘ floor

    """
    fit_peaks!(p::MultiPeaks;kwargs...)

Refits peaks returns optimization solution 
"""
function fit_peaks!(p::MultiPeaks;kwargs...)
        starting_vector =  fill_vector_with_pars!(Float64[],p,resizable=true)
        return fit_peaks(       nothing,
                                nothing;
                                p=p,
                                starting_vector = starting_vector,
                                optimizer = NelderMead(),
                                use_constraints =false,
                                kwargs...)
    end

    """
    fit_peaks(x::T,y::T ;
                        N::Union{Int,Nothing}=nothing,
                        PeakType::Type{<:AbstractPeak}=GaussPeak,
                        BaseLineType::Type{<:AbstractBaseLine}=LinearBaseLine,
                        starting_vector::Union{Nothing,AbstractVector} = nothing,
                        optimizer = nothing,
                        p::Union{Nothing,MultiPeaks} = nothing,
                        use_constraints::Bool=true,
                        constraints_expansion::Float64=1.0,
                        allow_negative_peaks_amplitude::Bool = true,
                        allow_negative_baseline_tangent::Bool = true) where T<:Union{Nothing,AbstractVector}

x - coordinates
y - curve
PeaksType  - all peaks in curve are of the same type
BaseLineType  - type of base line 
starting_vector - staring parameters can be provided esternally
optimizer - any supported optimizer from Optimization.jl package which does not use 
            the second derivative
p - can refit MultiPeaks object (if it is provided it prevails all other arguments which can be 
taken from it like x,y,starting_vector, peaks type and baseline type. 
use_constraints - true if box_constraints are used
constraints_expansion - coefficient of constraintes expansion (shrinks of dilates the box)
allow_negative_peaks_amplitude - if true all amplitude
allow_negative_baseline_tangent - if false tangent is supposed to gradually increase
"""
function fit_peaks(x::T,y::T ;
                        N::Union{Int,Nothing}=nothing,
                        PeakType::Type{<:AbstractPeak}=GaussPeak,
                        BaseLineType::Type{<:AbstractBaseLine}=LinearBaseLine,
                        starting_vector::Union{Nothing,AbstractVector} = nothing,
                        optimizer = nothing,
                        p::Union{Nothing,MultiPeaks} = nothing,
                        use_constraints::Bool=true,
                        constraints_expansion::Float64=1.0,
                        allow_negative_peaks_amplitude::Bool = true,
                        allow_negative_baseline_tangent::Bool = true) where T<:Union{Nothing,AbstractVector}
        
        is_p_provided = !isnothing(p)
        !is_p_provided || (N=peaknumber(p))
        #is_fixed_peaks_number = !isnothing(N)
        
        if isnothing(starting_vector)
            #pks = fix_peaks_number(y,N)
            #N = length(pks.indices)
            if !is_p_provided
                p = MultiPeaks(x,y,N,PeakType,BaseLineType)
            end
            starting = fill_starting_vector!(Float64[],p)
        else
            starting = starting_vector
        end
        fill_pars!(p,starting)
        optim_fun = OptimizationFunction(discr,grad=grad!,hess = hess!) #AutoFiniteDiff()
        
        if use_constraints
            box = box_constraints(p,
                            constraints_expansion=constraints_expansion,
                            allow_negative_peaks_amplitude=allow_negative_peaks_amplitude,
                            allow_negative_baseline_tangent = allow_negative_baseline_tangent)
            check_bounds_on_starting_vector!(starting,box)
            probl= OptimizationProblem(optim_fun, 
                                        starting,
                                        p,
                                        lb=box.lb,
                                        ub=box.ub,
                                        reltol=1e-8,
                                        maxiters = 1e4) 
        else
            probl= OptimizationProblem(optim_fun, 
                                        starting,
                                        p,
                                        reltol=1e-8,
                                        maxiters = 1e4)          
        end
        optimizer = isnothing(optimizer) ? default_optimizer[] : optimizer
        sol = solve(probl,optimizer)
        fill_pars!(p,sol.u)
        residual!(p)
        return (sol,p)
    end
    function check_bounds_on_starting_vector!(starting_vector,box)
        length(box.lb) == length(starting_vector) || resize!(starting_vector,length(box.lb))
        for (i,(lbi,si,ubi)) in enumerate(zip(box.lb,starting_vector,box.ub))
            box_span = ubi - lbi
            !(lbi < si < ubi) || continue     
            starting_vector[i] = si < lbi ?  lbi + 1e-6*box_span : ubi - 1e-6*box_span
        end
    end
    """
    split_peaks(mp::MultiPeaks{N,B,P};
    x::Union{Nothing,AbstractVector}=nothing,
    sort_by::Symbol=:μ,rev::Bool=true) where {N,B,P}

Splits each peak to a separated vector and puts all this peaks to the XY...Y matrix with 
X - independent coordinate, Y...Y each column for the peak intensity  
"""
function split_peaks(mp::MultiPeaks{N,B,P};
    x::Union{Nothing,AbstractVector}=nothing,
    sort_by::Symbol=:μ,rev::Bool=false) where {N,B,P} 
    x_calc = isnothing(x) ? mp.x : x
    o = Matrix{Float64}(undef,length(x_calc),N)
        if sort_by ∈ names(mp.peaks[1])
            inds = Vector{Int}(undef,N)
            sortperm!(inds,[getproperty(p,sort_by) for p in mp.peaks ], rev=rev)
        else
            inds = collect(1:N)
        end
        @views for i = 1:N 
            j = inds[i]
            v =  o[:,i] 
            map!(mp[j],v,x_calc)
        end
        return o,inds
    end
    """
    statistics(mp::MultiPeaks)

Evaluates basic statistics on peaks fitter structure
"""
function statistics(mp::MultiPeaks)
        N = length(mp.x) - 1 # degrees of freedom
        NmP = N - parnumber(mp) # reg degrees of freedom
        mean_y0 = mean(mp.y0)
        SSR = sum(t->(t-mean_y0)^2,mp.y) # regression sum of squares 
        SSE = sum(t->t^2,mp.r) # sum of square errors
        SST = sum(t->(t-mean_y0)^2,mp.y0) # total sum of squares
        return (V = sqrt(SSE/NmP),
                r=( 1 - SSE/SST),
                radj=(1 - (N/NmP)*(SSE/SST)),
                SSE=SSE,
                SST=SST,
                SSR=SSR,
                freedom_degrees = NmP,
                total_degrees = N )
    end
    """
    covariance_matrix(mp::MultiPeaks{N,B,P}) where {N,B<:AbstractBaseLine,P<:AbstractPeak}

Function to estimate the covariance as a pseudo-inverse of approximate hessian matrix:

Cov = Hₐ⁻¹*SSE/(N-P) where Hₐ = (J'*J), J is the Jacobian

"""
function covariance_matrix(mp::MultiPeaks{N,B,P}) where {N,B<:AbstractBaseLine,P<:AbstractPeak}
        params = fill_vector_with_pars!(Float64[],mp,resizable=true)
        npoints = length(mp.x)
        npars = parnumber(mp)
        if is_has_gradient(P)
            grad_fill! = (g,x) -> fill_from_tuples!(g,∇(x,1.0,mp)) 
        else # if gradient is not provided 
            #   TBW
 #           grad_fill! = (g,x) -> copyto!(g, ForwardDiff.gradient())
        end
        J = Matrix{Float64}(undef,npoints,npars)
        H = Matrix{Float64}(undef,npars,npars)
        for (i,xi) in enumerate(mp.x)
            ji = @view J[i,:]
            grad_fill!(ji,xi)
        end
        # calculating the approximate hessian J'*J
        mul!(H,transpose(J),J)
        stat_struct = statistics(mp)
        cov_mat = pinv(H)*stat_struct.SSE/stat_struct.freedom_degrees
        return (covariance = cov_mat,variance = diag(cov_mat),jacobian=J,hessian=H)
    end
    """
    box_constraints(p::MultiPeaks{N}) where N

Evaluates box constraints on the peaks parameters 
"""
function box_constraints(p::MultiPeaks{N,B,P};
            constraints_expansion::Float64=1.0,
            allow_negative_peaks_amplitude::Bool=true,
            allow_negative_baseline_tangent::Bool=true) where {N,B,P<:Union{GaussPeak,
                                                                LorentzPeak,VoigtPeak,WeibullPeak}}
        # method calculates box constraints 
        # of the feasible region (dumb version)
        npars = parnumber(p) # total parameters number in all peaks
        base_npars = parnumber(B)
        peaks_npars = npars-base_npars # peaks parameters number 

        lb = Vector{Float64}(undef,npars)
        ub = Vector{Float64}(undef,npars)

        mu_bnds = extrema(p.x)
        s_bnds = (0.0,p.x[end])
        a_bnds = extrema(p.y0)

        lb[1] = a_bnds[1]
        ub[1] = a_bnds[2]

        allow_negative_peaks_amplitude || a_bnds[1] >=0 || (a_bnds=(0.0,a_bnds[2])) 
        
        max_tan = abs(-(a_bnds...)/-(mu_bnds...))
        allow_negative_baseline_tangent ? fill!(view(lb,2:base_npars),-max_tan) : fill!(view(lb,2:base_npars),0.0)
        fill!(view(ub,2:base_npars),max_tan)

        _lb = @view lb[base_npars+1:npars]
        _ub = @view ub[base_npars+1:npars]
        if P<:VoigtPeak
            filler_fun_lower = _->(mu_bnds[1],s_bnds[1],a_bnds[1],0.0)
            filler_fun_upper = _->(mu_bnds[2],s_bnds[2],a_bnds[2],1.0)
        elseif P<:WeibullPeak
            filler_fun_lower = _->(mu_bnds[1],s_bnds[1],a_bnds[1],0.0)
            filler_fun_upper = _->(mu_bnds[2],s_bnds[2],a_bnds[2],10.0)
        else
            filler_fun_lower = _->(mu_bnds[1],s_bnds[1],a_bnds[1])
            filler_fun_upper = _->(mu_bnds[2],s_bnds[2],a_bnds[2])
        end
        fill_from_tuples!(_lb,ntuple(filler_fun_lower, Val(N)))
        fill_from_tuples!(_ub,ntuple(filler_fun_upper, Val(N)))
        return (lb=constraints_expansion*lb,ub=constraints_expansion*ub)
    end
    function distance_matrix(p_vect::AbstractVector)
        total_peaks_number = length(p_vect)#sum(peaknumber,mp)
        distance_mat = Matrix{Float64}(undef,total_peaks_number,total_peaks_number)
        for j in 1:total_peaks_number  # Iterate columns
            for i in 1:j # Iterate rows up to the current column
                distance_mat[i, j] = peaks_distance(p_vect[i],p_vect[j])[1]
            end
        end
        return Symmetric(distance_mat)
    end
    function collect_peaks(multipeaks_collection)
        for mp in multipeaks_collection
            
        end
    end
    #plot recipe
   include("PlotRecipes.jl")

end