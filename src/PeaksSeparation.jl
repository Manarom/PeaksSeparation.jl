module PeaksSeparation
    using Optimization, Peaks, StaticArrays,OptimizationOptimJL,ForwardDiff,RecipesBase
    export MultiPeaks,
            fit_peaks,
            fit_peaks!,
            split_peaks

    const default_optimizer = Ref(ParticleSwarm(n_particles=100))
    const fwhm = 2*sqrt(2*log(2.0)) # contant to evaluate the dispersion of Gaussian from peaks's fwhm
    include("Utils.jl")
    f_gauss(x,μ,σ,a) = a*exp(- ^(x-μ,2) * 0.5 * ^(σ,-2))
     
    sqr(x) = x^2
    function ∇gauss(x,μ,σ,a) 
        fval =  f_gauss(x,μ,σ,a)
        fa = isapprox(a,0.0,atol=1e-6) ? f_gauss(x,μ,σ,1) : fval/a  
        (-(x - μ)*fval/(σ^2), fval*( ^(x-μ,2)/^(σ,3)) ,fa )
    end 
    

    #∇mu_gauss(x,μ,σ,a) =  -(x - μ)*f_gauss(μ,σ,a,x)/σ^2
    #∇a_gauss(μ,σ,a,x) = f_gauss(μ,σ,1.0,x)
    #∇σ_gauss(μ,σ,a,x) = f_gauss(μ,σ,a,x) *( ^(x-μ,2)/^(σ,3))
    #∇mu_gauss!(g,μ,σ,a,x,y) = @. g= -(x - μ)*y/σ^2
    #∇σ_gauss(g,μ,σ,a,x,y) = @. g = ^(x-μ,2)* y /^(σ,3)



    # N - is the number of parameters of this component
    abstract type AbstractCurveComponent{N} end
    eval_to!(y,cc::AbstractCurveComponent,x) = map!(cc,y,x)
    """
    add_to!(y,cc::AbstractCurveComponent,x)

 Adds result of curve component (like peak or baseline) evaluation to the vector y
"""
add_to!(y,cc::AbstractCurveComponent,x) = @. y += cc(x)
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
        !(field == :b || field == :names || field == :flags) || return getfield(cc,field)
        nm = names(cc)
        for i in 1:N
            field != nm[i] || return parameters(cc)[i]
        end
        getfield(cc,field)
    end
    # baseline supertype
    abstract  type AbstractBaseLine{N}<:AbstractCurveComponent{N} end
    abstract type PolynomialBaseLine{N}<:AbstractBaseLine{N} end
    function ∇(x,::PolynomialBaseLine{N}) where N
        ntuple(i->^(x,i-1),Val(N))
    end
    """
    fill_pars!(::PolynomialBaseLine,x)

function to fill baseline polynomial coefficients from the iterable object x

"""
function fill_pars!(::PolynomialBaseLine,x) end
    for i in 1:5 # reserving function for up to five 
        @eval function fill_pars!(bl::PolynomialBaseLine{$i},x) 
                setfield!(bl,:parameters , @iterator_unpack($i,x))
                return bl
        end
    end
    (bl::PolynomialBaseLine)(x) = @evalpoly(x,getfield(bl,:parameters)...)
    mutable struct LinearBaseLine<:PolynomialBaseLine{2}
        parameters::NTuple{2,Float64}
        names::NTuple{2,Symbol}
        flags::NTuple{2,Bool}
        LinearBaseLine() = begin 
            new((0.0,0.0),(:b1,:b2),(false,false))
        end
    end
       # (bl::LinearBaseLine)(x) = bl.parameters[1] + bl.parameters[2]*x

    abstract type AbstractPeak{N}<:AbstractCurveComponent{N} end
    mutable struct GaussPeak<:AbstractPeak{3}
        parameters::NTuple{3,Float64}
        names::NTuple{3,Symbol}
        flags::NTuple{3,Bool}
        index::Int
        GaussPeak(i::Int=1,mu=0.0,sigma=0.0,a=0.0) = begin 
            new((mu,sigma,a),(:μ,:σ,:a),(false,false,false),i)
        end
    end
    """
    (gp::GaussPeak)(x)

When peak object is evaluated on variable it returns the value of Gaussian
"""
    (gp::GaussPeak)(x) = f_gauss(x,getfield(gp,:parameters)...)
    """
    ∇(x,p::GaussPeak)

Calculates the gradient of GaussPeak as a Tuple of (dμ,dσ,da) for 
a single value of x
"""
    ∇(x,p::GaussPeak) = ∇gauss(x,getfield(p,:parameters)...)
    function fill_pars!(p::GaussPeak,x)
            setfield!(p,:parameters, @iterator_unpack(3,x))
            return p
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

 #peak functions   
    mutable struct MultiPeaks{N,B<:AbstractBaseLine,P<:AbstractPeak,T<:AbstractVector}
        x::T
        y0::T
        y::T
        r::T
        peaks::SVector{N,P}
        flags
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
                            new{N,B,P,typeof(x)}(x,y0,
                                        similar(y0),similar(y0),
                                        SVector(ntuple(P,Val(N))),
                                        (false,false,false),
                                        B())
        end
    end
function Base.getindex(mp::MultiPeaks{N},i) where N
    (i >N || i < 0) ? error("Index out of bounds") : return i == 0 ? mp.baseline : mp.peaks[i]
end
parnumber(::MultiPeaks{N,B,P}) where {N,B,P} = parnumber(B)+ N*parnumber(P) # total parnumber of the curve 
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

Evaluates calculated values and fits y 
"""
function fill_y!(mp::MultiPeaks) 
        fill!(mp.y,0.0)
        add_baseline!(mp)
        add_peaks!(mp)
        return mp
    end
    """
    residual!(mp::MultiPeaks{N}) where N

Evaluates residual vector
"""
function residual!(mp::MultiPeaks{N}) where N
        fill_y!(mp)
        @. mp.r =mp.y - mp.y0
        return mp
    end
function grad!(g,x,mp::MultiPeaks{N}) where N

end
    """
    fill_pars!(mp::MultiPeaks{N},x) where N

Fills all curve components parameters from iterable object x
"""
function fill_pars!(mp::MultiPeaks{N},x) where N
        itr = Iterators.Stateful(x)
        fill_pars!(mp.baseline,itr)
        fill_pars!(mp.peaks,itr)
        return nothing
    end
    """
    discr(x,mp::MultiPeaks{N}) where N

Evaluates discrepancy value
"""
function discr(x,mp::MultiPeaks{N}) where N
        fill_pars!(mp,x)
        residual!(mp)
        return sum(sqr,mp.r)
    end
    discr_fun(x,p) = x->discr(x,p)
    #=function grad!(g::AbstractVector,t::Float64 ,e::MultiPeaks{N}) where N
        ∇!(t,e)
        residual!(e,t)
        if t!=e.Tgrad[]
            g[end]= - dot(e.ri,e.∇I) # filling gradient vector
            e.Tgrad[] = t
        end
        return nothing
    end=#
    """
    fill_from_tuples!(v,  args::Vararg;resizable::Bool=false)

Puts all content of args (tuple of tuple or vector of vectors et.c) into a single vector
"""

    peak_number(::MultiPeaks{N}) where N = N
    function fill_vector_with_pars!(v,p::MultiPeaks{N};resizable::Bool=false) where N 
        fill_from_tuples!(v,ntuple(i->parameters(p[i-1]), Val(N+1))...,resizable=resizable)
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
    param_estimate(::Type{GaussPeak},x,y,indices::Vector{Int})

Finds peaks and estimates 
"""
function param_estimate(::Type{GaussPeak},x,y,indices::Vector{Int})
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
    param_estimate(T::Type{GaussPeak},x,y,N::Union{Int,Nothing}) = param_estimate(T,x,y,maxima_indices(y,N))
    
    function fill_starting_vector!(starting_vector,mp::MultiPeaks{PeakNum,BaseType,PeakType}) where {PeakNum,BaseType,PeakType}
        is_resizable = length(starting_vector) != parnumber(mp)
        base_param_number = parnumber(BaseType)

        return fill_from_tuples!(starting_vector,  
                    (ntuple(t->0,base_param_number),
                    param_estimate(PeakType,mp.x,mp.y0,PeakNum)...)...
                    ;resizable=is_resizable)
    end
    int_floor = Int ∘ floor
    function fit_peaks!(p::MultiPeaks;kwargs...)
        starting_vector =  fill_vector_with_pars!(Float64[],p,resizable=true)
        return fit_peaks(       nothing,
                                nothing;
                                p=p,
                                starting_vector = starting_vector,optimizer = NelderMead(),
                                use_constraints =false,
                                kwargs...)
    end
    function fit_peaks(x::T,y::T ;
                        N::Union{Int,Nothing}=nothing,
                        starting_vector::Union{Nothing,AbstractVector} = nothing,
                        optimizer = nothing,
                        p::Union{Nothing,MultiPeaks} = nothing,
                        use_constraints::Bool=true,
                        constraints_expansion::Float64=1.0,
                        allow_negative_peaks_amplitude::Bool = true,
                        allow_negative_baseline_tangent::Bool = true) where T<:Union{Nothing,AbstractVector}
        
        is_p_provided = !isnothing(p)
        !is_p_provided || (N=peak_number(p))
        #is_fixed_peaks_number = !isnothing(N)
        
        if isnothing(starting_vector)
            #pks = fix_peaks_number(y,N)
            #N = length(pks.indices)
            if !is_p_provided
                p = MultiPeaks(x,y,N)
            end
            starting = fill_starting_vector!(Float64[],p)
        else
            starting = starting_vector
        end
        fill_pars!(p,starting)
        optim_fun = OptimizationFunction(discr,AutoFiniteDiff())
        
        if use_constraints
            box = box_constraints(p,
                            constraints_expansion=constraints_expansion,
                            allow_negative_peaks_amplitude=allow_negative_peaks_amplitude,
                            allow_negative_baseline_tangent = allow_negative_baseline_tangent)

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
        for i = 1:N 
            j = inds[i]
            v = @view o[:,i] 
            map!(mp[j],v,x_calc)
        end
        return o
    end
    """
    box_constraints(p::MultiPeaks{N}) where N

Evaluates box constraints on the 
"""
function box_constraints(p::MultiPeaks{N,B,P};
            constraints_expansion::Float64=1.0,
            allow_negative_peaks_amplitude::Bool=true,
            allow_negative_baseline_tangent::Bool=true) where {N,B,P<:GaussPeak}
        # method calculates box constraints 
        # of the feasible region (dumb version)
        mu_bnds = extrema(p.x)
        s_bnds = (0.0,p.x[end])
        a_bnds = extrema(p.y0)
        allow_negative_peaks_amplitude || a_bnds[1] >=0 || (a_bnds=(0.0,a_bnds[2])) 
        npars = parnumber(p) # total parameters number

        lb = Vector{Float64}(undef,npars)
        ub = Vector{Float64}(undef,npars)

        base_npars = parnumber(B)
        peaks_npars = npars-base_npars # peaks parameters number 
        
        lb[1] = a_bnds[1]
        ub[1] = a_bnds[2]

        max_tan = abs(-(a_bnds...)/-(mu_bnds...))
        allow_negative_baseline_tangent ? fill!(view(lb,2:base_npars),-max_tan) : fill!(view(lb,2:base_npars),0.0)
        fill!(view(ub,2:base_npars),max_tan)

        _lb = @view lb[base_npars+1:npars]
        _ub = @view ub[base_npars+1:npars]

        fill_from_tuples!(_lb,ntuple(i->(mu_bnds[1],s_bnds[1],a_bnds[1]), Val(N))...)
        fill_from_tuples!(_ub,ntuple(i->(mu_bnds[2],s_bnds[2],a_bnds[2]), Val(N))...)

        return (lb=constraints_expansion*lb,ub=constraints_expansion*ub)
    end
    #plot recipe
    @recipe function f(m::MultiPeaks)
        minorgrid--> true
        gridlinewidth-->2
        dpi-->600
        xlabel-->"Temperature"
        ylabel-->"DSC disgnal"
        label-->"y0"
        linewidth-->3
        markershape --> :diamond
        markersize -->3
        markeralpha-->0.5
        out_mat = split_peaks(m)
        # N=peak_number(m)
        for (i,c)in enumerate(eachcol(out_mat))
            @series begin 
                label:="p_$(i)"    
                linewidth:=2
                fillrange:=0
                fillalpha:=0.3
                markershape:=:none
                (m.x, c)
            end
        end
        @series begin 
            label:= "Σp_i"    
            linewidth:=2
            markershape:=:none
            (m.x, m.y)
        end
        return (m.x,m.y0)
    end 
    function replace_nans!(x,new_value)
        for i in eachindex(x)
            isnan(x[i]) || continue
            x[i] = new_value
        end
    end
end