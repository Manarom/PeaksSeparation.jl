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
    (out_mat,inds) = split_peaks(m)
    # N=peak_number(m)
    for (i,c)in enumerate(eachcol(out_mat))
        @series begin 
            label:="p_$(inds[i])"    
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
@recipe function f(x,m::AbstractCurveComponent)
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
    linewidth-->2
    fillrange-->0
    fillalpha-->0.3
    markershape:=:none
    return (x,m.(x))
end 