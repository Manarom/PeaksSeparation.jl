
macro iterator_unpack(N::Int, x)
    expr = Expr(:tuple)
    push!(expr.args, quote
        (val, state) = iterate($(esc(x)))
        val
    end)
    for i = 2:N
        push!(expr.args, quote
            (val, state) = iterate($(esc(x)), state)
            val
        end)
    end
    expr
end
function fill_from_tuples!(v,  args;resizable::Bool=false)
    counter = 0
    for arg in args
        for ei in arg
            !isnothing(ei) || continue
            counter+=1
            if !resizable
                v[counter] = ei
            else
                set_or_push!(v,counter,ei)
            end
        end
    end
    return v
end
function add_from_tuples!(v,  args ;resizable::Bool=false)
    counter = 0
    for arg in args
        for ei in arg
            !isnothing(ei) || continue
            counter+=1
            if !resizable
                v[counter] += ei
            else
                add_or_push!(v,counter,ei)
            end
        end
    end
    return v
end
function add_from_tuples!(to_add_to,args,bynary_operator::Function;resizable::Bool=false)
    counter = 0
    for arg in args
        for ei in arg
            !isnothing(ei) || continue
            counter+=1
            if !resizable
                v[counter] =bynary_operator(v[counter],ei)
            else
                add_or_push!(v,counter,ei,bynary_operator)
            end
        end
    end
    return v
end
function set_or_push!(v, i::Int, val)
    M=length(v)
    if 1 <= i <= M
        v[i] = val
    elseif i == M + 1
        push!(v, val)
    else i > M + 1
        resize!(v, i) 
        v[i] = val
    end
end
function add_or_push!(v, i::Int, val)
    M=length(v)
    if 1 <= i <= M
        v[i] += val
    elseif i == M + 1
        push!(v, val)
    else i > M + 1
        resize!(v, i) 
        v[i] = val
    end
end
function add_or_push!(v, i::Int, val,bynary_operator::Function)
    M=length(v)
    if 1 <= i <= M
        v[i] = bynary_operator(v[i],val)
    elseif i == M + 1
        push!(v, val)
    else i > M + 1
        resize!(v, i) 
        v[i] = val
    end
end

tuple_mult(tpl,factor) = tuple((x * factor for x in tpl)...)
function replace_nans!(x,new_value)
    for i in eachindex(x)
        isnan(x[i]) || continue
        x[i] = new_value
    end
end
