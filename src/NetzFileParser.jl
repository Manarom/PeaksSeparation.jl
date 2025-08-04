module NetzFileParser
    using CSV,DataFrames
    mutable struct NetzFile
        file_name::AbstractString
        unparsed_headers::Vector{String}
        headers ::Union{Nothing,AbstractDict}
        column_headers ::Vector{String}
        data ::Union{Nothing,AbstractDataFrame}
        data_start::Int
        NetzFile(file::AbstractString;delim=';',kwargs...) = begin 
            nf = new(file,String[],nothing,String[],nothing,1)
            parse_headers!(nf)
            fill_data!(nf;delim=delim,kwargs...)
            return nf
        end
    end
    function parse_headers!(np::NetzFile)   
        @assert isfile(np.file_name) "Not a file"
        !isnothing(np.headers) || (np.headers = Dict{String,String}())
        open(np.file_name) do io 
            run_flag = true
            np.data_start == 1 || begin 
                        np.headers = Dict{String,String}()  
                        np.unparsed_headers = Vector{String}()
                        np.data_start = 1
                    end
            while !eof(io)
                np.data_start +=1
                ln = readline(io)
                push!(np.unparsed_headers,ln)
                length(ln)>0 || continue
                run_flag = occursin("#",ln)
                !occursin("##",ln) ||begin 
                                        ln = strip(ln,'#')
                                        np.column_headers = String.(split(ln,";"))
                                    end
                if !run_flag 
                    resize!(np.unparsed_headers, length(np.unparsed_headers) - 1)
                    return nothing
                end
                if occursin(":",ln)
                    ln = strip(ln,'#')
                    str_spl = split(ln,":")
                    np.headers[str_spl[1]] = length(str_spl) > 2 ? join(str_spl[2:end],":") : str_spl[2]
                end
            end
        end
    end
    function write_like(nf::NetzFile,new_file; 
        data::Union{Nothing,AbstractArray,AbstractDataFrame}=nothing,
        delim=";",kwargs...)
        open(new_file,"w") do io 
            for l in nf.unparsed_headers
                println(io,l)
            end
        end
        if isnothing(data)
            CSV.write(new_file,nf.data;append=true,delim=delim,kwargs...)
        elseif isa(data,AbstractDataFrame)
            CSV.write(new_file,data;append=true,delim=delim,kwargs...)
        elseif isa(data,AbstractArray)
            CSV.write(new_file,DataFrame(data,:auto);append=true,delim=delim,kwargs...)
        end
    end
    function fill_data!(nf;kwargs...)
        nf.data = CSV.read(nf.file_name,DataFrame;
            ignoreemptyrows=true,
            skipto = nf.data_start
            ,kwargs...) #,types = Float64
        try 
            rename!(nf.data,nf.column_headers)
        catch ex 
            @show ex
        end
    end
end