using Eirene # Persistent Homology
using DelimitedFiles # Read CSV Files
using Plots # Visualize Datasets
using PersistenceDiagrams # Persistence Diagrams
using LinearAlgebra # Matrix Algebra
using Colors # Use Various Colors
using Graphs # Create Edges Between Vertices
using CSV # Save DataFrames to CSV Files
using DataFrames # Create and Use Dataframes

eul = 0.5772156649
subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))

function compare_rows(row1, row2)
    if row1[2] < row2[2]  # Compare by second column
        return true
    elseif row1[2] > row2[2] # Compare by second column
        return false
    else  # If values are equal, compare by first column in reverse order
        return row1[1] > row2[1]
    end
end

function barcode_dim(p::Dict{String, Any}, dim::Int64)
    # Compute barcode from PersHom
    barcode = Eirene.barcode(p, dim=dim)

    # Sort barcode by birth
    sorted_indices = sortperm(barcode[:, 2])
    sorted_barcode = barcode[sorted_indices, :]
    
    # Compute class representatives
    class_reps = [Eirene.classrep(p, dim=dim, class=i) for i in sorted_indices]
    
    # Compute persistence values
    persistences = -(sorted_barcode[:, 2], sorted_barcode[:, 1])
    
    return sorted_barcode, class_reps, persistences
end

function find_all(f, a) # Faster implementation of findall
    return [i for (i, x) in enumerate(a) if f(x)]
end

function any_equal_col(a::AbstractMatrix, b::AbstractMatrix) # Check if two matrices have at least one identical column
    return any(x -> in(x, eachcol(b)), eachcol(a))
end

function adjacency(a, b::AbstractMatrix) # Check for k-nearness
    return any(x -> any_equal_col(x, b), a)
end

function find_mergable_cycles(b, cycles) # Find first order merge clusters
    cluster = [[] for i in 1:length(cycles)]
    for i in 2:length(cycles)
        range = setdiff(1:i-1, union(cluster...))
        for j in range
            if any_equal_col(cycles[i], cycles[j]) == 1 && b[i, 1] < b[j, 2] 
                cluster[i] = vcat(cluster[i], [j])
            elseif findfirst(x -> any_equal_col(cycles[i], cycles[x]) == 1, cluster[j]) != nothing && b[i, 1] < b[j, 2] 
                cluster[i] = vcat(cluster[i], [j])
            end
        end
    end
    return cluster
end

function extract_second_column_values(b, merged_cycles) # Compute thresholds where merging instances happen
    return [b[cycle, 2] for cycle in merged_cycles]
end

function merging_process(b::Matrix{Float64}, c::Vector{Matrix{Int64}})
    # Assign number of homology classes
    s = size(b, 1)

    # Find first order merge clusters
    mergable_cycles = find_mergable_cycles(b, c)

    # Find merge thresholds
    t = extract_second_column_values(b, mergable_cycles)

    return mergable_cycles, t
end

function plot_centrality1(num, dom, time, bar, pers, mc, const1, const2, lambda)
    if isempty(mc[num])
        return [((bar[num, 1] < d <= bar[num, 2]) ? const1 * (d - bar[num, 1]) : (d > bar[num, 2]) ? const1 * pers[num] : 0) for d in dom]
    elseif length(mc[num]) == 1
        return [((bar[num, 1] < d < time[num][1]) ? const1 * (d - bar[num, 1]) : (time[num][1] <= d < bar[num, 2]) ? const1 * (d - bar[num, 1]) + const2 * sum(pers[mc[num]]) : (d >= bar[num, 2]) ? const1 * pers[num] + const2 * sum(pers[mc[num]]) : 0) for d in dom]
    else
        V = zeros(length(dom))
        for i in eachindex(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bar[num, 1])
            elseif (Y = findlast(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * sum(pers[mc[num]][1:Y])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * sum(pers[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1 * pers[num] + const2 * sum(pers[mc[num]])
                break
            end
        end
    end
    return V
end

function plot_centrality2(num, dom, time, bars, pers, class, const1, const2, lambda)
    if isempty(class[num])
        return [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d - bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
    elseif length(class[num]) == 1
        return [(bars[num, 1] < d < time[num][1]) ? const1 * (d - bars[num, 1]) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (d - bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bars[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bars[num, 1])
            elseif (Y = findlast(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bars[num, 1]) + dot(time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bars[num, 2]
                V[i] = const1 * (dom[i] - bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]])
            elseif dom[i] >= bars[num, 2]
                V[i:end] .= const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]])
                break
            end
        end
    end
    return V
end

function plot_centrality3(num, dom, time, bars, pers, class, const1, const2, p, lambda)
    if isempty(class[num])
        return [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d - bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
    elseif length(class[num]) == 1
        return [(bars[num, 1] < d < time[num][1]) ? const1 * (d - bars[num, 1]) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (d - bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bars[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bars[num, 1])
            elseif (Y = findlast(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bars[num, 1]) + dot(1 .- time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bars[num, 2]
                V[i] = const1 * (dom[i] - bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
            elseif dom[i] >= bars[num, 2]
                V[i:end] .= const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                break
            end
        end
    end
    return V
end

function transfer1(class, bar, pers)
    stack = [isempty(class[i]) ? 0 : sum(pers[class[i]]) for i in 1:length(class)]

    for i in 1:length(class)
        stack[i] = stack[i] + sum(stack[class[i]])
    end

    return stack
end

function plot_centrality4(num, dom, time, bar, pers, class, const1, const2)
    stack = transfer1(class, bar, pers)
    const1pers = const1 * pers[num]
    
    if isempty(class[num])
        return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
    end
    
    if length(class[num]) == 1
        return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d - bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d - bar[num, 1]) + const2 * stack[num]
            elseif d >= bar[num, 2]
                const1pers + const2 * stack[num]
            else
                0
            end
            for d in dom
        ]
    end

    V = zeros(length(dom))
    
    for i in 1:length(dom)
        if bar[num, 1] < dom[i] < time[num][1]
            V[i] = const1 * (dom[i] - bar[num, 1])
        elseif (local Y = findlast(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * (sum(pers[class[num][1:Y]]) + sum(stack[class[num][1:Y]]))
        elseif last(time[num]) <= dom[i] < bar[num, 2]
            V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * stack[num]
        elseif dom[i] >= bar[num, 2]
            V[i:end] .= const1pers + const2 * stack[num]
            break
        end           
    end

    return V
end

function transfer2(class, bar, pers, time)
    stack = [isempty(class[i]) ? 0 : dot(time[i]/bar[i, 2], pers[class[i]]) for i in 1:length(class)]

    for i in 1:length(class)
        stack[i] += sum(stack[class[i]])
    end

    return stack
end

function plot_centrality5(num, dom, time, bar, pers, class, const1, const2)
    stack = transfer2(class, bar, pers, time)
    const1pers = const1 * pers[num]
    
    if isempty(class[num])
        return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
    end
    
    if length(class[num]) == 1
        return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d - bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d - bar[num, 1]) + stack[num]
            elseif d >= bar[num, 2]
                const1pers + stack[num]
            else
                0
            end
            for d in dom
        ]
    end

    V = zeros(length(dom))

    for i in 1:length(dom)
        if bar[num, 1] < dom[i] < time[num][1]
            V[i] = const1 * (dom[i] - bar[num, 1])
        elseif (local Y = findlast(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            V[i] = const1 * (dom[i] - bar[num, 1]) + dot(time[num][1:Y] / bar[num, 2], pers[class[num][1:Y]]) + sum(stack[class[num][1:Y]])
        elseif last(time[num]) <= dom[i] < bar[num, 2]
            V[i] = const1 * (dom[i] - bar[num, 1]) + stack[num]
        elseif dom[i] >= bar[num, 2]
            V[i:end] .= const1pers + stack[num]
            break
        end       
    end

    return V
end

function transfer3(class, bar, pers, time)
    stack = [isempty(class[i]) ? 0 : dot(1 .- (time[i] / bar[i, 2]), pers[class[i]]) for i in 1:length(class)]
    
    for i in 1:length(class)
        stack[i] += sum(stack[class[i]])
    end
    
    return stack
end

function plot_centrality6(num, dom, time, bar, pers, class, const1, const2)
    stack = transfer3(class, bar, pers, time)
    const1pers = const1 * pers[num]
    
    if isempty(class[num])
        return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]    
    end
    
    if length(class[num]) == 1
        return [
                if (bar[num, 1] < d < time[num][1])
                    const1 * (d - bar[num, 1])
                elseif (time[num][1] <= d < bar[num, 2])
                    const1 * (d - bar[num, 1]) + stack[num]
                elseif d >= bar[num, 2]
                    const1pers + stack[num]
                else
                    0
                end
                for d in dom
            ]       
    end

    V = zeros(length(dom))

    for i in 1:length(dom)
        if bar[num, 1] < dom[i] < time[num][1]
            V[i] = const1 * (dom[i] - bar[num, 1])
        elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            V[i] = const1 * (dom[i] - bar[num, 1]) + dot(1 .- time[num][1:Y] / bar[num, 2], pers[class[num][1:Y]]) + sum(stack[class[num]][1:Y])
        elseif last(time[num]) <= dom[i] < bar[num, 2]
            V[i] = const1 * (dom[i] - bar[num, 1]) + stack[num]
        elseif dom[i] >= bar[num, 2]
            V[i:end] .= const1pers + stack[num]
            break
        end
    end

    return V
end

function centrality_plot(bars, pers, time, class, color_vec)
    x = collect(range(0, stop = maximum(bars[:, 2]), length = 2000))
    centrality_plots = []

    for i in 4:6
        centrality_measure = eval(Meta.parse("plot_centrality$i"))
        push!(centrality_plots, plot(x, hcat([centrality_measure(j, x, time, bars, pers, class, 1, 1) for j in 1:length(pers)]...), legend = false, palette = color_vec, ylabel = "J"*subscript(i)*"(σ,ϵ)", xlabel = "Threshold ϵ"))
    end

    return plot(centrality_plots..., layout = (1,3))
end

function load_data(data_path)
    data = transpose(readdlm(data_path, ',', Float64, '\n', header=true)[1])
    return data[1, :], data[2, :]
end

function save_scatter_plot(x, y, save_path)
    peach_fuzz_color = RGB(1.0, 0.75, 0.60)
    mintgreen_color = RGB(0.7333333333333333, 0.3803921568627451, 0.6784313725490196)
    plot_data = scatter(x, y, legend=false, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, grid=false, background_color=:transparent, foreground_color=:black, xlabel = "x", ylabel = "y")
    savefig(plot_data, save_path)
end

function process_data(data, max_dim=1)
    holes = eirene(data, model="pc", maxdim=max_dim)
    bars, reps, pers = barcode_dim(holes, max_dim)
    pers_mult = bars[:, 2] ./ bars[:, 1]
    el = sum(log.(log.(pers_mult)))/length(pers_mult)
    pers_trans = log.(log.(pers_mult)) .- eul .- el
    
    class, time = merging_process(bars, reps)
    diag = PersistenceDiagram([(bars[i, 1], bars[i, 2]) for i in 1:size(bars, 1)])
    
    return bars, reps, pers, pers_mult, pers_trans, class, time, diag
end

function save_diagram_and_measures(bars, pers, time, class, save_diagram_path, save_measures_path, diag, color_vec)
    cent_plots = centrality_plot(bars, pers, time, class, color_vec)
    savefig(cent_plots, save_measures_path)
end

function bootstrap_data(x, y, num_points)
    indices = rand(1:length(x), num_points)
    return x[indices], y[indices]
end