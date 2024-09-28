using Eirene # Persistent Homology
using DelimitedFiles # Read CSV Files
using Plots
using PersistenceDiagrams
using LinearAlgebra
using Colors
using Graphs
using CSV
using DataFrames

eul = 0.5772156649
subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))

function compare_rows(row1, row2)
    if row1[2] < row2[2]  # Compare by second column (y)
        return true
    elseif row1[2] > row2[2]
        return false
    else  # If y values are equal, compare by first column (x) in reverse order
        return row1[1] > row2[1]
    end
end

function barcode_dim(p::Dict{String, Any}, dim::Int64)
    # Compute barcode from PersHom
    barcode = Eirene.barcode(p, dim=dim)

    # Sort barcode by birth
    sorted_indices = sortperm(barcode[:, 2]) # by= x -> barcode[x, 2]
    sorted_barcode = barcode[sorted_indices, :] #barcode[sorted_indices, 1:2]
    
    # Compute class representatives
    class_reps = [Eirene.classrep(p, dim=dim, class=i) for i in sorted_indices]
    
    # Compute persistences
    persistences = -(sorted_barcode[:, 2], sorted_barcode[:, 1])
    
    return sorted_barcode, class_reps, persistences
end

function find_all(f, a)
    return [i for (i, x) in enumerate(a) if f(x)]
end

function any_equal_col(a::AbstractMatrix, b::AbstractMatrix)
    return any(x -> in(x, eachcol(b)), eachcol(a))
end

function adjacency(a, b::AbstractMatrix)
    return any(x -> any_equal_col(x, b), a)
end

#function find_equal_columns(cycles)
#    return vcat([[]]. [find_all(x -> any_equal_col(cycles[i], cycles[x]) == 1, i-1:length(cycles)-1) for i in 2:length(cycles)])
#end

function find_mergable_cycles(b, cycles)
    cluster = [[] for i in 1:length(cycles)]
    for i in 2:length(cycles)
        range = setdiff(1:i-1, union(cluster...))
        for j in range
            if any_equal_col(cycles[i], cycles[j]) == 1 # && b[:,2][i] > b[:,2][j] > b[:,1][i]
                cluster[i] = vcat(cluster[i], [j])
            elseif findfirst(x -> any_equal_col(cycles[i], cycles[x]) == 1, cluster[j]) != nothing
                cluster[i] = vcat(cluster[i], [j])
            end
        end
    end
    return cluster
end

# function merge_cycles!(merged_cycles, cycles, adjacency_func, b)
#     index = find_all(x -> x != [], merged_cycles)
    
#     changed = [[]]

#     for i in index
#         range = 2:length(merged_cycles[i])
#         range = setdiff(range, union(cluster...))
#         if (local p = range[find_all(x -> any_equal_col(cycles[i], cycles[x]) == 1 && b[:,2][i] > b[:,2][x] > b[:,1][i], range)]) != nothing
#             merged_cycles = vcat(merged_cycles[], p)
#         end
#         possible = filter(x -> x ∉ merged_cycles[i], 1:length(merged_cycles))
        
#         while (local n = possible[find_all(x -> adjacency_func(cycles[merged_cycles[i]], cycles[x]) == 1 && b[:,2][i] > b[:,2][x] > b[:,1][i], possible)]) != []
#             union!(merged_cycles[i], n)
#             possible = filter(x -> x ∉ merged_cycles[i], 1:length(merged_cycles))
# 	end
#     end
# end

# function remove_intersecting_cycles!(merged_cycles)
#     for i in 1:(length(merged_cycles) - 1)
#         [filter!((x) -> x ∉ intersect(merged_cycles[i], merged_cycles[j]), merged_cycles[j]) for j in i+1:length(merged_cycles)]
#     end
# end

# function sort_cycles!(merged_cycles, b)
#     [sort!(cycle, by = x -> b[x,2]) for cycle in merged_cycles]
# end

function extract_second_column_values(b, merged_cycles)
    return [b[cycle, 2] for cycle in merged_cycles]
end

function merging_process(b::Matrix{Float64}, c::Vector{Matrix{Int64}})
    s = size(b, 1)

#    equal_columns = find_equal_columns(c)
    mergable_cycles = find_mergable_cycles(b, c)

    # merge_cycles!(mergable_cycles, c, (x, y) -> adjacency(x, y), b)
    # remove_intersecting_cycles!(mergable_cycles)
    # sort_cycles!(mergable_cycles, b)

    t = extract_second_column_values(b, mergable_cycles)

    return mergable_cycles, t
end

function plot_centrality1(num, dom, time, bar, pers, mc, const1, const2, p, lambda)
    if isempty(mc[num])
        if p == 1
            return [((bar[num, 1] < d <= bar[num, 2]) ? const1 * (d - bar[num, 1]) : (d > bar[num, 2]) ? const1 * pers[num] : 0) for d in dom]
        # elseif p == 2
        #     return [((bar[num, 1] < d <= bar[num, 2]) ? const1 * (d / bar[num, 1]) : (d > bar[num, 2]) ? const1 * pers[num] : 0) for d in dom]
        # else
        #     return [((bar[num, 1] < d <= bar[num, 2]) ? const1 * (log(log(d / bar[num, 1])) - eul - lambda) : (d > bar[num, 2]) ? const1 * pers[num] : 0) for d in dom]
        end
    elseif length(mc[num]) == 1
        if p == 1
            return [((bar[num, 1] < d < time[num][1]) ? const1 * (d - bar[num, 1]) : (time[num][1] <= d < bar[num, 2]) ? const1 * (d - bar[num, 1]) + const2 * sum(pers[mc[num]]) : (d >= bar[num, 2]) ? const1 * pers[num] + const2 * sum(pers[mc[num]]) : 0) for d in dom]
        # elseif p == 2
        #     return [((bar[num, 1] < d < time[num][1]) ? const1 * (d / bar[num, 1]) : (time[num][1] <= d < bar[num, 2]) ? const1 * (d / bar[num, 1]) + const2 * sum(pers[mc[num]]) : (d >= bar[num, 2]) ? const1 * pers[num] + const2 * sum(pers[mc[num]]) : 0) for d in dom]
        # else
        #     return [((bar[num, 1] < d < time[num][1]) ? const1 * (log(log(d / bar[num, 1])) - eul - lambda) : (time[num][1] <= d < bar[num, 2]) ? const1 * (log(log(d / bar[num, 1])) - eul - lambda) + const2 * sum(pers[mc[num]]) : (d >= bar[num, 2]) ? const1 * pers[num] + const2 * sum(pers[mc[num]]) : 0) for d in dom]
        end
    else
        V = zeros(length(dom))
        for i in eachindex(dom)
            if p == 1
                if bar[num, 1] < dom[i] < time[num][1]
                    V[i] = const1 * (dom[i] - bar[num, 1])
                elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                    V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * sum(pers[mc[num]][1:Y])
                elseif last(time[num]) <= dom[i] < bar[num, 2]
                    V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * sum(pers[mc[num]])
                elseif dom[i] >= bar[num, 2]
                    V[i:end] .= const1 * pers[num] + const2 * sum(pers[mc[num]])
                    break
                end
            # elseif p == 2
            #     if bar[num, 1] < dom[i] < time[num][1]
            #         V[i] = const1 * (dom[i] / bar[num, 1])
            #     elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            #         V[i] = const1 * (dom[i] / bar[num, 1]) + const2 * sum(pers[mc[num]][1:Y])
            #     elseif last(time[num]) <= dom[i] < bar[num, 2]
            #         V[i] = const1 * (dom[i] / bar[num, 1]) + const2 * sum(pers[mc[num]])
            #     elseif dom[i] >= bar[num, 2]
            #         V[i:end] .= const1 * pers[num] + const2 * sum(pers[mc[num]])
            #         break
            #     end
            # else
            #     if bar[num, 1] < dom[i] < time[num][1]
            #         V[i] = const1 * (log(log(dom[i] / bar[num, 1])) - eul - lambda)
            #     elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            #         V[i] = const1 * (log(log(dom[i] / bar[num, 1])) - eul - lambda) + const2 * sum(pers[mc[num]][1:Y])
            #     elseif last(time[num]) <= dom[i] < bar[num, 2]
            #         V[i] = const1 * (log(log(dom[i] / bar[num, 1])) - eul - lambda) + const2 * sum(pers[mc[num]])
            #     elseif dom[i] >= bar[num, 2]
            #         V[i:end] .= const1 * pers[num] + const2 * sum(pers[mc[num]])
            #         break
            #     end
            end
        end
        return V
    end
end

function plot_centrality2(num, dom, time, bars, pers, class, const1, const2, p, lambda)
    if isempty(class[num])
        if p == 1
            V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d - bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        # elseif p == 2
        #     V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d / bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        # else
        #     V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (log(log(d/ bars[num, 1])) - eul - lambda) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        end
    elseif length(class[num]) == 1
        if p == 1
            V = [(bars[num, 1] < d < time[num][1]) ? const1 * (d - bars[num, 1]) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (d - bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        # elseif p == 2
        #     V = [(bars[num, 1] < d < time[num][1]) ? const1 * (d / bars[num, 1]) :
        #         (time[num][1] <= d < bars[num, 2]) ? const1 * (d / bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]]) :
        #         (d >= bars[num, 2]) ? const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        # else
        #     V = [(bars[num, 1] < d < time[num][1]) ? const1 * (log(log(d / bars[num, 1])) - eul - lambda) :
        #         (time[num][1] <= d < bars[num, 2]) ? const1 * (log(log(d / bars[num, 1])) - eul - lambda) + dot(time[num] / bars[num, 2], pers[class[num]]) :
        #         (d >= bars[num, 2]) ? const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        end
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if p == 1
                if bars[num, 1] < dom[i] < time[num][1]
                    V[i] = const1 * (dom[i] - bars[num, 1])
                elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                    V[i] = const1 * (dom[i] - bars[num, 1]) + dot(time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
                elseif last(time[num]) <= dom[i] < bars[num, 2]
                    V[i] = const1 * (dom[i] - bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]])
                elseif dom[i] >= bars[num, 2]
                    V[i:end] .= const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]])
                    break
                end
            # elseif p == 2
            #     if bars[num, 1] < dom[i] < time[num][1]
            #         V[i] = const1 * (dom[i] / bars[num, 1])
            #     elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            #         V[i] = const1 * (dom[i] / bars[num, 1]) + dot(time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
            #     elseif last(time[num]) <= dom[i] < bars[num, 2]
            #         V[i] = const1 * (dom[i] / bars[num, 1]) + dot(time[num] / bars[num, 2], pers[class[num]])
            #     elseif dom[i] >= bars[num, 2]
            #         V[i:end] .= const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]])
            #         break
            #     end
            # else
            #     if bars[num, 1] < dom[i] < time[num][1]
            #         V[i] = const1 * (dom[i] - bars[num, 1])
            #     elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
            #         V[i] = const1 * (log(log(dom[i] / bars[num, 1])) - eul - lambda) + dot(time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
            #     elseif last(time[num]) <= dom[i] < bars[num, 2]
            #         V[i] = const1 * (log(log(dom[i] / bars[num, 1])) - eul - lambda) + dot(time[num] / bars[num, 2], pers[class[num]])
            #     elseif dom[i] >= bars[num, 2]
            #         V[i:end] .= const1 * pers[num] + dot(time[num] / bars[num, 2], pers[class[num]])
            #         break
            #     end
            end
        end
    end
    return V
end

function plot_centrality3(num, dom, time, bars, pers, class, const1, const2, p, lambda)
    if isempty(class[num])
        if p == 1
            V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d - bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        elseif p == 2
            V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (d / bars[num, 1]) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        else
            V = [(bars[num, 1] < d <= bars[num, 2]) ? const1 * (log(log(d/ bars[num, 1])) - eul - lambda) : (d > bars[num, 2]) ? const1 * pers[num] : 0 for d in dom]
        end
    elseif length(class[num]) == 1
        if p == 1
            V = [(bars[num, 1] < d < time[num][1]) ? const1 * (d - bars[num, 1]) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (d - bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        elseif p == 2
            V = [(bars[num, 1] < d < time[num][1]) ? const1 * (d / bars[num, 1]) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (d / bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        else
            V = [(bars[num, 1] < d < time[num][1]) ? const1 * (log(log(d / bars[num, 1])) - eul - lambda) :
                (time[num][1] <= d < bars[num, 2]) ? const1 * (log(log(d / bars[num, 1])) - eul - lambda) + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) :
                (d >= bars[num, 2]) ? const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]]) : 0 for d in dom]
        end
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if p == 1
                if bars[num, 1] < dom[i] < time[num][1]
                    V[i] = const1 * (dom[i] - bars[num, 1])
                elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                    V[i] = const1 * (dom[i] - bars[num, 1]) + dot(1 .- time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
                elseif last(time[num]) <= dom[i] < bars[num, 2]
                    V[i] = const1 * (dom[i] - bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                elseif dom[i] >= bars[num, 2]
                    V[i:end] .= const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                    break
                end
            elseif p == 2
                if bars[num, 1] < dom[i] < time[num][1]
                    V[i] = const1 * (dom[i] / bars[num, 1])
                elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                    V[i] = const1 * (dom[i] / bars[num, 1]) + dot(1 .- time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
                elseif last(time[num]) <= dom[i] < bars[num, 2]
                    V[i] = const1 * (dom[i] / bars[num, 1]) + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                elseif dom[i] >= bars[num, 2]
                    V[i:end] .= const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                    break
                end
            else
                if bars[num, 1] < dom[i] < time[num][1]
                    V[i] = const1 * (dom[i] - bars[num, 1])
                elseif (Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                    V[i] = const1 * (log(log(dom[i] / bars[num, 1])) - eul - lambda) + dot(1 .- time[num][1:Y] / bars[num, 2], pers[class[num][1:Y]])
                elseif last(time[num]) <= dom[i] < bars[num, 2]
                    V[i] = const1 * (log(log(dom[i] / bars[num, 1])) - eul - lambda) + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                elseif dom[i] >= bars[num, 2]
                    V[i:end] .= const1 * pers[num] + dot(1 .- time[num] / bars[num, 2], pers[class[num]])
                    break
                end
            end
        end
    end
    return V
end

function transfer1(MC, Bar, Pers, c)
    CD = sort(1:length(MC), by = x -> Bar[x, 2])
    I = [isempty(MC[i]) ? sum(Pers[MC[i]]) : 0 for i in CD]

    for i in 1:length(CD)
        I[i] += sum(I[MC[CD[i]]])
    end

    C = I[sort(1:length(CD), by = x -> CD[x])]
    return C
end

function plot_centrality4(num, dom, time, bar, pers, mc, const1, const2, p, lambda)
    stack = transfer1(mc, bar, pers, time)
    const1pers = const1 * pers[num]
    
    if isempty(mc[num])
        if p == 1
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        elseif p == 2
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d / bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        else
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (log(log(d / bar[num,1])) - eul - lambda) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        end
    end
    
    if length(mc[num]) == 1
        if p == 1
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d - bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d - bar[num, 1]) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif d >= bar[num, 2]
                const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            else
                0
            end
            for d in dom
        ]
        elseif p == 2
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d / bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d / bar[num, 1]) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif d >= bar[num, 2]
                const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            else
                0
            end
            for d in dom
        ]
        else
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (log(log(d / bar[num,1])) - eul - lambda)
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (log(log(d / bar[num,1])) - eul - lambda) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif d >= bar[num, 2]
                const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            else
                0
            end
            for d in dom
        ]
        end
    end

    V = zeros(length(dom))
    for i in 1:length(dom)
        if p == 1
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * (sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] - bar[num, 1]) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
                break
            end
        elseif p == 2
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] / bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] / bar[num, 1]) + const2 * (sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] / bar[num, 1]) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
                break
            end
        else
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda)
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + const2 * (sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + const2 * (sum(pers[mc[num]]) + sum(stack[mc[num]]))
                break
            end
        end        
    end

    return V
end

function transfer2(MC, Bar, Pers, Time)
    CD = sort(1:length(MC), by = x -> Bar[x, 2])
    I = [isempty(MC[i]) ? dot(Time[i]/Bar[i, 2], Pers[MC[i]]) : 0 for i in CD]

    for i in 1:length(CD)
        I[i] += sum(I[MC[CD[i]]])
    end

    C = I[sort(1:length(CD), by = x -> CD[x])]
    return C
end

function plot_centrality5(num, dom, time, bar, pers, mc, const1, const2, p, lambda)
    stack = transfer2(mc, bar, pers, time)
    const1pers = const1 * pers[num]
    
    if isempty(mc[num])
        if p == 1
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        elseif p == 2
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d / bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        else
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (log(log(d / bar[num,1])) - eul - lambda) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        end
    end
    
    if length(mc[num]) == 1
        if p == 1
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d - bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d - bar[num, 1]) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif d >= bar[num, 2]
                const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            else
                0
            end
            for d in dom
        ]
        elseif p == 2
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (d / bar[num, 1])
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (d / bar[num, 1]) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif d >= bar[num, 2]
                const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            else
                0
            end
            for d in dom
        ]
        else
            return [
            if (bar[num, 1] < d < time[num][1])
                const1 * (log(log(d / bar[num,1])) - eul - lambda)
            elseif (time[num][1] <= d < bar[num, 2])
                const1 * (log(log(d / bar[num,1])) - eul - lambda) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif d >= bar[num, 2]
                const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            else
                0
            end
            for d in dom
        ]
        end
    end

    V = zeros(length(dom))
    for i in 1:length(dom)
        if p == 1
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bar[num, 1]) + dot(time[num][1:Y] / bar[num, 2], pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] - bar[num, 1]) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        elseif p == 2
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] / bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] / bar[num, 1]) + dot(time[num][1:Y] / bar[num, 2], pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] / bar[num, 1]) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        else
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda)
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + dot(time[num][1:Y] / bar[num, 2], pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        end        
    end

    return V
end

function transfer3(MC, Bar, Pers, Time)
    CD = sort(1:length(MC), by = x -> Bar[x, 2])
    I = [isempty(MC[i]) ? 0 : dot(1 .- (Time[i] / Bar[i, 2]), Pers[MC[i]]) for i in CD]
    
    for i in 1:length(CD)
        I[i] += sum(I[MC[CD[i]]])
    end
    
    return I[sort(1:length(CD), by = x -> CD[x])]
end

function plot_centrality6(num, dom, time, bar, pers, mc, const1, const2, p, lambda)
    stack = transfer3(mc, bar, pers, time)
    const1pers = const1 * pers[num]
    
    if isempty(mc[num])
        if p == 1
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d - bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        elseif p == 2
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (d / bar[num, 1]) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        else
            return [if (bar[num, 1] < d <= bar[num, 2]) const1 * (log(log(d / bar[num,1])) - eul - lambda) elseif d > bar[num, 2] const1pers else 0 end for d in dom]
        end    
    end
    
    if length(mc[num]) == 1
        if p == 1
            return [
                if (bar[num, 1] < d < time[num][1])
                    const1 * (d - bar[num, 1])
                elseif (time[num][1] <= d < bar[num, 2])
                    const1 * (d - bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                elseif d >= bar[num, 2]
                    const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                else
                    0
                end
                for d in dom
            ]
        elseif p == 2
            return [
                if (bar[num, 1] < d < time[num][1])
                    const1 * (d / bar[num, 1])
                elseif (time[num][1] <= d < bar[num, 2])
                    const1 * (d / bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                elseif d >= bar[num, 2]
                    const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                else
                    0
                end
                for d in dom
            ]
        else
            return [
                if (bar[num, 1] < d < time[num][1])
                    const1 * (log(log(d / bar[num,1])) - eul - lambda)
                elseif (time[num][1] <= d < bar[num, 2])
                    const1 * (log(log(d / bar[num,1])) - eul - lambda) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                elseif d >= bar[num, 2]
                    const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                else
                    0
                end
                for d in dom
            ]
        end        
    end

    V = zeros(length(dom))
    for i in 1:length(dom)
        if p == 1
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] - bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] - bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]][1:Y])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] - bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        elseif p == 2
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (dom[i] / bar[num, 1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (dom[i] / bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]][1:Y])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (dom[i] / bar[num, 1]) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        else
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda)
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x + 1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]][1:Y])
            elseif last(time[num]) <= dom[i] < bar[num, 2]
                V[i] = const1 * (log(log(dom[i] / bar[num,1])) - eul - lambda) + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num, 2]
                V[i:end] .= const1pers + dot(1 .- time[num] / bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
                break
            end
        end
    end

    return V
end

function centrality_plot(bars, pers, time, class, p, lambda, color_vec)
    x = collect(range(0, stop = maximum(bars[:, 2]), length = 2000))
    centrality_plots = []

    for i in 4:6
        centrality_measure = eval(Meta.parse("plot_centrality$i"))
        push!(centrality_plots, plot(x, hcat([centrality_measure(j, x, time, bars, pers, class, 1, 1, p, lambda) for j in 1:length(pers)]...), legend = false, palette = color_vec, ylabel = "J"*subscript(i)*"(σ,ϵ)"))
    end

    return plot(centrality_plots..., layout = 3)
end

function load_data(data_path)
    data = transpose(readdlm(data_path, ',', Float64, '\n', header=true)[1])
    return data[1, :], data[2, :]
end

function save_scatter_plot(x, y, save_path)
    peach_fuzz_color = RGB(1.0, 0.75, 0.60)
    mintgreen_color = RGB(0.7333333333333333, 0.3803921568627451, 0.6784313725490196)
    plot_data = scatter(x, y, legend=false, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, grid=false, background_color=:transparent, foreground_color=:black)
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

function save_diagram_and_measures(bars, pers, time, class, save_diagram_path, save_measures_path, p, lambda, diag, color_vec)
#    savefig(plot(diag), save_diagram_path)
    cent_plots = centrality_plot(bars, pers, time, class, p, lambda, color_vec)
    savefig(cent_plots, save_measures_path)
end

function bootstrap_data(x, y, num_points)
    indices = rand(1:length(x), num_points)
    return x[indices], y[indices]
end

data_path = "./Datasets/sierpinski_points.csv"
save_data_path = "./Datasets/SierpinskiTriangle.png"

data_x, data_y = load_data(data_path)
save_scatter_plot(data_x, data_y, save_data_path)

num_bootstraps = 1000
bootstrap_size = 800
save_bootstrap_dir = "./Sierpinski/BootstrapSamples/"
save_diagram_dir = "./Sierpinski/Diagrams/"
save_measure_dir = "./Sierpinski/Measures/"
save_mapping_dir1 = "./Sierpinski/Map0/"
save_mapping_dir2 = "./Sierpinski/Map1/"
save_mapping_dir3 = "./Sierpinski/Map2/"
peach_fuzz_color = RGB(1.0, 0.745, 0.596)
mintgreen_color = RGB(0.7333333333333333, 0.3803921568627451, 0.6784313725490196)

sgnl1_len = []
sgnl_len = []

for indx in 1:1
    println(indx)
    x, y = bootstrap_data(data_x, data_y, bootstrap_size)
    
    save_bootstrap_path = joinpath(save_bootstrap_dir, "BootstrapSample_$indx.csv")
    writedlm(save_bootstrap_path, hcat(x, y), ',')

    max_dim = 1
    data = transpose(hcat([x, y][1], [x, y][2]))
    bars, reps, pers, pers_mult, pers_trans, class, time, diag = process_data(data, max_dim)
    equiv = setdiff(1:length(pers), union(class...))
    lambda = sum(log.(log.(pers_mult)))/length(pers_mult)
    color_vec = distinguishable_colors(length(pers))
    
    save_diagram_path = joinpath(save_diagram_dir, "Diagram_$indx.png")
    save_measures_path = joinpath(save_measure_dir, "Measure_$indx.png")
    save_diagram_and_measures(bars, pers, time, class, save_diagram_path, save_measures_path, 1, lambda, diag, color_vec)
    savefig(plot(diag, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, grid=false, background_color=:transparent, foreground_color=:black), save_diagram_path)

    # Cycle Plotting Based on Barcode
    sgf = 0.05
    sgnl = find_all(x -> exp(-exp(x)) < (sgf/length(pers)), pers_trans)
    
    if sgnl != nothing
        push!(sgnl_len, length(sgnl))
        Plots.plot(data[1,:], data[2,:], seriestype = :scatter, legend = false, grid = false, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, background_color=:transparent, foreground_color=:black)
        for k in sgnl
            cycl = data[:, unique(reps[k])]
            Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
            for i in 1:size(reps[k],2)
                Plots.plot!(data[1, [reps[k][1, i], reps[k][2, i]]], data[2, [reps[k][1, i], reps[k][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
            end
        end
    else
        push!(sgnl_len, 0)
    end

    save_mapping_path = joinpath(save_mapping_dir1, "Mapping_$indx.png")
    savefig(save_mapping_path)
    
    dom = collect(range(0, stop = maximum(bars[:,2]), length = 2000))
    cent5 = hcat([plot_centrality5(j, dom, time, bars, pers, class, 1, 1, 1, lambda) for j in 1:length(pers)]...)
    max_cent = [maximum(cent5[:, i]) for i in 1:size(cent5, 2)]

    println(sgnl_len[indx])

    if sgnl_len[indx] > 0
        sgnl1 = setdiff(find_all(x -> x >= minimum(max_cent[sgnl]), max_cent), equiv)
    else
        sgnl1 = setdiff(find_all(x -> x >= 0.8*maximum(max_cent), max_cent), equiv)
    end 

    # if sgnl1 != nothing
    #     push!(sgnl1_len, length(sgnl1))
    #     for k in sgnl1
    #         j = k #cycl = data[:, unique(reps[k])]
    #         while findfirst(x -> x in class[j], 1:j-1) != nothing
    #             j = findfirst(x -> x in class[j], 1:j-1)
    #         end
    #         cycl = data[:, unique(reps[j])]
    #         Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
    #         for l in 1:size(reps[j],2)
    #             Plots.plot!(data[1, [reps[j][1, l], reps[j][2, l]]], data[2, [reps[j][1, l], reps[j][2, l]]], linecolor = color_vec[k], linewidth = 4.0)
    #         end
    #     end
    # else
    #     push!(sgnl1_len, 0)
    # end
    
    # save_mapping_path1 = joinpath(save_mapping_dir2, "Mapping_$indx.png")
    # savefig(save_mapping_path1)

    # if sgnl1 != nothing
    #     for k in sgnl1
    #         cycl = data[:, unique(reps[k])]
    #         Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
    #         for i in 1:size(reps[k],2)
    #             Plots.plot!(data[1, [reps[k][1, i], reps[k][2, i]]], data[2, [reps[k][1, i], reps[k][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
    #         end
    #     end
    # end
        
    # save_mapping_path2 = joinpath(save_mapping_dir3, "Mapping_$indx.png")
    # savefig(save_mapping_path2)

end

# writedlm("Signals.csv",  sgnl_len, ',')
# writedlm("Signals1.csv",  sgnl1_len, ',')

# max_dim = 1
# data = transpose(hcat([x, y][1], [x, y][2]))
# bars, reps, pers, pers_mult, pers_trans, class, time, diag = process_data(data, max_dim)
# print(class)
# peach_fuzz_color = RGB(1.0, 0.745, 0.596)
# mintgreen_color = RGB(0.7333333333333333, 0.3803921568627451, 0.6784313725490196)
# savefig(plot(diag, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, grid=false, background_color=:transparent, foreground_color=:black), save_diagram_path)
# lambda = sum(log.(log.(pers_mult)))/length(pers_mult)
# color_vec = distinguishable_colors(length(pers))
# save_diagram_and_measures(bars, pers, time, class, save_diagram_path, save_measures_path, 1, lambda)

# # Cycle Plotting Based on Barcode
# sgf = 0.05
# sgnl = find_all(x -> exp(-exp(x)) < (sgf/length(pers)), pers_trans)
# #println(sgnl)
# Plots.plot(data[1,:], data[2,:], seriestype = :scatter, legend = false, grid = false, markersize=4, markercolor=peach_fuzz_color, markerstrokecolor=mintgreen_color, background_color=:transparent, foreground_color=:black)

# for k in sgnl
#     cycl = data[:, unique(reps[k])]
#     Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
#     for i in 1:size(reps[k],2)
#         plt = Plots.plot!(data[1, [reps[k][1, i], reps[k][2, i]]], data[2, [reps[k][1, i], reps[k][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
#     end
# end

# savefig(save_mapping_path)

# x = collect(range(0, stop = maximum(bars[:,2]), length = 2000))
# cent5 = hcat([plot_centrality5(j, x, time, bars, pers, class, 1, 1, 1, lambda) for j in 1:length(pers)]...);
# y = [maximum(cent5[:, i]) for i in 1:size(cent5, 2)]
# sgnl1 = findall(x -> x >= minimum(y[sgnl]), y)

# for k in sgnl1
#     j = k #    cycl = data[:, unique(reps[k])]
#     #println(class[j])
#     while findfirst(x -> j in class[x], 1:j-1) != nothing
#         #println(j)
#         j = findfirst(x -> j in class[x], 1:j-1)
#     end
#     cycl = data[:, unique(reps[j])]
#     Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
#     for i in 1:size(reps[j],2)
#         plt = Plots.plot!(data[1, [reps[j][1, i], reps[j][2, i]]], data[2, [reps[j][1, i], reps[j][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
#     end
# end

# savefig(save_mapping_path1)

# for k in sgnl1
#     cycl = data[:, unique(reps[k])]
#     Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
#     for i in 1:size(reps[k],2)
#         plt = Plots.plot!(data[1, [reps[k][1, i], reps[k][2, i]]], data[2, [reps[k][1, i], reps[k][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
#     end
# end

# savefig(save_mapping_path2)