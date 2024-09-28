# Initialize Barcode (with Persistence) and Cycle Representatives (Specific Dimension)
function barcode_dim(p::Dict{String,Any}, d::Int64)
    B = Eirene.barcode(p, dim = d) # Barcode from PersHom
    OB = sort(1:size(B,1), by = x -> B[x,1]) # Sort by birth
    B = B[OB, 1:2] # Sorted Barcode by Birth
    return B, [classrep(p, dim = d, class = i) for i in OB], -(B[:,2], B[:,1])
end

# Merging Cycles
function find_all(f, a::Array{T, N}) where {T, N} # Faster and Efficient Altenative to 'findall' (c) Julia Language Forum
    j = 1
    b = Vector{Int}(undef, length(a))
     @inbounds for i in eachindex(a)
        @inbounds if f(a[i])
            b[j] = i
            j += 1
        end
    end
    resize!(b, j-1)
    sizehint!(b, length(b))
    return b
end

function any_equal_col(a::Array{Int64,2}, b::Array{Int64,2}) # Returns a Boolean 'true' if Some Rows of Two Arrays are Equal otherwise 'False'
    k = 0
    for i in 1:size(a, 2)
        for j in 1:size(b,2)
            if sort(a[:,i]) == sort(b[:,j])
                k = 1
                break
            end
        end
    end
    return k
end

function adjacency(a, b::Array{Int64,2}) # Returns a Boolean 'true' when Some Rows of Two Arrays are Equal
    k = 0
    for i in 1:length(a)
        if any_equal_col(a[i], b) == 1
            k = 1
            break
        end
    end
    return k
end

function merging_process(b::Matrix{Float64}, c::Vector{Matrix{Int64}}) # Returns Array where Coordinate i Corresponds to Cycles that Merge with cycle i
    s = size(b,1)
    a = [find_all(x -> any_equal_col(c[i], c[x]) == 1, collect(1:s)) for i in 1:s]
    m = [find_all(x -> any_equal_col(c[i], c[x]) == 1 && b[:,2][i] > b[:,2][x] > b[:,1][i], collect(1:s)) for i in 1:s]
    index = find_all(x -> x != [], m)
    for i in index
        possible = filter(x -> x ∉ m[i], collect(1:s))
        while (local n = possible[find_all(x -> adjacency(c[m[i]], c[x]) == 1 && b[:,2][i] > b[:,2][x] > b[:,1][i], possible)]) != []
            union!(m[i], n)
            possible = filter(x -> x ∉ m[i], collect(1:s))
        end
    end
    for i in 1:(length(m) - 1)
        [filter!((x) -> x ∉ intersect(m[i],m[j]), m[j]) for j in i+1:length(m)]
    end
    [sort!(i, by = x -> b[x,2]) for i in m]
    t = [b[i,2] for i in m]
    return m, t
end

#Centrality Measures
function transfer(MC, Bar, Time, Pers, num)
    CD = sort(1:length(MC), by = x -> Bar[x,2])
    if num == 4
        I = [if MC[i] != [] sum(Pers[MC[i]]) else 0 end for i in CD]
        for i in 1:length(CD)
            y = findall(x -> x ∈ MC[CD[i]], CD)
            # I[i] = I[i] + sum(I[MC[CD[i]]])
            I[i] = I[i] + sum(I[y])
        end
        C = I[sort(1:length(CD), by = x -> CD[x])]
    elseif num == 5
        I = [if MC[i] !=[] dot(Time[i]/Bar[i,2], Pers[MC[i]]) else 0 end for i in CD]
        for i in 1:length(CD)
            y = findall(x -> x ∈ MC[CD[i]], CD)
            # I[i] = I[i] + sum(I[MC[CD[i]]])
            I[i] = I[i] + sum(I[y])
            # I[i] = I[i] + sum(I[MC[CD[i]]])
        end
        C = I[sort(1:length(CD), by = x -> CD[x])]
    else
        I = [if MC[i] !=[] dot(1 .- (Time[i]/Bar[i,2]), Pers[MC[i]]) else 0 end for i in CD]
        for i in 1:length(CD)
            y = findall(x -> x ∈ MC[CD[i]], CD)
            # I[i] = I[i] + sum(I[MC[CD[i]]])
            I[i] = I[i] + sum(I[y])
            # I[i] = I[i] + sum(I[MC[CD[i]]])
        end
        C = I[sort(1:length(CD), by = x -> CD[x])]
    end
    return C
end

function plot_centrality(mc, time, bar, pers, cnum)
    if cnum == 1
        const1 = 1
        V = [if (mc[i] == Any[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + sum(pers[mc[i]]) else z - bar[i, 1] + const1*(sum(pers[mc[i][1:j - 1]])) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    elseif cnum == 2
        const1 = [if (mc[i] != Int64[]) time[i]/bar[i,2] end for i in 1:length(mc)]
        V = [if (mc[i] == Int64[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + dot(const1[i], pers[mc[i]])  else z - bar[i, 1] + dot(const1[i][1:j-1], pers[mc[i][1:j-1]]) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    elseif cnum == 3
        const1 = [if (mc[i] != Int64[]) 1 .- time[i]/bar[i,2] end for i in 1:length(mc)]
        V = [if (mc[i] == Int64[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + dot(const1[i], pers[mc[i]])  else z - bar[i, 1] + dot(const1[i][1:j-1], pers[mc[i][1:j-1]]) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    elseif cnum == 4
        stack = transfer(mc, bar, time, pers, cnum)
        const1 = 1
        V = [if (mc[i] == Int64[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + sum(pers[mc[i]]) + sum(stack[mc[i]])  else z - bar[i, 1] + sum(pers[mc[i][1:j-1]]) + sum(stack[mc[i][1:j-1]]) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    elseif cnum == 5
        stack = transfer(mc, bar, time, pers, cnum)
        const1 = [if (mc[i] != Int64[]) time[i]/bar[i,2] end for i in 1:length(mc)]
        V = [if (mc[i] == Int64[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + dot(const1[i], pers[mc[i]]) + sum(stack[mc[i]])  else z - bar[i, 1] + dot(const1[i][1:j-1], pers[mc[i][1:j-1]]) + sum(stack[mc[i][1:j-1]]) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    elseif cnum == 6
        stack = transfer(mc, bar, time, pers, cnum)
        const1 = [if (mc[i] != Int64[]) 1 .- time[i]/bar[i,2] end for i in 1:length(mc)]
        V = [if (mc[i] == Int64[])
        [if (k == 1) z - bar[i, 1] else pers[i] end for k in 1:2] else
        [if (j == 1) z - bar[i, 1] elseif (j == length(time[i])+2) pers[i] + dot(const1[i], pers[mc[i]]) + sum(stack[mc[i]])  else z - bar[i, 1] + dot(const1[i][1:j-1], pers[mc[i][1:j-1]]) + sum(stack[mc[i][1:j-1]]) end for j in 1:length(time[i]) + 2]
        end for i in 1:length(mc)]
    end
    return V
end

function bottleneck0(X, Y)
    if length(Y) < length(X)
        X0 = copy(X)
        X = Y
        Y = X0
    end
    X = -(sort(-X))
    Y = -(sort(-Y))
    d = 0
    N = length(X)
    Z = abs.(X[1:N]-Y[1:N])
    l = findfirst(x -> x == maximum(Z), Z)
    dtemp = Z[l]
    if N != length(Y) && dtemp < 0.5*Y[N]
        d = 0.5*Y[N]
    elseif l >= 1 && length(Z) > 1
        while length(Z) > 1
            k = 0.5*max(X[l], Y[l])
            if maximum(Z[setdiff(1:length(Z), l)]) < k && k < dtemp
                d = k
                break
            elseif maximum(Z[setdiff(1:length(Z), l)]) >= k
                if length(Z[findall(x -> x >= k, Z)]) == length(findall(y -> y >= l, findall(x -> x >= k, Z)))
                    d = k
                    break
                else
                    Z = Z[1:(l-1)]
                    X = X[1:(l-1)]
                    Y = Y[1:(l-1)]
                    l = findfirst(x -> x == maximum(Z), Z)
                    dtemp = Z[l]
                    if length(Z) == 1
                        d = min(dtemp, 0.5*max(X[l], Y[l]))
                        break
                    end
                end
            else
                d = dtemp
                break
            end
        end
    else
        d = min(dtemp, 0.5*max(X[1], Y[1]))
    end
    return d
end

function metric(cp2, mt2, bar2, area1)
    # int2 = [sort(union(mt2[i], bar2[i,1:2])) for i in 1:length(mt2)]
    # area2 = [sum([0.5*(int2[i][j] - int2[i][j - 1])^2 + (subs(cp2[i][j-1], z, int2[i][j-1])*(int2[i][j] - int2[i][j - 1]))  for j in 2:length(int2[i])]) for i in 1:length(int2)]
    area2 = og_area(mt2, bar2, cp2)
    return bottleneck0(area1, area2)
end

function og_area(mt1, bar1, cp1)
    int1 = [sort(union(mt1[i], bar1[i,1:2])) for i in 1:length(mt1)]
    area1 = [sum([0.5*(int1[i][j] - int1[i][j - 1])^2 + (subs(cp1[i][j-1], z, int1[i][j-1])*(int1[i][j] - int1[i][j - 1]))  for j in 2:length(int1[i])]) for i in 1:length(int1)]
    return convert(Vector{Float64}, area1)
end

function stability(data, dim, area, cnum)
    phom = eirene(data, model = "pc", maxdim = dim)
    bar, crep, pers = barcode_dim(phom, dim)
    mc, mt = merging_process(bar, crep)
    cp = plot_centrality(mc, mt, bar, pers, cnum)
    return metric(cp, mt, bar, area)
end

function centrality_map(mt, cp1, bars, pers)
    x = sort(union(mt...))
    if length(x) > 1
        x0 = minimum([x[i] - x[i - 1] for i in 2:length(x)])
        xs = minimum(bars[:,1]):x0:maximum(bars[:,2])
        if ((maximum(bars[:,2]) - minimum(bars[:,1])) % x0) != 0
            xs = union(xs, last(xs) + x0)
        end
    elseif length(x) == 1
        xs =  [minimum(bars[:,1]), x[1], maximum(bars[:,2])]
    else
        xs = [minimum(bars[:,1]), maximum(bars[:,2])]
    end
    ys = ["σ"*subscript(i) for i = 1:length(mt)]
    CM = [if (mt[i] == [])
            [if (j < bars[i, 1])
                0
            elseif (j >= bars[i, 2])
                last(cp1[i])
            else
                subs(cp1[i][1], z, j)
            end for j in xs]
        else
            [if (j < bars[i, 1])
                0
            elseif (j >= bars[i, 2])
                last(cp1[i])
            elseif (last(mt[i]) <= j < bars[i, 2])
                subs(cp1[i][length(cp1[i]) - 1], z, j)
            else
                subs(cp1[i][findfirst(x -> j < mt[i][x], 1:length(mt[i]))], z, j)
            end for j in xs]
        end for i in 1:length(mt)]
    # return CM
    return Plots.heatmap(ys, xs, transpose(vcat([transpose(i) for i in CM]...)), background_color = :transparent, foreground_color = :black, c=[:black, :white, :aliceblue, :lightskyblue, :yellow1, :darkred])
    # return heatmap(ys, xs, CM, c=[:black, :white, :aliceblue, :lightskyblue, :yellow1, :darkred])
end

function plot_centrality_map(mt, cp, bar, pers)
    cm1 = centrality_map(mt, cp[1], bar, pers);
    cm2 = centrality_map(mt, cp[2], bar, pers);
    cm3 = centrality_map(mt, cp[3], bar, pers);
    cm4 = centrality_map(mt, cp[4], bar, pers);
    cm5 = centrality_map(mt, cp[5], bar, pers);
    cm6 = centrality_map(mt, cp[6], bar, pers);
    return HM = Plots.plot(cm1,cm2,cm3,cm4,cm5,cm6, layout = 6, title = ["J"*subscript(1)*"(σ,ϵ)" "J"*subscript(2)*"(σ,ϵ)" "J"*subscript(3)*"(σ,ϵ)" "J"*subscript(4)*"(σ,ϵ)" "J"*subscript(5)*"(σ,ϵ)" "J"*subscript(6)*"(σ,ϵ)"]);
end

function plot_box(cpd)
    b1 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[1]...), leg = false);
    b2 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[2]...), leg = false);
    b3 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[3]...), leg = false);
    b4 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[4]...), leg = false);
    b5 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[5]...), leg = false);
    b6 = StatsPlots.boxplot(reshape([if i == 1 string(0.1) else string(0.25*(i-1)) end for i in 1:41], (1,41)), hcat(cpd[6]...), leg = false);
    return Plots.plot(b1,b2,b3,b4,b5,b6, layout = 6, title=["J₁" "J₂ (late)" "J₂ (early)" "J₃ (f = 1)" "J₃ (late)" "J₃ (early)"], size = (1200, 800))
end

function bottleneck_experiment(data, pd, dim)
    phom = eirene(data, model = "pc", maxdim = dim)
    bar, crep, pers = barcode_dim(phom, dim)
    pd1 = PersistenceDiagram([(bar[i,1], bar[i,2]) for i in 1:size(bar,1)])
    return Bottleneck()(pd, pd1)
end

# Centrality Plots (Original Representation)
function plot_centrality1(num, dom, time, bar, pers, mc, const1, const2)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + const2*sum(pers[mc[num]])
            elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + const2*sum(pers[mc[num]]) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + const2*sum(pers[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + const2*sum(pers[mc[num]])
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + const2*sum(pers[mc[num]])], length(dom) - i + 1)
                break
            end
        end
    end
    return V
end

## Centrality 2
function plot_centrality2(num, dom, time, bar, pers, mc, const1)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + ((time[num]/bar[num,2]) ⋅ pers[mc[num]])
            elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + ((time[num]/bar[num,2]) ⋅ pers[mc[num]]) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + ((time[num][1:Y]/bar[num,2]) ⋅ pers[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + ((time[num]/bar[num,2]) ⋅ pers[mc[num]])
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + ((time[num]/bar[num,2]) ⋅ pers[mc[num]])], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

## Centrality 3
function plot_centrality3(num, dom, time, bar, pers, mc, const1)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + ((1 .-(time[num]/bar[num,2])) ⋅ pers[mc[num]])
            elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + ((1 .-(time[num]/bar[num,2])) ⋅ pers[mc[num]])  else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + ((1 .- (time[num][1:Y]/bar[num,2])) ⋅ pers[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + ((1 .-(time[num]/bar[num,2])) ⋅ pers[mc[num]])
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + ((1 .-(time[num]/bar[num,2])) ⋅ pers[mc[num]])], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

## Centrality 4
function plot_centrality4(num, dom, time, bar, pers, mc, const1, const2)
    stack = transfer1(mc, bar, pers, const2)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
    elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]])) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

function plot_centrality41(num, dom, time, bar, pers, mc, const1, const2)
    stack = transfer1(mc, bar, pers, const2)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i]/bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i]/bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i]/bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
    elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]])) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i]/bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i]/bar[num,1]) + const2*(sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i]/bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

function plot_centrality42(num, dom, time, bar, pers, mc, const1, const2)
    stack = transfer1(mc, bar, pers, const2)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
    elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]])) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]]))
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + const2*(sum(pers[mc[num]]) + sum(stack[mc[num]]))], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

function transfer1(MC, Bar, Pers, c)
    CD = sort(1:length(MC), by = x -> Bar[x,2])
    I = [if MC[i] != [] sum(Pers[MC[i]]) else 0 end for i in CD]
    for i in 1:length(CD)
        I[i] = I[i] + sum(I[MC[CD[i]]])
    end
    C = I[sort(1:length(CD), by = x -> CD[x])]
    return C
end

## Centrality 5
function plot_centrality5(num, dom, time, bar, pers, mc, const1)
    stack = transfer2(mc, bar, pers, time)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + dot(time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
    elseif (dom[i] >= bar[num,2]) const1*(pers[num]) + dot(time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])  else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + dot(time[num][1:Y]/bar[num, 2], pers[mc[num][1:Y]]) + sum(stack[mc[num][1:Y]])
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + dot(time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + dot(time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

function transfer2(MC, Bar, Pers, Time)
    CD = sort(1:length(MC), by = x -> Bar[x,2])
    I = [if MC[i] !=[] dot(Time[i]/Bar[i,2], Pers[MC[i]]) else 0 end for i in CD]
    for i in 1:length(CD)
        I[i] = I[i] + sum(I[MC[CD[i]]])
    end
    C = I[sort(1:length(CD), by = x -> CD[x])]
    return C
end

## Centrality 6
function plot_centrality6(num, dom, time, bar, pers, mc, const1)
    stack = transfer3(mc, bar, pers, time)
    if mc[num] == Int64[]
        V = [if (bar[num,1] < dom[i] <= bar[num,2]) const1*(dom[i] - bar[num,1]) elseif (dom[i] > bar[num,2]) const1*(pers[num]) else 0 end for i in 1:length(dom)]
    elseif length(mc[num]) == 1
        V = [if (bar[num,1] < dom[i] < time[num][1]) const1*(dom[i] - bar[num,1]) elseif (time[num][1] <= dom[i] < bar[num,2]) const1*(dom[i] - bar[num,1]) + dot(1 .- time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
    elseif (dom[i] >= bar[num,2]) const1*pers[num] + dot(1 .- time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]]) else 0 end for i in 1:length(dom)]
    else
        V = zeros(length(dom))
        for i in 1:length(dom)
            if bar[num, 1] < dom[i] < time[num][1]
                V[i] = const1*(dom[i] - bar[num,1])
            elseif (local Y = findfirst(x -> time[num][x] <= dom[i] < time[num][x+1], 1:length(time[num]) - 1)) != nothing
                V[i] = const1*(dom[i] - bar[num,1]) + dot(1 .- time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]][1:Y])
            elseif last(time[num]) <= dom[i] < bar[num,2]
                V[i] = const1*(dom[i] - bar[num,1]) + dot(1 .- time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])
            elseif dom[i] >= bar[num,2]
                V[i:length(dom)] = repeat([const1*(pers[num]) + dot(1 .- time[num]/bar[num, 2], pers[mc[num]]) + sum(stack[mc[num]])], length(dom)-i+1)
                break
            end
        end
    end
    return V
end

function transfer3(MC, Bar, Pers, Time)
    CD = sort(1:length(MC), by = x -> Bar[x,2])
    I = [if MC[i] !=[] dot(1 .- (Time[i]/Bar[i,2]), Pers[MC[i]]) else 0 end for i in CD]
    for i in 1:length(CD)
        I[i] = I[i] + sum(I[MC[CD[i]]])
    end
    C = I[sort(1:length(CD), by = x -> CD[x])]
    return C
end

function centrality_plot(bar, pers, mt, mc, p)
    x = collect(range(0, stop = maximum(bar[:,2]), length = Int(ceil(maximum(bar[:,2])/minimum(p)*1))))
    p1 = Plots.plot(x, hcat([plot_centrality1(i, x, mt, bar, pers, mc, 1, 1) for i in 1:length(pers)]...), legend = false,  ylabel = "J"*subscript(1)*"(σ,ϵ)", palette = palette(:tab20c, size(bar, 1)))
    p2 = Plots.plot(x, hcat([plot_centrality2(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...), legend = false, palette = palette(:tab20c, size(bar, 1)), ylabel = "J"*subscript(2)*"(σ,ϵ)")
    p3 = Plots.plot(x, hcat([plot_centrality3(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...), legend = false, palette = palette(:tab20c, size(bar, 1)),  ylabel = "J"*subscript(3)*"(σ,ϵ)")
    p4 = Plots.plot(x, hcat([plot_centrality4(i, x, mt, bar, pers, mc, 1, 1) for i in 1:length(pers)]...), legend = false, palette = palette(:tab20c, size(bar, 1)),  ylabel = "J"*subscript(4)*"(σ,ϵ)")
    p5 = Plots.plot(x, hcat([plot_centrality5(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...), legend = false, palette = palette(:tab20c, size(bar, 1)),  ylabel = "J"*subscript(5)*"(σ,ϵ)")
    p6 = Plots.plot(x, hcat([plot_centrality6(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...), legend = false, palette = palette(:tab20c, size(bar, 1)),  ylabel = "J"*subscript(6)*"(σ,ϵ)")
    return Plots.plot(p1,p2,p3,p4,p5,p6, layout = 6)
end

#RdYlGn_11 previous
#cubehelix current

function maximum_centrality(bar, pers, mt, mc)
    x = collect(range(0, stop = maximum(bar[:,2]), length = 10))
    p1 = hcat([plot_centrality1(i, x, mt, bar, pers, mc, 1, 1) for i in 1:length(pers)]...)
    p2 = hcat([plot_centrality2(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...);
    p3 = hcat([plot_centrality3(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...);
    p4 = hcat([plot_centrality4(i, x, mt, bar, pers, mc, 1, 1) for i in 1:length(pers)]...);
    p5 = hcat([plot_centrality5(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...);
    p6 = hcat([plot_centrality6(i, x, mt, bar, pers, mc, 1) for i in 1:length(pers)]...);
    y1 = [maximum(p1[:, i]) for i in 1:size(p1, 2)]
    y2 = [maximum(p2[:, i]) for i in 1:size(p2, 2)]
    y3 = [maximum(p3[:, i]) for i in 1:size(p3, 2)]
    y4 = [maximum(p4[:, i]) for i in 1:size(p4, 2)]
    y5 = [maximum(p5[:, i]) for i in 1:size(p5, 2)]
    y6 = [maximum(p6[:, i]) for i in 1:size(p6, 2)]
    return [y1 y2 y3 y4 y5 y6]
end

function landscape(bar, dim)
    A = sortslices(sortslices(bar, dims = dim, by = x -> x[2], rev = true), dims = 1, by = x -> x[1])
    k = 1
    L = []
    while size(A, 1) != 0
        b = A[1,1]
        d = A[1,2]
        A = A[1:size(A, 1) .!= 1, :]
        p = 1
        LP = [[-Inf, 0], [b, 0], [(b + d)/2, (d - b)/2]]
        while LP[length(LP)] != [Inf, 0]
            if (size(A, 1)) == 0 || (p >= size(A,1))
                LP = vcat(LP, [[d, 0], [Inf, 0]])
                break
            end
            if d >= maximum(A[p:size(A,1), 2])
                LP = vcat(LP, [[d, 0], [Inf, 0]])
            else
                p = collect(p:size(Barcode, 1))[findfirst(x -> A[x, 2] > d, p:size(A, 1))]
                b1 = A[p,1]
                d1 = A[p,2]
                A = A[1:size(A, 1) .!= p, :]
                if b1 > d
                    LP = vcat(LP, [[d, 0]])
                end
                if b1 >= d
                    LP = vcat(LP, [[b1, 0]])
                else
                    LP = vcat(LP, [[(b1 + d)/2, (d - b1)/2]])
                    A = vcat(A, [b1 d])
                    A = vcat(A[1:p-1, :], sortslices(sortslices(A[p:size(A,1), :], dims = 1, by = x -> x[2], rev = true), dims = 1, by = x -> x[1]))
                end
                LP = vcat(LP, [[(b1 + d1)/2, (d1 - b1)/2]])
                b = b1
                d = d1
            end
        end
        L = vcat(L, [LP])
    end
    return L
end

function landscape_plot(L)
    land = hcat(L[1]...)
    P = plot(land[1, :], land[2, :], legend = false)
    if length(L) > 1
        for i in 1:length(L)
            land = hcat(L[i]...)
            P = plot!(land[1, :], land[2, :], legend = false)
        end
    end
    xlabel!("x");
    ylabel!("λ(x)");
    return P
end

function q_zero(mc, bar, num)
    if num == 1
        k = maximum([length(i) for i in mc])
    else
        k = count_cluster(mc,bar)
    end
    return k
end

function count_cluster(m, b)
    cd = sort(1:length(m), by = x -> b[x,2])
    I = [if m[i] != [] length(m[i]) else 0 end for i in cd]
    for i in 1:length(cd)
        y = findall(x -> x ∈ m[cd[i]], cd)
        I[i] = I[i] + sum(I[y])
    end
    return I[sort(1:length(cd), by = x -> cd[x])]
end

# function metric(cp1, cp2, mt1, mt2, bar1, bar2)
#     ind1 = sort(1:length(cp1), by = x -> last(cp1[x]))
#     ind2 = sort(1:length(cp2), by = x -> last(cp2[x]))
#     interv = [if (i <= min(length(ind1), length(ind2))) sort(union(mt1[ind1[i]],mt2[ind2[i]],bar1[ind1[i],1:2], bar2[ind2[i],1:2])) elseif (i > min(length(ind1), length(ind2)) && length(ind2) > length(ind1)) sort(union(mt2[ind2[i]], bar2[ind2[i],1:2])) elseif (i > min(length(ind1), length(ind2)) && length(ind1) > length(ind2)) sort(union(mt1[ind1[i]], bar1[ind1[i],1:2])) end for i in 1:max(length(ind1), length(ind2))]
#     f1 = [[if (interv[i][j] < bar1[ind[i],1]) z elseif (interv[i][j] > bar1[ind[i],2]) last(cp1[ind[i]]) else cp1[ind1][findfirst(x -> interv[i][j] < mt1[ind[i]][x], 1:length(mt1[ind[i]]))] end for j in 1:length(interv[i])-1] for i in 1:length(ind1)]
#     f2 = [[if (interv[i][j] < bar2[ind[i],1]) z elseif (interv[i][j] > bar2[ind[i],2]) last(cp2[ind[i]]) else cp2[ind1][findfirst(x -> interv[i][j] < mt2[ind[i]][x], 1:length(mt2[ind[i]]))] end for j in 1:length(interv[i])-1] for i in 1:length(ind2)]
# end
