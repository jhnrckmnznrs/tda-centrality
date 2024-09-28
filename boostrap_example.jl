include("main.jl")

data_path = "./Datasets/sierpinski_points.csv"
save_data_path = "./Datasets/SierpinskiTriangle.png"

data_x, data_y = load_data(data_path)
save_scatter_plot(data_x, data_y, save_data_path)

num_bootstraps = 1000
bootstrap_size = 800
save_bootstrap_dir = "./Sierpinski/BootstrapSamples/"
save_bootstrap_plots_dir = "./Sierpinski/BootstrapPlots/"
save_diagram_dir = "./Sierpinski/Diagrams/"
save_measure_dir = "./Sierpinski/Measures/"
save_mapping_dir1 = "./Sierpinski/Map0/"
save_mapping_dir2 = "./Sierpinski/Map1/"
save_mapping_dir3 = "./Sierpinski/Map2/"
peach_fuzz_color = RGB(1.0, 0.745, 0.596)
mintgreen_color = RGB(0.7333333333333333, 0.3803921568627451, 0.6784313725490196)

sgnl_len = []

for indx in 1:1000
    println(indx)
    x, y = bootstrap_data(data_x, data_y, bootstrap_size)
    
    save_bootstrap_path = joinpath(save_bootstrap_dir, "BootstrapSample_$indx.csv")
    save_bootstrap_data_path = joinpath(save_bootstrap_plots_dir, "BootstrapSample_$indx.png")
    save_scatter_plot(x, y, save_bootstrap_data_path)
    writedlm(save_bootstrap_path, hcat(x, y), ',')

    max_dim = 1
    data = transpose(hcat([x, y][1], [x, y][2]))
    bars, reps, pers, pers_mult, pers_trans, class, time, diag = process_data(data, max_dim)
    lambda = sum(log.(log.(pers_mult)))/length(pers_mult)
    color_vec = distinguishable_colors(length(pers))
    
    save_diagram_path = joinpath(save_diagram_dir, "Diagram_$indx.png")
    save_measures_path = joinpath(save_measure_dir, "Measure_$indx.png")
    save_diagram_and_measures(bars, pers, time, class, save_diagram_path, save_measures_path, diag, color_vec)
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
    cent5 = hcat([plot_centrality5(j, dom, time, bars, pers, class, 1, 1) for j in 1:length(pers)]...)
    max_cent = [maximum(cent5[:, i]) for i in 1:size(cent5, 2)]

    sgnl1 = find_all(x -> x >= 0.8*maximum(max_cent), max_cent)
    
    if sgnl1 != nothing
        for k in sgnl1
            j = k
            while findfirst(x -> x in class[j], 1:j-1) != nothing
                j = findfirst(x -> x in class[j], 1:j-1)
            end
            cycl = data[:, unique(reps[j])]
            Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
            for l in 1:size(reps[j],2)
                Plots.plot!(data[1, [reps[j][1, l], reps[j][2, l]]], data[2, [reps[j][1, l], reps[j][2, l]]], linecolor = color_vec[k], linewidth = 4.0)
            end
        end
    end
    
    save_mapping_path1 = joinpath(save_mapping_dir2, "Mapping_$indx.png")
    savefig(save_mapping_path1)

    if sgnl1 != nothing
        for k in sgnl1
            cycl = data[:, unique(reps[k])]
            Plots.plot!(cycl[1,:], cycl[2,:], seriestype = :scatter, legend = false, markersize = 5, markercolor = color_vec[k])
            for i in 1:size(reps[k],2)
                Plots.plot!(data[1, [reps[k][1, i], reps[k][2, i]]], data[2, [reps[k][1, i], reps[k][2, i]]], linecolor = color_vec[k], linewidth = 4.0)
            end
        end
    end
        
    save_mapping_path2 = joinpath(save_mapping_dir3, "Mapping_$indx.png")
    savefig(save_mapping_path2)

end

writedlm("Signals.csv",  sgnl_len, ',')