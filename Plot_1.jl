using DelimitedFiles
using Plots
using Interpolations
using StatsBase

# --- Função para transformar histograma em função contínua (interpolada)
function histogram_to_function(data::Vector{Float64}; nbins::Int64 = 50)
    h = fit(Histogram, data; nbins=nbins, closed=:left)

    # Converte os pesos para Float64 antes de normalizar
    weights = Float64.(h.weights)
    bin_width = diff(h.edges[1])[1]
    weights ./= sum(weights) * bin_width

    bin_centers = 0.5 .* (h.edges[1][1:end-1] .+ h.edges[1][2:end])
    itp = LinearInterpolation(bin_centers, weights, extrapolation_bc=Flat())
    return itp, extrema(bin_centers)
end

# --- Carrega os dados dos arquivos txt
energy = Float64.(readdlm("/home/alberto/Downloads/Dados/energy30.txt"))
#energy = Float64.(readdlm("/home/alberto/Downloads/Dados/inner_product.txt"))

# --- Divide em três amostras de 3000 pontos cada
data3 = energy[1:30000]
data2 = energy[30001:60000]
data1 = energy[60001:90000]

# --- Converte histogramas em funções interpoladas
f1, range1 = histogram_to_function(data1, nbins=50)
f2, range2 = histogram_to_function(data2, nbins=50)
f3, range3 = histogram_to_function(data3, nbins=50)

# --- Define intervalo comum para plotagem
x_min = minimum([range1[1], range2[1], range3[1]])
x_max = maximum([range1[2], range2[2], range3[2]])
x = range(x_min, x_max, length=1000)

# --- Cores em tons de cinza
gray1 = RGBA(0.2, 0.2, 0.2, 1.0)
gray2 = RGBA(0.6, 0.6, 0.6, 1.0)
gray3 = RGBA(0.8, 0.8, 0.8, 1.0)

# --- Plot final
plot(legendfontsize=14,
     guidefontsize=12,   # labels dos eixos
     tickfontsize=10,    # números dos ticks
     titlefontsize=16)   # título
plot!(x, f1.(x), label="Radius 0.1", linewidth=3, color=gray1)
plot!(x, f2.(x), label="Radius 0.2", linewidth=3, color=gray2)
plot!(x, f3.(x), label="Radius 0.3", linewidth=3, color=gray3)

# --- Eixos, título, legenda
xlabel!("eV")
ylabel!("Density")
title!("Energy Density Function")
plot!(legend=:topright,legendfontsize=12,guidefontsize=14)

display(current())
savefig("Plot.png")