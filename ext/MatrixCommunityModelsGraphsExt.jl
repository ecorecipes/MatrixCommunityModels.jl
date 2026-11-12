module MatrixCommunityModelsGraphsExt

using MatrixCommunityModels: InteractionNetwork
import Graphs

"""
    Graphs.SimpleDiGraph(net::InteractionNetwork)

The links of an `InteractionNetwork` as a Graphs.jl directed graph, with vertices in
label order and an edge `j → i` where node `j` affects node `i`. Use it for the network
measures Graphs.jl provides (centrality, components, cycles, ...).
"""
function Graphs.SimpleDiGraph(net::InteractionNetwork)
    g = Graphs.SimpleDiGraph(length(net.labels))
    for j in axes(net.adjacency, 2), i in axes(net.adjacency, 1)
        net.adjacency[i, j] && Graphs.add_edge!(g, j, i)
    end
    return g
end

end
