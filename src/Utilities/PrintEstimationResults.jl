using PrettyTables
function PrintEstimationResults(results, names)
    results = [names[:] results]
    clabels = ["parameter", "estimate", "st. err", "t-stat", "p-value"]
    h1 = Highlighter(
                      (results, j) -> results[j,5] .< 0.01,
                      foreground = :green
                     );
    h5 = Highlighter(
                      (results, j) -> (results[j,5] .< 0.05) & (results[j,4] .> 0.01) ,
                      foreground = :blue
                     );
    h10 = Highlighter(
                       (results, j) -> (results[j,5] .< 0.1) & (results[j,4] .> 0.05) ,
                      foreground = :yellow
                     );
    pretty_table(results, clabels; highlighters=(h1, h5, h10))
end
