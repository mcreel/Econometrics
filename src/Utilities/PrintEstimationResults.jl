using PrettyTables
function PrintEstimationResults(results, names)
    clabels = ["estimate", "st. err", "t-stat", "p-value"]
    h1 = Highlighter(
                      (results, j) -> results[j,4] .< 0.01,
                      foreground = :green
                     );
    h5 = Highlighter(
                      (results, j) -> (results[j,4] .< 0.05) & (results[j,4] .> 0.01) ,
                      foreground = :blue
                     );
    h10 = Highlighter(
                       (results, j) -> (results[j,4] .< 0.1) & (results[j,4] .> 0.05) ,
                      foreground = :yellow
                     );
    prettyprint(results, clabels, names; highlighters=(h1, h5, h10))
end
