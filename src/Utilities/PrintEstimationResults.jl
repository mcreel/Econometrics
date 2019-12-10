using PrettyTables
function PrintEstimationResults(results, names)
    results = [names[:] results]
    clabels = ["parameter", "estimate", "st. err", "t-stat", "p-value"]
    pretty_table(results, clabels)
end
