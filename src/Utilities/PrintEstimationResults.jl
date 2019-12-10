using PrettyTables
function PrintEstimationResults(results, names)
    results = [names[:] results]
    clabels = ["parameter", "estimate", "st. err", "t-stat", "p-value"]
    pretty_table(results, clabels; formatter=ft_printf("%12.5f"))
end
