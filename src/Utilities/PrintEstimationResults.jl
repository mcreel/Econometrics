using PrettyTables
function PrintEstimationResults(results, names)
    results = [names[:] results]
    clabels = ["parameter", "estimate", "st. err", "t-stat", "p-value"]
    h10 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.1), crayon"yellow");
    h5 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.05), crayon"blue");
    h1 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.01), crayon"green");
    pretty_table(results, clabels; formatters=ft_printf("%12.5f"), highlighters=(h1, h5, h10))
end
