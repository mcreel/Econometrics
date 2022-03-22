using PrettyTables, Term
function PrintEstimationResults(results, names, title="")
    results = [names[:] results]
    clabels = ["parameter", "estimate", "st. err", "t-stat", "p-value"]
    h10 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.1), crayon"yellow");
    h5 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.05), crayon"blue");
    h1 = Highlighter((results,i,j) -> (j==5) && (results[i,5] < 0.01), crayon"green");
    Panel(
        pretty_table(results; header=clabels, formatters=ft_printf("%12.5f"), highlighters=(h1, h5, h10)),
    width=80, justify=:center, style="blue", box=:DOUBLE, title=title, title_style="white")
end
