using PrettyTables

function prettyprint(a, cnames="", rnames="")
if rnames !=""
    rnames = rnames[:]
    a = [rnames a]
    if cnames != ""
        cnames = cnames[:]
        cnames = vcat("", cnames)
    end    
end
if cnames !=""
    pretty_table(a, cnames; formatters=ft_printf("%12.5f"))
else
    pretty_table(a; formatters=ft_printf("%12.5f"))
end
end

#=
#using Printf
# formatted print of array, with column names
function prettyprint(a, cnames="", rnames="",digits=12, decimals=5)
    # TBD: try to use this to allow using specified digits and decimals
    #fmt = @sprintf("%d",digits)"."@sprintf("%d",decimals)"%f"
    #@eval dofmt(x) = @sprintf($fmt, x)
    
    # print column names
    if cnames != ""
        for i = 1:size(a,2)
            pad = digits
            if rnames != "" && i==1
                pad = 2*digits
            end    
            @printf("%s", lpad(cnames[i],pad," "))
        end
        @printf("\n")
    end    
    # print the rows
    for i = 1:size(a,1)
        if rnames != ""
            @printf("%s", lpad(rnames[i],digits," "))
        end
        for j = 1:size(a,2)
            # TBD: use fmt defined above to print array contents
            @printf("%12.5f",(a[i,j]))
        end
        @printf("\n")
    end
    return
end
=#
