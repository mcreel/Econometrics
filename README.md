# Econometrics
First year graduate level econometrics notes with embedded examples using the Julia language.

The notes cover linear regression models in the first half, and move on to maximum likelihood and GMM estimation of potentially nonlinear models in the second half. After these core methods, there are several chapters which introduce more specialized topics such as panel data, quantile regression, nonparametric method, etc. 

To simply view the notes, see the econometrics.pdf document. However, for the links to code that are in the pdf to work, you should download the whole repository, as is described below. Much of the material is presented with the expectation that the code will be examined, run, and experimented with.

NEWS: there is a MakeSysimage.jl script that will precompile all of the packages and  some of the commonly used functions. If you use VS Code you need to edit the Julia extension settings to enable “use an existing custom sysimage when starting the REPL”.


## To make full use of the notes and examples
You need to install this repository as a Julia project. Do this as follows:

1. Download a zip file, or git clone this repo, putting the contents at some convenient place on your storage media

2. open a terminal, and change directories to the location of the Econometrics code
3. start Julia, using ```julia --proj``` 
4. do ```] instantiate```  This will install a last few packages.

5. (optional, only if you're curious) basic testing of the package has been added. Do ```using Pkg; Pkg.test("Econometrics")```

6. then do ```using Econometrics``` in Julia to use the package. You can see some examples by typing 
   ```julia
   ols()
   mleresults()
   gmmresults()
   mcmc()
   npreg()
   samin()
   ```
  Running the script ```warmpup.jl``` will run these and other code bits, which will pre-compile many of the functions that are used in the examples. This will make the examples run more quickly. It is best to run the warmup script, and then keep the Julia session open, to run the examples when you like.

7. To see the source code for those examples, type ```edit(ols,())```, to see the OLS source code.

8. To run examples, cd into the relevant subdirectory of Econometrics/Examples, and then just include the script you would like to run.

## There are a couple of unusual thing about these notes:
- they are available in editable form (econometrics.lyx), so that you can modify them to suit your needs: see the first chapter for more information, and get LyX from  www.lyx.org. 
- they contain links that point to example programs using the Julia language. The examples show how to use the methods and illustrate properties of estimators. The example code can be modified to allow exploration.

## To integrate the notes and the example files, so that links are ready to run
- the file econometrics_local.lyx allows you to create a pdf version with links that open using the example files installed on your computer. See the search and replace information at the top of the document. If you configure your system to open .jl files in your favorite editor, then the links in the pdf will open with that editor. For example, on Linux, I set KDE Plasma file associations for .jl files to open using ```code -a ~/Mystuff/Econometrics -r %F```, where code is the VScode editor executable, and ~/Mystuff/Econometrics is the location where I have the archive installed. When I click a link in econometrics_local.pdf, the link opens in the editor, with the appropriate environment, ready to run.
To get an idea of how this works, the following figure, on the R, shows an explanation in the pdf version of the notes, with a link to an example. The code of the example is visible in the upper L, and the output of running the example in Julia is at the lower L. The plot appears in the VScode plot frame.

![example](https://github.com/mcreel/Econometrics/blob/master/example.png)
