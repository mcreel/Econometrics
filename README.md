# Econometrics
Graduate level econometrics notes with embedded examples using the Julia language.

NEWS: there is a MakeSysimage.jl script that will precompile all of the packages and 
some of the commonly used functions. Works from REPL, but not from VS Code. For VS Code,
use its system for making a system image.

To get just the notes, click on econometrics.pdf, and then on Download, at the upper R of the page, which will download only the pdf. Links in the pdf point to files here on github, and will open in your browser.

## To run the examples
You need to install this repository as a Julia package. Do this as follows:

1. git clone this repo to some convenient place on your storage media

2. add the package to Julia by doing ```] add <path where you downloaded the repo>```, for example, if you cloned it into the git directory of your home directory, you would do ```] add ~/git/Econometrics```  This will install many supporting packages, so be patient.

3. go to the location of the Econometrics code, start Julia, and do ```] activate .; instantiate```  This will install a last few packages.

4. (optional, only if you're curious) basic testing of the package has been added. Do ```using Pkg; Pkg.test("Econometrics")```

5. then do ```using Econometrics``` in Julia to use the package. You can see some examples by typing 
   ```julia
   ols()
   mleresults()
   gmmresults()
   mcmc()
   npreg()
   samin()
   ```
  Running the script ```warmpup.jl``` will run these and other code bits, which will pre-compile many of the functions that are used in the examples. This will make the examples run more quickly. It is best to run the warmup script, and then keep the Julia session open, to run the examples when you like.

6. To see the source code for those examples, type ```edit(ols,())```, to see the OLS source code.

7. To run examples, cd into the relevant subdirectory of Econometrics/Examples, and then just include the script you would like to run.

## There are a couple of unusual thing about these notes:
- they are available in editable form (econometrics.lyx), so that you can modify them to suit your needs: see the first chapter for more information, and get LyX from  www.lyx.org. 
- they contain links that point to example programs using the Julia language. The examples show how to use the methods and illustrate properties of estimators. The example code can be modified to allow exploration.

## To integrate the notes and the example files, so that links are ready to run
- the file econometrics_local.lyx allows you to create a pdf version with links that open using the example files installed on your computer. See the search and replace information at the top of the document. If you configure your system to open .jl files in your favorite editor, then the links in the pdf will open with that editor. For example, on Linux, I set KDE Plasma file associations for .jl files to open using ```code -a ~/Mystuff/Econometrics -r %F```, where code is the VScode editor executable, and ~/Mystuff/Econometrics is the location where I have the archive installed. When I click a link in econometrics_local.pdf, the link opens in the editor, with the appropriate environment, ready to run.
To get an idea of how this works, the following figure, on the R, shows an explanation in the pdf version of the notes, with a link to an example. The code of the example is visible in the upper L, and the output of running the example in Julia is at the lower L. The plot appears in the VScode plot frame.

![example](https://github.com/mcreel/Econometrics/blob/master/example.png)
