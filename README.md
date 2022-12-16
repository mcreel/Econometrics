[![testing](https://github.com/mcreel/Econometrics/actions/workflows/testing.yml/badge.svg)](https://github.com/mcreel/Econometrics/actions/workflows/testing.yml) (testing is done on latest stable release of Julia, for Linux, Windows, and MacOS)

# Econometrics
First year graduate level econometrics notes with embedded examples using the Julia language.

The notes cover linear regression models in the first half (about 30 hours of class time). The second half (another 30-40 hours of class time) moves on to nonlinear optimization, maximum likelihood and GMM estimation of potentially nonlinear models. After these core methods, there are several chapters which introduce more specialized topics such as panel data, quantile regression, nonparametric methods, etc.

The notes are in the file econometrics.pdf. It is not recommended to use this file by itself, as one will not have access to the linked code and examples. The notes were prepared with the expectation that the example code will be examined, run, and experimented with. See below for how to make full use of the materials.

## To make full use of the notes and examples

### How to videos (ctrl-click for new tab)
[The playlist with all videos](https://www.youtube.com/watch?v=I5lclHKlQnA&list=PLA7k_EnUgZs2U2Wh0bxhcabfKrGQQOCyn)

The basic introductory videos:
1. [Introduction](https://youtu.be/I5lclHKlQnA): short description of what's here
2. [Installation](https://youtu.be/N_aWT7OiX4k): how to install this as a Julia project
3. [Speed](https://youtu.be/a-_ZNTBeLCw): how to compile a system image to make it fast
4. [Comfort](https://youtu.be/Nbhmq4VWVJU): how to set up VS Code for comfortable use

Practical summary videos:
1. [Ch. 12, Numerical Optimization](https://www.youtube.com/watch?v=u7XQOusjZ3c): explains
   how to use the practical summaries, and does the summary for Ch. 12.
### Old School Instructions:
You need to install this repository as a Julia project. Do this as follows: 

1. download the code:
    (a) download a zip of the repo, and uncompress it in a convenient directory, or
    (b) git clone the repository to the desired location

2. Go to that directory and start Julia using ``julia --proj``

3. In Julia, the first time you use the files, do ``using Pkg; Pkg.instantiate()`` This will take some time, as Econometrics relies on a number of other packages.

4. then do ``using Econometrics`` in Julia to use the package. The first time you do this, it will take a **long** time, maybe 15 minutes or so. *Don't worry*, this is normal. All of the packages that were downloaded are being compiled for the first time. We will be able to make this go *much, much* faster when we want to use the code.

5. To run examples, cd into the relevant subdirectory of Econometrics/Examples, and then just include the script you would like to run.

6. Once this is done, you can use the code at any time by repeating steps 2 and 4.

7. You can see some examples by typing 
   ```julia
   ols()
   mleresults()
   gmmresults()
   mcmc()
   npreg()
   samin()
   ```
   when in Julia.

8. To run examples, cd into the relevant subdirectory of Econometrics/Examples, and then just include the script you would like to run.

9. I recommend setting your operating system to open .jl files with your favorite editor,
   to make the links in econometrics.pdf open in that editor.
10. To speed things up **a lot**, type ``ìnclude("MakeSysimage.jl")`` from Julia, having started Julia with ``julia --proj`` in the directory where you installed the notes. This will compile all of the packages, and keep the compiled images for re-use. If you do this, when you start Julia in the future, use ``julia --proj -J JuliaSysimage.so`` (*caveat*: I have only tried this using Linux. If you're a Windows or OSX user, please open an issue if these instructions don't work) 

## There are a couple of unusual thing about these notes:
- they are available in editable form (econometrics.lyx), so that you can modify them to suit your needs: see the first chapter for more information, and get LyX from  www.lyx.org. 
- they contain links that point to example programs using the Julia language. The examples show how to use the methods and illustrate properties of estimators. The example code can be modified to allow exploration.

To get an idea of how this works, the following figure, on the R, shows an explanation in the pdf version of the notes, with a link to an example. The code of the example is visible in the upper L, and the output of running the example in Julia is at the lower L. The plot appears in the VScode plot frame.

![example](https://github.com/mcreel/Econometrics/blob/main/example.png)
