## Resubmission

This is a resubmission. Thanks very much for the comments from Gregor Seyer. The issues pointed out are fixed. The responses are as follows:

>Please always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'worldclim'

As a dataset, WorldClim and CMCC-BioClimInd are written in single quotes.

>If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
  authors (year) <arXiv:...>
  authors (year, ISBN:...)
or if those are not available: <https:...>
  with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

All necessary references are added in DESCRIPTION. The narrative is modified a bit correspondingly.

>Please make sure that you do not change the user's options, par or 
working directory. If you really have to do so within functions, please 
ensure with an *immediate* call of on.exit() that the settings are reset 
when the function is exited. e.g.:
...
old <- options()         # code line i
on.exit(options(old))     # code line i+1
...
options(timeout = 1e5)
...
e.g.: cmcc_bioclim.R, future_cmcc_bioclim.R, future_worldclim2.R, 
worldclim2.R
If you're not familiar with the function, please check ?on.exit. This 
function makes it possible to restore options before exiting a function 
even if the function breaks. Therefore it needs to be called immediately 
after the option change within a function.

Oops, this is a rookie mistake. Deleted these settings from functions. It is good to know the safe way to change user's options within package functions.

## Test environments

1. Local macOS Monterey 12.1, R version 4.0.2

2. Github actions

- Microsoft Windows Server 2019 10.0.17763, R version 4.1.2 (2021-11-01)
- Ubuntu 20.04.3 LTS, R version 4.1.2 (2021-11-01)
- Ubuntu 20.04.3 LTS, R Under development (unstable) (2022-01-01 r81419)
- Ubuntu 20.04.3 LTS, R version 4.0.5 (2021-03-31)

3. rhub/devtools check

- Windows Server 2022, R-devel, 64 bit
- Fedora Linux, R-devel, clang, gfortran
- Ubuntu Linux 20.04.1 LTS, R-release, GCC
- x86_64-w64-mingw32 (64-bit), R version 4.0.5 (2021-03-31)
- x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2022-01-03 r81439 ucrt)
- x86_64-w64-mingw32 (64-bit), R version 4.1.2 (2021-11-01)

## R CMD check results

There is no error or warning. It got one note related to CRAN release:

```
New submission

Possibly misspelled words in DESCRIPTION:
  Brunner (12:9)
  Caporaso (28:38)
  Fick (26:66)
  Guha (14:5)
  Hariri (11:51)
  Hijmans (27:14)
  Kononenko (21:56)
  Liu (10:45, 12:63)
  Lundberg (18:60)
  Mishra (14:15)
  Noce (28:28)
  SCiForest (9:54)
  Santini (28:56)
  Schrijvers (14:39)
  Shapley (16:52, 18:5)
  Zhou (11:5, 13:20)
  bioclimatic (26:5)
  iForest (8:50)
  itsdm (22:45)
  trumbelj (21:39)
```

- This is a new release.
- These possibly misspelled words are not misspelled. 

## Others notes

On rhub check: Windows Server 2022, R-devel, 64 bit, it got an extra note:

```
Found the following files/directories:
'lastMiKTeXException'
```
It seems related to knit. I can't reproduce this note anywhere else.

I pre-compiled the package vignettes to save checking time. They already be tested on multiple platforms.
