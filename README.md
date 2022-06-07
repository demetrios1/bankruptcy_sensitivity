# bankruptcy_sensitivity

This is the github for the paper **Do forecasts of bankruptcy cause bankruptcy? A machine learning sensitivity analysis.**, which can be accessed additionally on [Arxiv](https://arxiv.org/pdf/2106.04503.pdf). The paper pdf file is on the main page of the repo and in the tex file is the latex code to produce the main file as well as a supplemental file with some additional results.  The figures folder contains all the figures made for the project, codes includes R scripts to create the plots and tables.

How to get started.  
The first thing you want to  do is 

```devtools::install_github("jaredsmurray/monbart", ref='main')```
This will install one of the necessary packages from [Jared Murray](https://jaredsmurray.github.io/), called [monbart](https://github.com/jaredsmurray/monbart).  This implements the BART with monotonicity constraint described in the paper.

In particular, the script ```methods_demo_script.R``` will provide a quick demonstration of our methodology.  In order to run it, you must first download [monbart](https://github.com/jaredsmurray/monbart/tree/main/src), which provides a monotone bart implementation.
