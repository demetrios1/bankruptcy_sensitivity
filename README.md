# bankruptcy_sensitivity

This is the github for the paper **Do forecasts of bankruptcy cause bankruptcy? A machine learning sensitivity analysis.**, which can be accessed additionally on [Arxiv](https://arxiv.org/pdf/2106.04503.pdf). The paper pdf file is on the main page of the repo and in the tex file is the latex code to produce the main file as well as a supplemental file with some additional results.  The figures folder contains all the figures made for the project, codes includes R scripts to create the plots and tables, and the monotone_bart folder includes a tar file called ```fastbart_2.0.tar.gz``` that includes a modified version of BART that will need to be installed in order to run the other R scripts.  

In particular, the script ```methods_demo_script.R``` will provide a quick demonstration of our methodology.  In order to run it, you must first download [monbart](https://github.com/jaredsmurray/monbart/tree/main/src), which provides a monotone bart implementation.
