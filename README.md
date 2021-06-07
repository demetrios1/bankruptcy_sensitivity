# bankruptcy_sensitivity

This is the github for the paper **Do forecasts of bankruptcy cause bankruptcy? A machine learning sensitivity analysis.**.  In the tex file is the latex code to produce the main file, figures are all the figures made for the project, codes includes R scripts to create the plots and tables, and the monotone_bart folder includes a tar file called ```fastbart_2.0.tar.gz``` that includes a modified version of BART that will need to be installed in order to run the other R scripts.  

In particular, the script ```methods_demo_script.R``` will provide a quick demonstration of our methodology.  In order to run it, you must first uncomment the line 
```#install.packages(path_to_file, repos = NULL, type="source")```, where the file of interest is in the monotone_bart folder and is the tar file ```fastbart_2.0.tar.gz```
