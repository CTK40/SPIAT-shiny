# SPIAT-shiny

This is a Shiny App that assists formatting spatial proteomics/transciptomics data into an SPE (SpatialExperiment) object. Accepted data format includes inForm and HALO for spatial proteomics data and Visium, Xenium, MERSCOPE, and CosMX for spatial transcriptomics data. All data format that was not listed here can be read via "general" format. 

To use the Shiny App, navigate to your local folder and use `git` to clone the `SPIAT-shiny` repo.
```
git clone https://github.com/TrigosTeam/SPIAT-shiny.git
```

If using RStudio, open the `app.R` file and click `Run App` to start the application;
If using command line, make sure you are in the SPIAT-shiny folder and enter `Rscript app.R`, and then paste the url into a browser to use the application. You can follow the video [here](https://youtu.be/gBYkKk2fuA4) to get started.

More video tutorials:
[Formatting inForm data](https://youtu.be/0VVJ9mWuXZY)
[Formatting HALO data](https://youtu.be/eb_hw5u4nRA)
[Formatting data from csv files ("general" option)](https://youtu.be/x45AVHv36aU)
