library(learnr)
library(shiny)

# Define the path to the Rmd tutorial file
tutorial_file <- "Tutorial.Rmd"

# Launch the tutorial
rmarkdown::run(tutorial_file)

