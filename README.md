# *Scorecard Diplomacy* code repository

[Judith G. Kelley](https://sanford.duke.edu/people/faculty/kelley-judith) • [Sanford School of Public Policy](https://sanford.duke.edu/) • [Duke University](https://www.duke.edu/)

---

> Judith G. Kelley. 2017. *Scorecard Diplomacy: Grading States to Influence their Reputation and Behavior.* Cambridge: Cambridge University Press. (ISBN 978-1316649138)

This repository contains all the code and data for Judith Kelley's [*Scorecard Diplomacy*](http://www.cambridge.org/se/academic/subjects/politics-international-relations/international-relations-and-international-organisations/scorecard-diplomacy-grading-states-influence-their-reputation-and-behavior) book project. 

The commit history for each file is truncated, however. Heavy development on the code occurred in a [separate repository](https://github.com/andrewheiss/jk_misc) between August 2015 and September 2016, but that repository was poorly structured, could not be run cleanly, and contained code for other projects. This repository was created as a simpler, fully reproducible, standalone repository. Those interested in the full development history of each file can visit the [`jk_misc`](https://github.com/andrewheiss/jk_misc) repository.

The book also uses code and data two other repositories (though these are not needed to generate the figures and tables for the book):

- [andrewheiss/From-the-Trenches-Anti-TIP-NGOs-and-US](https://github.com/andrewheiss/From-the-Trenches-Anti-TIP-NGOs-and-US): Repository for Andrew Heiss and Judith G. Kelley. 2016. "From the Trenches: A Global Survey of Anti-TIP NGOs and their Views of US Efforts." *Journal of Human Trafficking*. doi: 10.1080/23322705.2016.1199241
- [andrewheiss/Human-trafficking-NGO-survey](https://github.com/andrewheiss/Human-trafficking-NGO-survey): Repository for cleaning up the NGO survey

## Prerequisites

All the results, figures, and tables in the book can be recreated using [R 3.3](https://www.r-project.org/) (preferably within [RStudio](https://www.rstudio.com/), since the code is structured as an [RStudio project](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects)). 

You'll also need to install the following R packages from CRAN:

- `tidyverse`: Hadley Wickham's universe of tidy-friendly packages (automatically installs `dplyr`, `tidyr`, `ggplot2`, and friends)
- `ggstance`: Adds horizontal `geom`s for `ggplot2` 
- `ggrepel`: Automatically repel overlapping labels in `ggplot2` plots
- `gridExtra`: Layout multiple `ggplot2` plots
- `pander`: A convenient interface for [Pandoc](http://pandoc.org/)
- `stargazer`: Create pretty regression tables
- `Cairo`: Nicer graphics creation library that correctly embeds fonts into PDF plots and correctly generates high resolution PNG plots
- `countrycode`: Convert between country names and ISO/COW codes
- `WDI`: Access the World Bank's [World Development Indicator API](http://data.worldbank.org/)
- `maptools`: Read spatial objects
- `rgdal`: Use spatial objects with [gdal](http://www.gdal.org/)
- `testthat`: Perform unit tests
- `zoo`: Various time series-related functions
