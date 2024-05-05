
* ggstatsplot

* t-test / anova，non-parametric， corr， regression

* Central tendency measure
Parametric	mean 
Non-parametric	median 
Robust	trimmed mean 
Bayesian	MAP (maximum a posteriori probability) estimate 


* https://indrajeetpatil.github.io/ggstatsplot/index.html


#install.packages("ggstatsplot")
library(ggstatsplot)


# I. ggbetweenstats

* disproportionate: Welch's ANOVA

# for reproducibility
set.seed(123)
head(iris)

# plot
ggbetweenstats(
  data = iris,
  x = Species,
  y = Sepal.Length,
  title = "Distribution of sepal length across Iris species"
)


# II. ggwithinstats

* similar to ggbetweenstats
* pair-wise test

# for reproducibility and data
set.seed(123)
library(WRS2) # for data
library(afex) # to run anova
#BiocManager::install('ggthemes')
head(WineTasting)
table(WineTasting$Taster) #22对
# plot
ggwithinstats(
  data = WineTasting,
  x = Wine,
  y = Taste,
  title = "Wine tasting",
  caption = "Data source: `WRS2` R package",
  ggtheme = ggthemes::theme_fivethirtyeight()
)


# III. gghistostats 

* test.value 

# for reproducibility
set.seed(123)

# plot
gghistostats(
  data = ggplot2::msleep, # dataframe from which variable is to be taken
  x = awake, # numeric variable whose distribution is of interest
  title = "Amount of time spent awake", # title for the plot
  caption = substitute(paste(italic("Source: "), "Mammalian sleep data set")),
  test.value = 12, # default value is 0
  normal.curve = TRUE, # 
  normal.curve.args = list(color = "red", size = 1),
  binwidth = 1, # 
  ggtheme = hrbrthemes::theme_ipsum_tw(),
  bf.message = FALSE ## 
)


# IV. ggdotplotstats 
* type = "robust" :Bootstrap-t method for one-sample test

# plot
ggdotplotstats(
  data = dplyr::filter(.data = gapminder::gapminder, continent == "Asia"),
  y = country,
  x = lifeExp,
  test.value = 55,
  type = "robust",
  title = "Distribution of life expectancy in Asian continent",
  xlab = "Life expectancy",
  caption = substitute(
    paste(
      italic("Source"),
      ": Gapminder dataset from https://www.gapminder.org/"
    )
  )
)

# V. ggscatterstats 

1. marginal.type：
  * histograms
  * boxplots
  * density
  * violin
  * densigram (density + histogram)
  
2. test：
  * Parametric	Pearson’s correlation coefficient
  * Non-parametric	Spearman’s rank correlation coefficient
  * Robust	Winsorized Pearson correlation coefficient
  * Bayesian	Pearson’s correlation coefficient
  
3. remove bayes: bf.message = FALSE 


data <- dplyr::filter(movies_long, genre == "Action")
head(data)
# plot
ggscatterstats(
  data = data,
  x = budget,
  y = rating,
  type = "robust", # Winsorized Pearson correlation coefficient
  xlab = "Movie budget (in million/ US$)", # label for x axis
  ylab = "IMDB rating", # label for y axis
  label.var = title, 
  label.expression = rating < 5 & budget > 100, 
  title = "Movie budget and IMDB rating (action)", # title text for the plot
  caption = expression(paste(italic("Note"), ": IMDB stands for Internet Movie DataBase")),#右下方注释
  ggtheme = hrbrthemes::theme_ipsum_ps(), # choosing a different theme
  # turn off `ggstatsplot` theme layer
  marginal.type = "densigram",
  xfill = "pink", # color fill for x-axis marginal distribution
  yfill = "#009E73" # color fill for y-axis marginal distribution
)

# VI. ggcorrmat 

# for reproducibility
set.seed(123)

# as a default this function outputs a correlation matrix plot
ggcorrmat(
  data = ggplot2::msleep,
  type = "Parametric",
  colors = c("#B2182B", "white", "#4D4D4D"),
  title = "Correlalogram for mammals sleep dataset",
  subtitle = "sleep units: hours; weight units: kilograms"
)

# output = "dataframe" 
head(
  ggcorrmat(
  data = ggplot2::msleep,
  type = "Parametric",
  colors = c("#B2182B", "white", "#4D4D4D"),
  title = "Correlalogram for mammals sleep dataset",
  subtitle = "sleep units: hours; weight units: kilograms",
  output = "dataframe"
))

# VII. ggpiestats 


# for reproducibility
set.seed(123)

# plot
ggpiestats(
  data = mtcars,
  x = am,
  y = cyl,
  package = "wesanderson",
  palette = "Royal1", 
  title = "Dataset: Motor Trend Car Road Tests", # title for the plot
  legend.title = "Transmission", # title for the legend
  caption = substitute(paste(italic("Source"), ": 1974 Motor Trend US magazine"))
)


# VIII. ggbarstats 

# plot
ggbarstats(
  data = movies_long,
  x = mpaa,
  y = genre,
  title = "MPAA Ratings by Genre",
  xlab = "movie genre",
  legend.title = "MPAA rating",
  ggtheme = hrbrthemes::theme_ipsum_pub(),
  palette = "Set2"
)


