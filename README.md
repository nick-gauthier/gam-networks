---
title: "Generalized additive mixed models for spatial networks"
author: "Nick Gauthier"
output: 
  html_document: 
    highlight: zenburn
    keep_md: yes
    theme: cosmo
    toc: yes
---


GAMMs are a flexible form of regression model well-matched to the complexities of the archaeological record, including non-normal distributions in the form of counts or proportions, non-independent observations with correlated errors, and non-linear functional relationships. Using two case studies -- an ethnographic marriage network and an archaeological assemblage similarity network -- I illustrate how this approach can lead to unbiased parameter estimates and more robust comparisons of competing hypotheses.

Distance is a fundamental constraint on human social interaction. This basic principle motivates the use of spatial interaction models for estimating flows of people, information, and resources on spatial and social networks. These models have both valid dynamical and statistical interpretations, a key strength well supported by theory and data from geography, economics, ecology, and genetics. To date, archaeologists have primarily relied on the dynamical approach because the idiosyncrasies of archaeological data make the wholesale adoption of statistical approaches from other fields impractical.

# Setup


```r
library(raster) # DEM for map
library(tidyverse) # data processing and plotting functions
library(sf) # tidy spatial data and mapping
library(mgcv) # fit GAMs and GAMMs
library(archdata) # source for the Oxford Pots data
library(tidygraph) # tidy network processing
# devtools::install_github('thomasp85/ggraph') # need dev version of ggraph 
library(ggraph) # network plots
# devtools::install_github('nspope/corMLPE')
library(corMLPE) # network correlation structure for undirected networks
# library(ggmap) # for geocoding, but not used by default (see below)
library(ggspatial) # scale bar for maps
library(maps) # country outlines for maps
```

# Case Studies

## Oxford Pots
First, we'll look at a dataset of Late Romano-British pottery available in the `archdata` package. \emph{\small{David L. Carlson and Georg Roth (2018). archdata: Example Datasets from Archaeological Research. R package version 1.2.}}

```r
data("OxfordPots")
```

Use `tidyr` and `dplyr` to get the data into a tidy format.

```r
pots <- OxfordPots %>%
  as_tibble %>%
  rename(to = Place) %>%
  gather(key, value, OxfordPct:NewForestDst) %>%
  separate(key, c('from', 'measure'), sep = -3) %>%
  spread(measure, value) %>%
  rename(percent = Pct, distance = Dst) %>%
  mutate(from = str_replace(from, 'NewForest', 'New Forest'),
         water = as.factor(WaterTrans == 1)) %>%
  select(to, from, percent, distance, water) %>%
  filter(!is.na(percent))
```


```
## # A tibble: 43 x 5
##    to               from       percent distance water
##    <chr>            <chr>        <dbl>    <dbl> <fct>
##  1 Alchester        Oxford       22.5        11 TRUE 
##  2 Bath             New Forest    1.25       45 TRUE 
##  3 Bath             Oxford       21.2        55 TRUE 
##  4 Brough-on-Humber Oxford        1.5       140 FALSE
##  5 Caerwent         Oxford       19          74 TRUE 
##  6 Caister          Oxford        4         135 FALSE
##  7 Canterbury       Oxford       17.5       105 TRUE 
##  8 Chichester       New Forest   13          40 FALSE
##  9 Chichester       Oxford        6.75       65 FALSE
## 10 Cholchester      Oxford        7          95 FALSE
## # … with 33 more rows
```

These data don't come with any location information, so we have to do some geocoding using the `ggmap` package. This requires you to get and register your own Google Maps API key, so by default this code chunk is not run and we instead use precomputed values. If you do want to run this chunk, follow the instructions in `?ggmap::register_google` and replace the line `<your-api-key-here>` with your API key.


```r
ggmap::register_google(key = "<your-api-key-here>")

pots_sites <- pots %>%
  select(from, to) %>% 
  gather %>%
  pull(value) %>%
  unique %>% 
  tibble(site = .) %>%
  # adjust some names 
  mutate(search_term = str_replace_all(site, c('Clausentum' = 'Bitterne', 
                                               'Mildenhall' = 'Mildenhall, Wiltshire')),
         search_term = paste0(search_term, ', UK')) %>%
  mutate_geocode(search_term) %>%
  select(-search_term) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  bind_cols(., as_tibble(st_coordinates(.))) %>%
  rename(x = X, y = Y)

saveRDS(pots_sites, 'oxfordpots_locations.RDS')
```

Join the Oxford Pots data to the geocoded locations.

```r
pots_sites <- readRDS('data/oxfordpots_locations.RDS') %>%
  left_join(dplyr::select(pots, to), by = c('site' = 'to'))
```

Make a network object using ``tidygraph` for plotting, and an edgelist tibble for modeling.

```r
pots_net <- tbl_graph(nodes = pots_sites, edges = pots, directed = TRUE) %E>%
  filter(!is.na(percent)) %>%
  mutate(similarity = percent / 100)

pots_dat <- pots_net %E>%
  mutate(x = .N()$x[to],
         y = .N()$y[to]) %>% 
  as_tibble %>%
  mutate(from = as.factor(from),
         to = as.factor(to))

pots_net
```

```
## # A tbl_graph: 45 nodes and 43 edges
## #
## # A directed acyclic simple graph with 14 components
## #
## # Edge Data: 43 x 6 (active)
##    from    to percent distance water similarity
##   <int> <int>   <dbl>    <dbl> <fct>      <dbl>
## 1     2     3   22.5        11 TRUE      0.225 
## 2     1     4    1.25       45 TRUE      0.0125
## 3     2     4   21.2        55 TRUE      0.212 
## 4     2     6    1.5       140 FALSE     0.015 
## 5     2     7   19          74 TRUE      0.19  
## 6     2     8    4         135 FALSE     0.04  
## # … with 37 more rows
## #
## # Node Data: 45 x 4
##   site           x     y             geometry
##   <chr>      <dbl> <dbl>          <POINT [°]>
## 1 New Forest -1.63  50.9 (-1.631463 50.87652)
## 2 Oxford     -1.26  51.8 (-1.257726 51.75202)
## 3 Alchester  -1.87  52.2 (-1.867605 52.21531)
## # … with 42 more rows
```




```r
knitr::knit_exit()
```
















































