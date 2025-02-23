# ---------------------------------------- #
# Cumulative SARS-CoV-2 in Georgia, U.S.A.
# ---------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: April 11, 2020
#
# Recently modified by: @idblr
# Recently modified on: April 13, 2020
#
# Notes:
# A) 04/11/2020 (IB) - Basic data importation, management, and visualization of SARS-CoV-2 data from Johns Hopkins University
# B) 04/11/2020 (IB) - Data visualization of cumulative SARS-CoV-2 cases in Georgia January 22, 2020 to April 10, 2020 (raw and per capita)
# C) 04/11/2020 (IB) - Data visualizations are static, greyscale, and set at sextiles (not to be widely disseminated)
# D) 04/11/2020 (IB) - Population from 2010 census and GIS shapefile of GA counties from (must use your own key)
# E) 04/11/2020 (IB) - Added code for interactive leaflet plot (not to be widely disseminated)
# F) 04/13/2020 (IB) - Cleaned up code and streamlined data management
# G) 04/13/2020 (IB) - Added calculation for percent non-caucasian householders from 2010 Census Short Form
# ---------------------------------------- #

####################
# SETTINGS & PATHS #
####################

# the URL of the Johns Hopkins University data page
# "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
jhu_url <- paste("https://raw.githubusercontent.com/CSSEGISandData/", 
                 "COVID-19/master/csse_covid_19_data/", 
                 "csse_covid_19_time_series/",
                 "time_series_covid19_confirmed_US.csv",
                 sep = ""
                 )

# U.S. Census API key
## Obtain one at http://api.census.gov/data/key_signup.html
census_key <- "XXXXX" # INSERT PERSONAL CENSUS KEY HERE

############
# PACKAGES #
############

library(broom)
library(dplyr)
library(ggplot2)
library(htmlwidgets)
library(leaflet)
library(readr)
library(sp)
library(tigris)
library(tidycensus)

####################
# DATA IMPORTATION #
####################

# SARS-CoV-2 daily case data in the United States
covid <- readr::read_csv(jhu_url)
names(covid)[10] <- "Lon" # fix typo in column named "Long_"

# County Shapefiles (Based on 2018 County Boundaries)
tigris_spdf <- tigris::counties()
#str(tigris_spdf)

# 2010 Population U.S. Cenus Bureau
#census2010 <- tidycensus::load_variables(2010, "sf1")
#census2010[20:30,]

GA_pop_2010 <- tidycensus::get_decennial(geography = "county",
                                         variables = "P001001", # total population
                                         year = 2010,
                                         state = c("Georgia"),
                                         geometry = TRUE,
                                         key = census_key
)
GA_pop_2010 <- tibble(GA_pop_2010)[,-c(2,3,5)] # drop geometry, variable code, and county name
names(GA_pop_2010)[2] <- "Pop2010"

GA_cauc_2010 <- tidycensus::get_decennial(geography = "county",
                                          variables = "H006002", # total householder caucasian-only
                                          year = 2010,
                                          state = c("Georgia"),
                                          geometry = TRUE,
                                          key = census_key
)
GA_cauc_2010 <- tibble(GA_cauc_2010)[,-c(2,3,5)] # drop geometry, variable code, and county name
names(GA_cauc_2010)[2] <- "Cauc2010"

GA_hous_2010 <- tidycensus::get_decennial(geography = "county",
                                          variables = "H001001", # total households
                                          year = 2010,
                                          state = c("Georgia"),
                                          geometry = TRUE,
                                          key = census_key
)
GA_hous_2010 <- tibble(GA_hous_2010)[,-c(2,3,5)] # drop geometry, variable code, and county name
names(GA_hous_2010)[2] <- "Hous2010"

# Merge 2010 Census data
GA_2010 <- merge(GA_pop_2010, GA_cauc_2010, by = "GEOID")
GA_2010 <- merge(GA_2010, GA_hous_2010, by = "GEOID")

# Calculate percent non-caucasian householders
GA_2010$pernoncauc <- (GA_2010$Hous2010 - GA_2010$Cauc2010)/GA_2010$Hous2010*100

GA_pop_2010 <- NULL # conserve memory
GA_cauc_2010 <- NULL # conserve memory
GA_hous_2010 <- NULL # conserve memory

####################
# DATA MANAGEMENT #
####################

# Fix FIPS in covid dataset to include leading 0 for states with STATE FIPS < 10
library(stringr)
covid$FIPS <- stringr::str_pad(covid$FIPS, 5, "left", "0")

# Reformat column names for dates to match world standard
names(covid)[12:ncol(covid)] <- format(as.Date(names(covid[12:ncol(covid)]), format = "%m/%d/%y"),
                                       format = "%d/%m/%y"
                                       )

# Merge case data with county shapefile
geoCoV <- sp::merge(tigris_spdf, covid, by.x = "GEOID", by.y = "FIPS", all.x = TRUE, all.y = TRUE)
tigris_spdf <- NULL # conserve memory, remove full county-level tigris data
covid <- NULL # conserve memory

# Subset for complete data
geoCoV_na <- geoCoV[is.na(geoCoV$UID),] # jurisdictions with missing case data are U.S. territories 
geoCoV <- geoCoV[!is.na(geoCoV$UID),] # omit jurisdictionss with missing case data

# Calculate cumulative cases
geoCoV$cumulative <- rowSums(geoCoV@data[29:ncol(geoCoV)], na.rm = TRUE) 

# Subset for Georiga
CoV_GA <- geoCoV[geoCoV@data$Province_State == "Georgia",]
geoCoV <- NULL # conserve memory

# Merge case data with population data
CoV_GA_pop <- sp::merge(CoV_GA, tibble(GA_2010), "GEOID", "GEOID")
CoV_GA <- NULL # conserve memory
GA_2010 <- NULL # conserve memory

# Calculate cumulative cases per capita
CoV_GA_pop$cumpercap <- CoV_GA_pop$cumulative/CoV_GA_pop$Pop2010

# Calculate cumulative cases per 100,000
CoV_GA_pop$cumrate <- CoV_GA_pop$cumpercap*100000

# Geographic Projection
## NAD83 / UTM zone 17N (Georgia, USA)
CoV_GA_UTM17N <- sp::spTransform(CoV_GA_pop, sp::CRS("+init=EPSG:26917")) 

######################
# DATA VISUALIZATION #
######################

# Cumulative SARS-CoV-2 Cases in Georgia (22/01/2020 - 04/10/2020)

f <- 1 # expansion factor
## The scale argument sets length of bar in map units
text1 <- list("sp.text", c(80000+10000, 3359000+15000), "0 km", cex = 1*f) # lower limit set just above scale
text2 <- list("sp.text", c(80000+110000, 3359000+15000), "100 km", cex = 1*f) # upper limit set just above scale

## Custom scale
scale <- list("SpatialPolygonsRescale",
              sp::layout.scale.bar(),
              scale = 100000, # 100 km
              fill = c("transparent", "black"), # colors of scale
              offset = c(80000+10000, 3359000+2000) # location in plot
              #offset = c(0,0) # location in plot
              )

## Custom compass rose
arrow <- list("SpatialPolygonsRescale", 
              sp::layout.north.arrow(),
              scale = 25000, # size of plot 
              offset = c(80000+150000, 3359000+1000) # location in plot (just above and middle of scale)
              )

## Plot of cumulative cases (unadjusted)
### Colorkey scaled by count per zipcode
at_break <- seq(from = 0, 
                to = max(CoV_GA_UTM17N$cumulative),
                by = max(CoV_GA_UTM17N$cumulative)/6
                ) # set breaks in colorkey (must be one more than default to fit min and max values)
at_names <- round(seq(from = 0, 
                      to = max(CoV_GA_UTM17N$cumulative),
                      by = max(CoV_GA_UTM17N$cumulative)/6
                      ),
                  digits = 0
                  )  # set name of breaks in colorkey
### Plot
grDevices::png(file = "figures/COVID_Georgia_Cumulative.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, # data
           "cumulative", # column name
           col.regions = gray.colors(length(at_break)), # color palette
           at = at_break, 
           par.settings = list(axis.line = list(col =  'transparent')), # remove default box around plot
           colorkey = list(labels = list(at = c(at_break[-length(at_break)],
                                                at_break[length(at_break)]-0.000000001
                                                ), # set top break a little less than default (to include max value on colorkey)
                                         labels = at_names, # set name of breaks to match data
                                         cex = 2*f, # size of labels of colorkey
                                         fontface = 1, # set font style (normal)
                                         legend = "frequency" # label of colorkey
                                         )
                           ),
           sub = list(label = "Copywrite Ian Buller", cex = 1*f), # subtitle for credit
           sp.layout = list(scale, text1, text2, arrow) # additions: scale and north arrow
           )
dev.off()

## Plot of cumulative per capita
### Colorkey scaled by count per zipcode
at_break <- seq(from = 0,
                to = max(CoV_GA_UTM17N$cumpercap),
                by = max(CoV_GA_UTM17N$cumpercap)/6
                ) # set breaks in colorkey (must be one more than default to fit min and max values)
at_names <- format(round(seq(from = 0,
                             to = max(CoV_GA_UTM17N$cumpercap),
                             by = max(CoV_GA_UTM17N$cumpercap)/6
                             ), digits = 2
                         ), nsmall = 1
                   )   # set name of breaks in colorkey
### Plot
grDevices::png(file = "figures/COVID_Georgia_Cumulative_percapita.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, # data
           "cumpercap", # column name
           col.regions = gray.colors(length(at_break)), # color palette
           at = at_break, 
           par.settings = list(axis.line = list(col =  'transparent')), # remove default box around plot
           colorkey = list(labels = list(at = c(at_break[-length(at_break)],
                                                at_break[length(at_break)]-0.000000001
                                                ), # set top break a little less than default (to include max value on colorkey)
                                         labels = at_names, # set name of breaks to match data
                                         cex = 2*f, # size of labels of colorkey
                                         fontface = 1, # set font style (normal)
                                         legend = "frequency" # label of colorkey
                                         )
                           ),
           sub = list(label = "Copywrite Ian Buller", cex = 1*f), # subtitle for credit
           sp.layout = list(scale, text1, text2, arrow) # additions: scale and north arrow
           )
dev.off()

## Using ggplot
### helpful material: https://cengel.github.io/rspatial/4_Mapping.nb.html
### Data preparation, ggplot2 requires a data.frame
CoV_GA_df <- broom::tidy(CoV_GA_UTM17N) # convert to tidy data frame
CoV_GA_UTM17N$polyID <- sapply(slot(CoV_GA_UTM17N, "polygons"), function(x) slot(x, "ID")) # preserve polygon id
CoV_GA_df <- merge(CoV_GA_df, CoV_GA_UTM17N, by.x = "id", by.y="polyID") # merge data

### Plot of cumulative cases per capita
f <- 1 # exansion factor
png(file = "figures/COVID_Georgia_Cumulative_percapita_ggplot.png", height = 1000*f, width = 1000*f)
ggplot2::ggplot() +                                      # initialize ggplot object
  ggplot2::geom_polygon(                                 # make a polygon
    data = CoV_GA_df,                                    # data frame
    ggplot2::aes(x = long, y = lat, group = group,       # coordinates, and group them by polygons
        fill = ggplot2::cut_number(cumpercap, 6)),       # variable to use for filling
    colour = "black") +                                  # color of polygon borders
  ggplot2::scale_fill_brewer("Cumulative cases per capita",# title of colorkey 
                             palette = "Greys",          # fill with brewer colors 
                             direction = -1,             # reverse colors in colorkey
                             guide = ggplot2::guide_legend(reverse = T)) +  # reverse order of colokey
  ggplot2::ggtitle("Cumulative SARS-CoV-2 cases per capita (January 22, 2020 - April 10, 2020)", # add title
                   subtitle = "Copywrite Ian Buller") +  # add subtitle
  ggplot2::theme(line = ggplot2::element_blank(),        # remove axis lines
        axis.text = ggplot2::element_blank(),            # remove tickmarks
        axis.title = ggplot2::element_blank(),           # remove axis labels
        panel.background = ggplot2::element_blank(),     # remove background gridlines
        text = ggplot2::element_text(size = 15*f)) +     # set font size
  ggplot2::coord_equal()                                 # both axes the same scale
dev.off()

### Plot of cumulative cases per 100,000
f <- 1 # exansion factor
png(file = "figures/COVID_Georgia_Cumulative_Rate_ggplot.png", height = 1000*f, width = 1000*f)
ggplot2::ggplot() +                                      # initialize ggplot object
  ggplot2::geom_polygon(                                 # make a polygon
    data = CoV_GA_df,                                    # data frame
    ggplot2::aes(x = long, y = lat, group = group,       # coordinates, and group them by polygons
                 fill = ggplot2::cut_number(cumrate, 6)),# variable to use for filling
    colour = "black") +                                  # color of polygon borders
  ggplot2::scale_fill_brewer("Cumulative rate",          # title of colorkey 
                             palette = "Greys",          # fill with brewer colors 
                             direction = -1,             # reverse colors in colorkey
                             guide = ggplot2::guide_legend(reverse = T)) +  # reverse order of colokey
  ggplot2::ggtitle("Cumulative SARS-CoV-2 rate per 100,000 (January 22, 2020 - April 10, 2020)", # add title
                   subtitle = "Copywrite Ian Buller") +  # add subtitle
  ggplot2::theme(line = ggplot2::element_blank(),        # remove axis lines
                 axis.text = ggplot2::element_blank(),            # remove tickmarks
                 axis.title = ggplot2::element_blank(),           # remove axis labels
                 panel.background = ggplot2::element_blank(),     # remove background gridlines
                 text = ggplot2::element_text(size = 15*f)) +     # set font size
  ggplot2::coord_equal()                                 # both axes the same scale
dev.off()

## Leaflet Plot
### Work with unprojected spatialpolygonsdataframe CoV_GA_pop
### Project to WGS84 EPSG:4326
CoV_GA_WGS84 <- sp::spTransform(CoV_GA_pop, CRS("+init=epsg:4326"))

### Popups
CoV_GA_WGS84$popup1 <- paste(CoV_GA_WGS84$NAME, "County:", format(round(CoV_GA_WGS84$cumulative, digits = 0), big.mark = ",", trim = T), "total cases", sep = " ")
CoV_GA_WGS84$popup2 <- paste(CoV_GA_WGS84$NAME, "County:", format(round(CoV_GA_WGS84$cumrate, digits = 0), big.mark = ",", trim = T), "per 100,000", sep = " ")

### Palettes
pal_cum <- leaflet::colorNumeric(palette = "Greys", domain = CoV_GA_WGS84$cumulative)
pal_rate <- leaflet::colorNumeric(palette = "Greys", domain = CoV_GA_WGS84$cumrate)

### Create leaflet plot
ga_m1 <- leaflet::leaflet(CoV_GA_WGS84) %>% # initial data = Georiga county-level COVID-19 data
  leaflet::setView(lng = -83, lat = 32, zoom = 7) %>% # starting location (Georgia, USA)
  leaflet::addProviderTiles(providers$OpenStreetMap.BlackAndWhite, group = "OSM (Default)") %>% # default basemap
  leaflet::addProviderTiles(providers$Esri.WorldTopoMap, group = "Terrain") %>% # additional basemap
  leaflet::addProviderTiles(providers$OpenTopoMap, group = "Topographic") %>% # additional basemap
  leaflet::addPolygons(data = CoV_GA_WGS84, color = "black", weight = 1, smoothFactor = 0.5, opacity = 1,
                       fillOpacity = 0.5, fillColor = ~pal_cum(cumulative), popup = ~popup1, # stroke = FALSE
                       highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = TRUE),
                       group = "Cumulative Cases"
  ) %>%
  leaflet::addPolygons(data = CoV_GA_WGS84, color = "black", weight = 1, smoothFactor = 0.5, opacity = 1,
                       fillOpacity = 0.5, fillColor = ~pal_rate(cumrate), popup = ~popup2, # stroke = FALSE
                       highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = TRUE),
                       group = "Cumulative Rate"
  ) %>%
  leaflet::addLayersControl(baseGroups = c("OSM (Default)", "Terrain", "Topographic"),
                            overlayGroups = c("Cumulative Cases", "Cumulative Rate"),
                            options = layersControlOptions(collapsed = FALSE)
  ) %>% # layer selection
  addLegend("topright", pal = pal_cum, values = ~cumulative,
            title = "Cumulative Cases",
            opacity = 1,
            group = "Cumulative Cases") %>%
  addLegend("topright", pal = pal_rate, values = ~cumrate,
            title = "Cumulative Rate per 100,000",
            opacity = 1,
            group = "Cumulative Rate") %>%
  leaflet::hideGroup(c("Cumulative Cases", "Cumulative Rate")) %>% # no data shown (default)
  leaflet::addMiniMap() # add mini map
ga_m1 # display leaflet plot

###############
# Exportation #
###############

# export leaflet plot
htmlwidgets::saveWidget(ga_m1, file = "figures/COVID_Georgia_Leaflet.html", selfcontained = TRUE)

# spatial COVID-19 data in Georgia 
save(CoV_GA_pop, CoV_GA_UTM17N, CoV_GA_WGS84, CoV_GA_df, file = "data/CoV_GA.Rdata")
# -------------------- END OF CODE -------------------- #