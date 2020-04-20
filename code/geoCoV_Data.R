# ---------------------------------------- #
# SARS-CoV-2 open data with geography
# ---------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: April 19, 2020
#
# Recently modified by: 
# Recently modified on: 
#
# Notes:
# 1) Curated data sets for US SARS-CoV-2 Data
# ------------------------------- #

############
# PACKAGES #
############

loadedPackages <- c("readr")
invisible(lapply(loadedPackages, require, character.only = T))

####################
# DATA IMPORTATION #
####################

# the URL for the Johns Hopkins University data page (U.S. Counties, confirmed cases)
jhu_url <- paste("https://raw.githubusercontent.com/CSSEGISandData/", 
                 "COVID-19/master/csse_covid_19_data/", 
                 "csse_covid_19_time_series/",
                 "time_series_covid19_confirmed_US.csv",
                 sep = "")

# the URL for the New York Times data page (U.S. Counties, confirmed cases)
nyt_url <- paste("https://raw.githubusercontent.com/nytimes/", 
                 "covid-19-data/master/us-counties.csv",
                 sep = "")

# the URL for the GitHub Data Packaged Core Datasets (U.S. Counties, confirmed cases)
dpc_url <- paste("https://raw.githubusercontent.com/datasets/", 
                 "covid-19/master/data/us_confirmed.csv",
                 sep = "")

# the URL for the COVID Tracking Project (US COVID-19 testing and mortality by state)
ctp_url <- "https://covidtracking.com/api/v1/states/current.csv"


covid_jhu <- readr::read_csv(jhu_url)
covid_nyt <- readr::read_csv(nyt_url)
covid_dpc <- readr::read_csv(dcp_url)
covid_ctp <- readr::read_csv(ctp_url)


head(covid_jhu)
head(covid_nyt)
head(covid_dpc)
head(covid_ctp)

# -------------------- END OF CODE -------------------- #