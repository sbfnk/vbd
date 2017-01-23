library('tidyverse')
library('magrittr')
library('readxl')

## read data
incidence <- read_excel("table\ for\ modeling\ ZIKV.xlsx", "Notified cases bahia salv", skip=2)

###################
## data cleaning ##
###################

incidence %<>% setNames(make.names(tolower(names(incidence)), unique=TRUE))
names(incidence)[1:2] <- c("year", "week") # date headers

incidence %<>%
  select(-starts_with("NA")) %>% # remove columns without header
  filter(!is.na(week)) %>% # remove empty rows
  mutate(year=as.integer(gsub("[^0-9]", "", year))) %>% # convert year column to integer
  fill(year) %>% # fill NAs
  mutate(date=as.Date(paste(year, week, 1, sep="-"), "%Y-%W-%u")) # make dates from Mondays of each week

## make long data set
incidence %<>%
  gather(region, incidence, grep("(bahia|salvador)", colnames(incidence)))

## bring back infections / outcomes
incidence %<>%
  mutate(infection=sub("^[a-z\\.]*", "", region)) %>%
  mutate(region=sub("\\..*$", "", region)) %>%
  mutate(infection=
           case_when(.$infection == "" ~ "zikv",
                     .$infection == "1" ~ "chikv",
                     .$infection == "2" ~ "denv",
                     .$infection == "3" ~ "exantemic",
                     .$infection == "4" ~ "mn3_exantemic",
                     .$infection == "5" ~ "mn5_exantemic"))

##################################
## data preparation for fitting ##
##################################

case_data <- incidence %>%
  filter(region == "salvador" & infection == "zikv") %>%
  filter(date < "2015-12-10") %>%
  mutate(time=1:n()) %>%
  rename(value=incidence) %>%
  select(time, value)

serology_data <- data.frame(time=max(case_data$time), 
                            value=427)

saveRDS(list(Incidence=case_data, Serology=serology_data), "fit_data.rds")
