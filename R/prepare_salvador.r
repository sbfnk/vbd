library('tidyverse')
library('magrittr')
library('readxl')
library('lubridate')

## read data
incidence <- read_excel("table\ for\ modeling\ ZIKV.xlsx", "Notified cases bahia salv", skip=2)

###################
## data cleaning ##
###################

incidence %<>% setNames(make.names(tolower(names(incidence)), unique=TRUE))
names(incidence)[1:2] <- c("year", "week") # date headers

incidence %<>%
  select(-starts_with("NA")) %>% # remove columns without header
  dplyr::filter(!is.na(week)) %>% # remove empty rows
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
  dplyr::filter(region == "salvador" & infection == "zikv") %>%
  dplyr::filter(date < "2015-11-01")

max_case_date <- max(case_data$date)

case_data %<>%
  mutate(time=1:n()) %>%
  rename(value=incidence) %>%
  select(time, value)

serology_time <- max(case_data$time) +
  round(as.integer(as.Date("2016-04-01") - max_case_date) / 7)

serology_data <- data.frame(time=serology_time, value=401)

saveRDS(list(Incidence=case_data, Serology=serology_data), "fit_data.rds")

##############
## serology ##
##############

serology <- read_excel("table\ for\ modeling\ ZIKV.xlsx", "ZIKV 2016 table", col_names=FALSE, skip=1) %>%
  rename(sample=X__2, date=X__3, zika=X__4) %>%
  mutate(date=as.Date(date),
         month=date-mday(date)+1)

sero_summary <- serology %>%
  group_by(month) %>%
  summarise(pos=sum(zika == "positivo"),
            neg=sum(zika == "negativo"),
            n=n()) %>%
  ungroup %>%
  filter(!is.na(month))

errors <- binom.confint(sero_summary$pos, sero_summary$n, method="wilson")

sero_summary <- cbind(sero_summary, errors %>% select(-x, -n, -method))

ggplot(sero_summary,
       aes(x=month, y=mean, ymin=lower, ymax=upper)) +
  geom_point() +
  geom_errorbar()

ggplot(serology, aes(x=month, fill=sample))+geom_bar(position="stack")+scale_fill_brewer(palette="Set1")
