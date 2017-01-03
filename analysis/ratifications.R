# ----------------
# Load libraries
# ----------------
library(tidyverse)
library(lubridate)
library(rvest)
library(countrycode)
library(magrittr)
library(stringr)

source(file.path(PROJHOME, "shared_functions.R"))


# -----------
# Load data
# -----------
# Get ratification data directly from the UN
treaty.url <- "https://web.archive.org/web/20160403063433/https://treaties.un.org/Pages/ViewDetails.aspx?src=TREATY&mtdsg_no=XVIII-12-a&chapter=18&lang=en"

ratifications <- read_html(treaty.url) %>%
  html_nodes(xpath='//*[@id="ctl00_ContentPlaceHolder1_tblgrid"]') %>%
  html_table(header=TRUE) %>% bind_rows() %>%
  set_colnames(c("participant", "signature", "ratification")) %>%
  mutate(participant = countrycode(participant, "country.name", "country.name")) %>%
  mutate_each(funs(str_extract(str_replace(., "\\t", " "), 
                               "\\d{1,2}\\s+\\w{3}\\s+\\d{4}")),
              c(signature, ratification)) %>%
  mutate_each(funs(date = dmy(.)), c(signature, ratification)) %>%
  mutate_each(funs(year = year(.)), c(signature_date, ratification_date)) %>%
  filter(!is.na(participant))

# Summarize count of ratifications by year
df.rat <- ratifications %>%
  filter(!is.na(ratification_date_year)) %>%
  group_by(ratification_date_year) %>%
  summarise(num.ratified = n()) %>%
  mutate(num.ratified.cum = cumsum(num.ratified),
         ratification_date_year = ymd(paste0(ratification_date_year, "-01-01")))


# -----------
# Plot data
# -----------
fig.ratified <- ggplot(df.rat, aes(x=ratification_date_year, y=num.ratified.cum)) + 
  geom_line(size=0.75) + 
  labs(x=NULL, y="Number of ratifying countries") +
  coord_cartesian(ylim=c(0, 175)) +
  theme_clean(10)


# ------------------------------------
# Export underlying plot data to CSV
# ------------------------------------
to.csv <- df.rat %>%
  mutate(year = as.integer(year(ratification_date_year))) %>%
  select(year, num_ratified = num.ratified, num_ratifited_cum = num.ratified.cum)

write_csv(to.csv, path=file.path(PROJHOME, "data", "processed",
                                 "ratifications_count.csv"))
