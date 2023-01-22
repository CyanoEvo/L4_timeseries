# metadata processing scripts

#### import ####

# dependencies -------------------------------------------------------------------
library(viridis)
library(tidyverse)
library(lubridate)
library(colorspace)
library(ggplot2)
library(Hmisc)
library(MBA)
library(reshape2)
library(egg)
library(httr)
library(purrr)

theme_ro <- function(){
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 12))
}

theme_gss <- function() {
  theme_light() + theme(panel.grid.major.y = element_line(size = 0.05, colour = "grey40"),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.border = element_blank(), 
                        axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
                        axis.ticks.y = element_blank(),
                        text = element_text(colour = "grey40", size = 10),
                        axis.line.x = element_line(size = 0.3, colour = "grey40"),
                        plot.subtitle = element_text(size = 10, colour = "grey40"),
                        plot.title = element_text(size = 14, colour = "grey20", face = "bold"),
                        axis.title.x = element_text(size = 10, colour = "grey40"),
                        axis.title.y = element_blank(),
                        strip.background = element_blank(),
                        strip.text = element_text(colour = "grey40")) 
}

# end

# river flow import (api) -------------------------------------------------
## station information
api_rivers <- 'https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info?station=47001,47003,47011,47022,47004,47009&format=json-object&fields=station-information,gdf-start-date,gdf-mean-flow'
api_response_rivers <- httr::GET(url = api_rivers, timeout(10)) 
all_stations <- content(api_response_rivers, "text") %>% jsonlite::fromJSON()
glimpse(all_stations)

## for loop data api calls
station_list <- c("47001", "47003","47011", "47022", "47004", "47009") ## these are the river stations of interest
all_gdf <- data.frame()
for (i in 1:length(station_list)) {
  gdf <- jsonlite::fromJSON(paste0("https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?format=json-object&data-type=gdf&station=",as.numeric(station_list)[i],""))$`data-stream`[c(F,T)] ## call the gauged river flow (i.e. odd rows)
  date <- jsonlite::fromJSON(paste0("https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?format=json-object&data-type=gdf&station=",as.numeric(station_list)[i],""))$`data-stream`[c(T,F)] ## call the date (i.e. even rows)
  id <- rep(as.numeric(station_list)[i],length(gdf)) ## set the station information
  bound <- as.data.frame(cbind(date, gdf, id), col.names = c("date", "gdf", "id"))
  all_gdf <- rbind(all_gdf, bound)
}
all_gdf$date <- as_date(all_gdf$date)
all_gdf$gdf <- as.numeric(as.character(all_gdf$gdf))
all_gdf <- as_tibble(all_gdf) %>% mutate(id = as.character(id))
glimpse(all_gdf)

## all_stations (station metadata) and all_gdf (gauged daily flow time series) and be linked through the 'id' parameter

## Inspecting end date of each time-series
tail(all_gdf %>% filter(id == "47001")) # 2019-09-30
tail(all_gdf %>% filter(id == "47003")) # 1980
tail(all_gdf %>% filter(id == "47011")) # 2019-09-30
tail(all_gdf %>% filter(id == "47022")) # 2018-09-30
tail(all_gdf %>% filter(id == "47004")) # 2019-09-30
tail(all_gdf %>% filter(id == "47009")) # 2019-09-30

#end


# rainfall import ---------------------------------------------------------
## rainfall at gunnislake
tamar_rain <- as_tibble(read.csv("/home/ro/r_projects/l4_project/server_export/tamar_rain_gunnis.csv")) ## path to csv
tamar_rain$date <- as.Date(tamar_rain$date, format = "%d/%m/%Y")
# end
# ctd data import ---------------------------------------------------------
## basic function test 
test_a <- read.table("/home/ro/r_projects/l4_project/server_export/l4_metadata/ctd/l4_asc/040607_l4.asc", header = T, sep = ";") %>% as_tibble() %>% select("Tv290C", "Sal00", "DepSM")
test_a$file <- gsub("_l4.asc", "", list.files(pattern = ".asc")[1])
head(test_a)
## end basic function test

## loop implementation
ctd_con <- tibble()
for(i in 1:length(list.files(path = "/home/ro/r_projects/l4_project/server_export/l4_metadata/ctd/l4_asc/", pattern = ".asc", full.names = T))){
  ctd_dat <- read.table(list.files(path = "/home/ro/r_projects/l4_project/server_export/l4_metadata/ctd/l4_asc/", pattern = ".asc", full.names = T)[i], header = T, sep = ";", row.names = NULL) %>% as_tibble() %>% select("Tv290C", "Sal00", "DepSM") #extract only temperature, salinity, and depth
  ctd_dat$file <- gsub("_l4.asc", "", list.files(path = "/home/ro/r_projects/l4_project/server_export/l4_metadata/ctd/l4_asc/", pattern = ".asc")[i]) # extracts date from file name and adds column (dont require 'full.names' argument here)
  ctd_con <- bind_rows(ctd_con, ctd_dat)
} ## This isn't working because the working directory is not set appropriately 


## some raw files required manual curation. These files could be identified using tail(ctd_con) when loops failed, to find the final .asc file which successfully ran and pinpoint the issue to the following files. Issues included files being tab delimited rather than semi-colon delimited as is the normal. This could be sorted using sed 's/[[:space:]]\{1,\}/;/g' input > output. Issues also arose where there were more data columns than column headers. This could be resolved by adding dummy headers to specific asc files. However this has all been fixed with the asc files now on the server.

## convert file names to dates
ctd_con$date <- as.Date(ctd_con$file, format = "%y%m%d")

# end

# nutrient data import ----------------------------------------------------
nuts <- as_tibble(read.csv("/home/ro/r_projects/l4_project/server_export/l4_metadata/nutrients/l4_nutrients_2000-2019.csv"))[,1:11] #all values in µM; path to csv location on server
nuts$Date <- dmy(nuts$Date)
nuts <- rename(nuts, date = Date, nutrient_depth = Depth, nutrients_frozen = Frozen, nutrients_filtered = Filtered, nutrient_qc = QC, nitrite = NITRITE, nitrate_nitrite = NITRATE.NIT, ammonia = AMMONIA, silicate = SILICATE, phosphate = PHOSPHATE)
nuts$nitrite <- as.numeric(nuts$nitrite) # note that this manipulation is because R reads these columns as characters due to missing values - a warning that NAs are introduced will be issued but this is expected
nuts$ammonia <- as.numeric(nuts$ammonia)
nuts$phosphate <- as.numeric(nuts$phosphate)
#end
# sample data import ------------------------------------------------------
dna_sam <- as_tibble(read.csv("/home/ro/r_projects/l4_project/server_export/L4_time_series_database_RA_final.csv", header = T)[,c(3,4,17,18,19)]) %>% filter(type == 0.45) %>% mutate(date_collected = as.Date(date_collected, format = "%d/%m/%Y"))
colnames(dna_sam) <- c("date", "filter", "ext_conc", "ext_260280", "ext_260230")

# end




#### analysis ####

# river flow analysis -----------------------------------------------------
## gauged daily flow from relevant rivers across full time series
gdf_plot <- ggplot(all_gdf %>% filter(date >= as.Date("2001-01-01")), aes(x = date, y = gdf, fill = id)) + 
  geom_bar(stat = "identity") + 
  scale_x_date(expand = c(0,0), date_labels = "%Y", date_breaks = "year") + 
  scale_fill_brewer(palette = "Paired", labels = c("Tamar", "Lynher", "Tiddy","Plym","Tory Brook")) + 
  ylim(0,350) + # problem is this axis- values also dont match to original
  theme_gss() +
  ggtitle("", subtitle = expression(GDF~m^3~s^-1)) +
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") # note Tory Brook series ended in 2018

## additive total gauged daily flow from all relevant rivers (i.e. same as above but a single line)
gdf_join <- list(all_gdf %>% filter(id == "47001"), all_gdf %>% filter(id == "47004"), all_gdf %>% filter(id == "47009"),all_gdf %>% filter(id == "47011")) %>% reduce(full_join, by = "date") %>% mutate(sum_gdf = gdf.x+gdf.y+gdf.x.x+gdf.y.y)

gdf_sum_plot <- ggplot(gdf_join %>% filter(date >= as.Date("2001-01-01")), aes(x = date, y = sum_gdf)) + 
  geom_line(colour = "black", size = 0.1) + 
  scale_x_date(expand = c(0,0), date_labels = "%Y", date_breaks = "year") + 
  theme_gss() +
  theme(axis.title.x = element_blank()) # note Tory Brook series ended in 2018 and doesn't cover full length of time-series

## autocorrelation function for summative flow
acf_plot <- acf(gdf_join %>% filter(date >= as.Date("2001-01-01")) %>% select(sum_gdf), lag.max = 365*10) # note that a lag max of 10 years is chosen, but this can be increased up to the full length of the time series (19 years)

gam_poly_acf <- ggplot(as_tibble(cbind(acf_plot$lag, acf_plot$acf)), aes(x = V1, y = V2)) + 
  geom_bar(stat = "identity", colour = "black") + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 21), colour = "firebrick1", fill = NA, size = 0.75) +
  xlab("Lag (Days)") + 
  ggtitle("", subtitle = "ACF (River Flow)") +
  theme_gss() # note that the 'k' selected for GAM essentially indicates the number of oscillations you are expecting, i calculate this are the 1 + (n(years)*2), so for 10 years, k is 21, for 5 years, k is 11. 

# end

# relating salinity and river flow ----------------------------------------
## join salinity (ctd) and river flow data by date
gdf_sal <- full_join(ctd_con %>% filter(DepSM == 4, Sal00 >= 33), gdf_join %>% select(date, sum_gdf), by= "date") # note that outlier salinity values were removed here, and ctd data from only 4 m were kept (good coverage of all casts, prevents extra points from all depths of each ctd cast)

gdf_sal_plot <- ggplot(gdf_sal, aes(x = sum_gdf, y = Sal00)) + 
  geom_point(colour = "steelblue4", fill = "steelblue4", alpha = 0.2) + 
  geom_smooth(method = "lm", colour = "firebrick1", fill = "firebrick1", alpha = 0.1, size = 0.5) + 
  xlab(expression(GDF~m^3~s^-1)) +
  ggtitle("",subtitle = "Salinity (PSU)") + 
  scale_x_continuous(limits = c(0, 170)) +
  annotate(geom = "text", x = 150, y = 35.65, label = expression(paste(R^2:~0.17)), size = 3, colour = "grey40") +
  theme_gss()

summary(lm(gdf_sal$sum_gdf ~ gdf_sal$Sal00)) ## R^2: 0.17, p < 0.0001 - this is not a great fit, and can be interpretted as 17% of variation in salinity at l4 being explained by guaged daily flow. There are other important factors influencing salinity at L4, including tide.

## Technical note is that days of maximal flow (summative flow rate >= 200), there is no CTD data - this is perhaps an important data bias resulting from conditions where it is appropriate to sail being negatively associated with rainfall conditions. Buoy data may help to correct for these biases. 

#end

# seasonal stratification analysis -------------------------------------------------
## aggregating ctd data across all depths (i.e. group by date ctd was collected)
ctd_group <- group_by(ctd_con, date)
temp_strat <- summarise(ctd_group %>% filter(between(DepSM, 4, 40)), mean = mean(Tv290C), sd  = sd(Tv290C)) %>% mutate(julian = yday(date)) %>% mutate(year = year(date)) ## adds a julian day column using lubridate

## thermal stratification index is the standard devidation of seawater temperature throughout the water column
strat_plot <-ggplot(temp_strat %>% filter(mean >= 4, sd <= 3), aes(x = julian, y = sd, fill = mean)) + geom_point(shape = 21, alpha = 1, stroke = 0.5) + theme_gss() + scale_fill_distiller(palette = "RdYlBu", ) + xlab("Julian Day") + ggtitle("Thermal Stratification Index", subtitle = "sd(Temperature ºC)") + labs(fill = "Mean\nTemp ºC") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 7), size = 0.25, colour = "black", fill = NA, linetype = 1)

# end





# monthly metadata analysis -----------------------------------------------
## surface tempearture
temp_group <- group_by(ctd_con %>% mutate(month = month(date)) %>% filter(date >= as.Date("2001-01-01")) %>% filter(DepSM == 4) %>% filter(between(Tv290C, 2, 20)) , month) %>% summarise(mean_temp = mean(Tv290C), sd_temp = sd(Tv290C))

surf_temp_monthly <- ggplot(temp_group, aes(x = month, y = mean_temp)) +
  geom_errorbar(aes(x = month, ymin = mean_temp-sd_temp, ymax = mean_temp+sd_temp), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "firebrick1", size = 2) +
  scale_x_discrete(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Temperature (ºC)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## surface salinity
sal_group <- group_by(ctd_con %>% mutate(month = month(date)) %>% filter(date >= as.Date("2001-01-01")) %>% filter(DepSM == 4) %>% filter(between(Sal00, 29, 38)), month) %>% summarise(mean_sal = mean(Sal00), sd_sal = sd(Sal00))

surf_sal_monthly <- ggplot(sal_group, aes(x = month, y = mean_sal)) +
  geom_errorbar(aes(x = month, ymin = mean_sal-sd_sal, ymax = mean_sal+sd_sal), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "firebrick1", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Salinity (PSU)") + 
  scale_y_continuous(limits = c(34,36)) +
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## rainfall
rain_group <- group_by(tamar_rain %>% mutate(month = month(date)), month) %>% summarise(mean = mean(catch_daily_rain_mm), sd = sd(catch_daily_rain_mm))                           

rainfall_monthly <- ggplot(rain_group, aes(x = month, y = mean)) +
  geom_errorbar(aes(x = month, ymin = mean-sd, ymax = mean+sd), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "lightblue", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Rainfall (mm)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## river flow
river_group <- all_gdf %>% mutate(month = month(date)) %>% filter(date >= as.Date("2001-01-01")) %>% group_by(month) %>% summarise(mean = mean(gdf), sd = sd(gdf)) 

river_monthly <- ggplot(river_group, aes(x = month, y = mean)) +
  geom_errorbar(aes(x = month, ymin = mean-sd, ymax = mean+sd), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "lightblue", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = expression(GDF~m^3~s^-1)) + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## nutrients
nuts_group <- group_by(nuts %>% mutate(month = month(date)), month) %>% summarise(mean_nit = mean(nitrate_nitrite, na.rm = T), sd_nit = sd(nitrate_nitrite, na.rm = T), mean_pho = mean(phosphate, na.rm = T), sd_pho = sd(phosphate, na.rm = T), mean_sil = mean(silicate, na.rm = T), sd_sil = sd(silicate, na.rm = T), mean_amm = mean(ammonia, na.rm = T), sd_amm = sd(ammonia, na.rm = T))                                                      

nit_monthly <- ggplot(nuts_group, aes(x = month, y = mean_nit)) +
  geom_errorbar(aes(x = month, ymin = mean_nit-sd_nit, ymax = mean_nit+sd_nit), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "orange", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Nitrate + Nitrite (µM)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pho_monthly <- ggplot(nuts_group, aes(x = month, y = mean_pho)) +
  geom_errorbar(aes(x = month, ymin = mean_pho-sd_pho, ymax = mean_pho+sd_pho), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "orange", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Phosphate (µM)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


sil_monthly <- ggplot(nuts_group, aes(x = month, y = mean_sil)) +
  geom_errorbar(aes(x = month, ymin = mean_sil-sd_sil, ymax = mean_sil+sd_sil), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "orange", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Silicate (µM)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


amm_monthly <- ggplot(nuts_group, aes(x = month, y = mean_amm)) +
  geom_errorbar(aes(x = month, ymin = mean_amm-sd_amm, ymax = mean_amm+sd_amm), size = 0.2, width = 0.2) + 
  geom_point(shape = 21, alpha = 1, stroke = 0.5, fill = "orange", size = 2) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  ggtitle("", subtitle = "Ammonia (µM)") + 
  theme_ro() +
  theme(axis.title.x = element_blank()) + theme(aspect.ratio=0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## grid plot
monthly_env <- ggarrange(surf_temp_monthly, rainfall_monthly, nit_monthly, pho_monthly, surf_sal_monthly, river_monthly, sil_monthly, amm_monthly, ncol = 4)
monthly_env_sig <- ggarrange(surf_temp_monthly, nit_monthly, river_monthly, ncol = 1)

# end


# relation


# metadata sample relation ------------------------------------------------
## here we aim to link all metadata for a given sampling day to the sample itself, creating an object which can be bound to the sample_data of the phyloseq object. 

## set ctd data object (ctd_con) aligns very poorly witht the dates associated with the sampling data. By definition, there must be ctd data associated with each sample, the question is whether a) the recorded dates of the ctd casts are sometimes erroneous, b) the recorded dates of the filters are erroneous (or if, for example they represent the actual date of filtration rather than collection, and in some cases filtration did not occur on the same date as collection), or c) the ctd dataset pulled from the WCO website is incomplete, and there are further ctd profiles to be obtained. 
sam_tib <- as_tibble(sample_data(l4_clean))

nearest_ctd <- vector(length = length(sam_tib$date))
for (i in 1:length(sam_tib$date)){
  x <- sam_tib$date
  y1 <- c(ctd_con$date, as.Date("2000-01-01")) # adding a dummy date as there are samples which exist before first CTD breaking loop (this is dealt with using NA)
  y <- y1[y1 <= x[i]]
  nearest_ctd[i] <- as.character(y[which.min(abs(y - x[i]))])
  nearest_ctd[nearest_ctd == "2000-01-01"] <- NA
}
cbind(as.character(sam_tib$date), nearest_ctd)
min(sam_tib$date)
min(ctd_con$date)

sam_tib$nearest_ctd <- as.Date(nearest_ctd)
sam_tib %>% select(date, nearest_ctd) # noice.
sam_tib$shift <- as.numeric(sam_tib$nearest_ctd - sam_tib$date) # this outputs the a column of the absolute difference between the ctd date and the filtration date (so can be used to discard data with a bigger difference)
sam_tib %>% select(date, nearest_ctd, shift)

## some values are very high, currently we will use a threshold of 2.5 (i.e. if there is more than 2.5 days difference between the ctd date and the sample date, the ctd data will be disregarded, otherwise it will be retained). This will be done by converting the nearest_ctd column to NA where the shift value is greater than abs(2.5), so that when the sample data and the ctd data are joined using 'left_join' ctd data will not be pulled for samples with a high shift value. IMPORTANT - this is not a terminal solution to the ctd data challenge, resolving the *reason* behind the misalignment of dates is important. 

## there are 23 samples with unreasonable gaps between ctd data and filtration - CTD data are not available of WCO portal for any 'nearer' dates. 
ctd_issue_dates <- sam_tib %>% select(date, nearest_ctd, shift) %>% filter(shift <= -2.5)

## manually updating nearest CTD for samples which may have been mislablled
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2005-04-11"] <- "2005-04-26"
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2009-10-19"] <- "2009-10-26"
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2010-08-11"] <- "2010-08-16"
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2013-01-21"] <- "2013-02-07"
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2013-01-21"] <- "2013-02-07"
sam_tib$nearest_ctd[sam_tib$nearest_ctd == "2013-04-19"] <- "2018-04-19"

sam_tib$nearest_ctd <- replace(sam_tib$nearest_ctd, abs(sam_tib$shift) > 2.5, NA)
sam_tib <- sam_tib %>% select(-Tv290C, -DepSM, -Sal00)

## now join the ctd data to the sample data (sam_tib) based on the nearest_ctd column in sam_tib
### add a dummy column with nearest_ctd name for the join
ctd_con$nearest_ctd <- ctd_con$date
### join sample data and ctd data based on ctd readings at 4m depth (again, you may want to revise this, the reason the 'surface' value isn't chosen is that the first depth of reasonable measurement varies between casts).
sam_tib_join <- left_join(sam_tib, ctd_con %>% filter(DepSM == 4) %>% select(-date), by = "nearest_ctd") # note selection of 4m depth, this could be changed.


## gauged daily flow data can now be added to the sample sample data object. In doing this we want to separate flow data from each river, and ultimately the cumulative flow, thus there is some initial wrangling to be done. 
unique(all_gdf$id)
only_47001 <- all_gdf %>% filter(id == "47001")
colnames(only_47001) <- c("date", "gdf_47001", "id")
only_47003 <- all_gdf %>% filter(id == "47003")
colnames(only_47003) <- c("date", "gdf_47003", "id")
only_47011 <- all_gdf %>% filter(id == "47011")
colnames(only_47011) <- c("date", "gdf_47011", "id")
only_47022 <- all_gdf %>% filter(id == "47022")
colnames(only_47022) <- c("date", "gdf_47022", "id")
only_47004 <- all_gdf %>% filter(id == "47004")
colnames(only_47004) <- c("date", "gdf_47004", "id")
only_47009 <- all_gdf %>% filter(id == "47009")
colnames(only_47009) <- c("date", "gdf_47009", "id")

anchor_date <- all_gdf %>% filter(date >= "2001-02-05") %>% select(date) 
colnames(anchor_date) <- "date"

## perform a nested left join
sep_flow <- left_join(anchor_date, only_47001 %>% select(date, gdf_47001), by="date") %>%
  left_join(., only_47003 %>% select(date, gdf_47003), by="date") %>%
  left_join(., only_47011 %>% select(date, gdf_47011), by = "date") %>%
  left_join(., only_47022 %>% select(date, gdf_47022), by="date") %>%
  left_join(., only_47004 %>% select(date, gdf_47004), by="date") %>%
  left_join(., only_47009 %>% select(date, gdf_47009), by="date") 

sep_flow$total_gdf <- sep_flow$gdf_47001 + sep_flow$gdf_47004 + sep_flow$gdf_47009 + sep_flow$gdf_47011 # adds composite river flow score

sep_flow <- distinct(sep_flow)

## join to sample data object - this metadata object can now be appended  to the phyloseq object
sam_tib_join_gdf <- left_join(sam_tib_join, sep_flow, by = "date") %>% mutate(month = factor(as.character(month), levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")))     

## add nutrient data
## this must first be summarised by date, as replicates currently have individual rows in the nuts dataframe
nuts_sum <- nuts %>% 
  select(date, nitrite, nitrate_nitrite, ammonia, silicate, phosphate) %>%
  group_by(date) %>%
  summarise(mean_nitrite = mean(nitrite),
            mean_nirate_nitrite = mean(nitrate_nitrite),
            mean_ammonia = mean(ammonia),
            mean_silicate = mean(silicate),
            mean_phosphate = mean(phosphate))
  
sam_tib_join_gdf_nuts <- left_join(sam_tib_join_gdf, nuts_sum %>% mutate(nearest_ctd = nuts_sum$date) %>% select(-date) , by = "nearest_ctd") # this needs to be done by nearest ctd again.


#binding to a phyloseq sample data table
sam_tib_join_gdf_nuts_df <- as.data.frame(sam_tib_join_gdf_nuts)
rownames(sam_tib_join_gdf_nuts_df) <- rownames(sample_data(l4_clean))
l4_clean_sd <- l4_clean
sample_data(l4_clean_sd) <- sam_tib_join_gdf_nuts_df
# note object is renames l4_clean_sd to prevent l4_clean overwrite