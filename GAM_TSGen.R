## Statistical analysis for manuscript "Shifting Phenology as a Key Driver of Shelf Zooplankton Population Variability"

# R Code for generating complete timeseries of Calanus finmarchicus
# Isabel Honda
# May 2024

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("mgcv")
})


# Importing Data
data <- # Read in MARMAP/EcoMon Data

selection <- data %>% select(,c("strata_46","month", "doy","year","lon","lat","calfin_100m3")) # Selecting variables (calfin_100m3 == C. finmarchicus abundance per 100m^3)

all_log <- selection
all_log$calfin_100m3 <- log((selection$calfin_100m3)/100 + 1) # Log transforming the data 

resps <- c("calfin_100m3")
knots <- list(doy = c(0.5, 366.5))
strata = c(37, 42, 48) # Strata 37 == Wilkinson Basin | Strata 42 == Jordan Basin | Strata 48 == Georges Basin

for (resp in resps) {
  resp_est <- data.frame(year = rep(min(data$year):max(data$year),each=365), doy = 1:365)
  
  for (s_input in strata) {
    
    if (s_input == 48) {
      Xvar <- all_log %>% filter(strata_46 == 38 | strata_46 == 39) 
    }  else {
      Xvar <- all_log %>% filter(strata_46 == s_input) 
    }
    Xvar <- na.omit(Xvar)
    
    preds_unres <- list(
      "s(doy, k = 12, bs = 'cc')",
      "s(year, k = 10)",
      "ti(doy, year, bs = c('cc', 'tp'), k = c(12, 10))"
    )
    
    gam_terms = data.frame(year = rep(min(data$year):max(data$year),each=365), julian = 1:365)
    
    curGam <- bam(formula(paste(resp, "~", paste(preds_unres,collapse="+"))), data=Xvar, method = 'fREML', 
                  discrete = TRUE, select = TRUE, na.action="na.omit", knots = knots)
    terms <- predict(curGam,gam_terms,type='response',se.fit=TRUE)
    resp_est <- cbind(resp_est,terms)
    colnames(resp_est)[length(names(resp_est))-1] <- paste0('s',s_input)
    colnames(resp_est)[length(names(resp_est))] <- paste0('se_',s_input)
  }
   write.csv(resp_est,row.names=FALSE) # Save file as .csv
}
