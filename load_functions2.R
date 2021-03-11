## R-script with useful functions for Testscript

#Removes the last n-characters from string x
substr_min_right <- function(x, n) {
  substr(x, 1, nchar(x)-n)
}

#Extracts last n-characters from string x, until and with n2
#Count starts at 1. (i.e., n2 <- 1 would lead to inclusion of the last character)
substr_right <- function(x, n, n2) {
     substr(x, nchar(x)-n+1, nchar(x)-n2+1)
}

#Editing the fibro_id
edit_id <- function(dt, controls) {
     dt$fibro_id <- substr_min_right(dt$fibro_id, 4)
}

#Removes measurements where the OCR was 0 (i.e. background) or lower than background.
remove_wells <- function(dt, method = "without transformation", background = c("A", "H")) {
     a <- nrow(dt)
     print(paste("Input dt has", a, "rows (measurements)"))
     dt2 <- dt[!(substr_min_right(well1, 2) %in% background)]
     b <- a - nrow(dt2)
     print(paste0(b, " background measurements removed (", b/16, " wells)")) #16 is number of measurement intervals
     if (method == "without transformation") {
          dt3 <- subset.data.frame(dt2, dt2$ocr >= 0)
          c <- nrow(dt2) - nrow(dt3)
          d <- c/nrow(dt2) *100 
          print(paste0(c, " measurements with negative OCR removed (", d, "%)"))
     } else if (method == "with transformation") {
          abs_min_ocr <- abs(min(dt2$ocr))
          dt2$ocr <- dt2$ocr + abs_min_ocr + 0.0001
          print(paste0("All OCR values transformed by adding abs(min(dt$ocr)) and 0.0001"))
     } else {
          print("Error: not a valid input for method-variable in remove_wells")
     }
     dt3
}


#Calculates the median OCR per sample_id
calculate_int_median <- function(dt, group_by_col = "id_int_F") { #Enter dataframe to be analyzed
  dt %>%
    group_by(.data[[!!group_by_col]]) %>%
    summarise(median_int = median(ocr))
}

#Calculates the average of each interval (distinguishing different FCCP concentrations) for each cell culture
average_cc_int_F <- function(dt) {
  dt %>%
    group_by(cc_int_F) %>%
    summarise(cc_int_F_average = mean(median_int))
}

#Calculates the maximum OCR for Int3
calculate_max_ocr <- function(dt) {
  dt %>%
    filter(fccp_conc != "FCCP0") %>%
    group_by(cc_int) %>%
    summarise(Max_ocr = max(cc_int_F_average))
}

#Remove the FCCP concentrations were the OCR was not maxed out
remove_other_fccp <- function(dt) {
  dt %>%
    filter(fccp_conc == "FCCP0" | cc_int_F_average == Max_ocr)
}

#Pivot interval measurements into columns in a summarised version of the table
pivot_ints_to_cols <- function(dt) {
  dt %>%
    group_by(id_int_F) %>%
    summarise(median_int = median(median_int), interval = first(interval), sample_id = first(sample_id)) %>%
    pivot_wider(id_cols = sample_id, names_from = interval, values_from = median_int)
}

#Find those wells that showed the normal mito-stress test response
find_stress_passes <- function(dt) {
  dt %>%
    filter(Int3 > Int1,
           Int1 > Int2,
           Int2 > Int4)
}

#Find those wells that showed the normal mito-stress test response
#Allow wells that had NA at Int3
find_stress_passes_adj <- function(dt) {
  dt %>%
    filter((Int3 > Int1) | is.na(Int3), #OCR values < 0 resulted in NA in remove_wells()
           Int1 > Int2,
           Int2 > Int4)
}

#Removing non-stress-test adhering wells from the original data-table
remove_stress_fails <- function(dt_original, dt_stress_passes) {
  subset(dt_original, dt_original$sample_id %in% dt_stress_passes$sample_id)
}

#Calculate bioenergetic measures and add these in new columns
add_cols_dt6 <- function(dt){
  dt %>%
    mutate(basal_resp = Int1 - Int4) %>%
    mutate(ATP_production = Int1 - Int2) %>%
    mutate(proton_leak = Int2 - Int4) %>%
    mutate(max_resp = Int3 - Int4) %>%
    mutate(spare_capacity = Int3 - Int1) %>%
    mutate(non_mito_resp = Int4)
}

#Calculate bioenergetic measures and add these in new columns
add_bioenergetic_cols1 <- function(dt_for_max_resp){
     dt_for_max_resp %>%
          mutate(max_resp = Int3 - Int4) %>%
          mutate(spare_capacity = Int3 - Int1) %>%
          select(-Int1, -Int2, -Int3, -Int4)
}

#Calculate bioenergetic measures and add these in new columns
add_bioenergetic_cols2 <- function(dt_for_rest){
     dt_for_rest %>%
          mutate(basal_resp = Int1 - Int4) %>%
          mutate(ATP_production = Int1 - Int2) %>%
          mutate(proton_leak = Int2 - Int4) %>%
          mutate(non_mito_resp = Int4)
}



#-------------------------summarizing functions--------------------------
summarize1 <- function(dt) {
  dt %>%
    group_by(id_int_F) %>%
    summarise(median_int = median(median_int), 
              interval = first(interval), 
              well = first(well),
              well_n = first(well_n),
              sample_id = first(sample_id), 
              plate_id = first(plate_id), 
              fibro_id = first(fibro_id), 
              cell_culture = first(cell_culture),
              max_resp = first(max_resp), 
              spare_capacity = first(spare_capacity), 
              log_max_resp = log(first(max_resp)), 
              reciproc_max_resp = 1/first(max_resp))
}


summarize2 <- function(dt) {
  dt %>%
    group_by(sample_id) %>%
    summarise(median_int = median(median_int), 
              well = first(well), 
              well_n = first(well_n),
              plate_id = first(plate_id), 
              fibro_id = first(fibro_id), 
              cell_culture = first(cell_culture),
              max_resp = first(max_resp),
              spare_capacity = first(spare_capacity),
              log_max_resp = log(first(max_resp)), 
              reciproc_max_resp = 1/first(max_resp))
}


summarize3 <- function(dt) {
  dt %>%
    group_by(cell_culture) %>%
    summarise(mean_max_resp = mean(max_resp),
              plate_id = first(plate_id),
              fibro_id = first(fibro_id))
}

summarize4 <- function(dt) {
     dt %>%
          group_by(id_int_F) %>%
          summarise(median_int = median(median_int), 
                    interval = first(interval), 
                    well = first(well),
                    well_n = first(well_n),
                    sample_id = first(sample_id), 
                    plate_id = first(plate_id), 
                    fibro_id = first(fibro_id), 
                    cell_culture = first(cell_culture),
                    max_resp = first(max_resp), 
                    basal_resp = first(basal_resp),
                    ATP_production = first(ATP_production), 
                    proton_leak = first(proton_leak), 
                    spare_capacity = first(spare_capacity), 
                    non_mito_resp = first(non_mito_resp),
                    Int1 = first(Int1), 
                    Int2 = first(Int2), 
                    Int3 = first(Int3), 
                    Int4 = first(Int4),
                    log_max_resp = log(first(max_resp)), 
                    reciproc_max_resp = 1/first(max_resp))
}

summarize5 <- function(dt) {
     dt %>%
          group_by(sample_id) %>%
          summarise(median_int = median(median_int), 
                    well = first(well), 
                    well_n = first(well_n),
                    plate_id = first(plate_id), 
                    fibro_id = first(fibro_id), 
                    cell_culture = first(cell_culture),
                    max_resp = first(max_resp), 
                    basal_resp = first(basal_resp),
                    ATP_production = first(ATP_production), 
                    proton_leak = first(proton_leak), 
                    spare_capacity = first(spare_capacity), 
                    non_mito_resp = first(non_mito_resp),
                    Int1 = first(Int1), 
                    Int2 = first(Int2), 
                    Int3 = first(Int3), 
                    Int4 = first(Int4),
                    log_max_resp = log(first(max_resp)), 
                    reciproc_max_resp = 1/first(max_resp))
}

summarize6 <- function(dt) {
     dt %>%
          group_by(cell_culture) %>%
          summarise(mean_max_resp = mean(max_resp),
                    mean_basal_resp = mean(basal_resp),
                    mean_ATP_production = mean(ATP_production),
                    mean_proton_leak = mean(proton_leak),
                    mean_spare_capacity = mean(spare_capacity),
                    mean_non_mito_resp = mean(non_mito_resp),
                    plate_id = first(plate_id),
                    fibro_id = first(fibro_id))
}


#-----------------------unused functions---------------------------
filter1 <- function(dt) {
     dt %>%
          filter(fibro_id %in% control_ids)
}