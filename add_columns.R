#Adds to dt the following necessary columns: 
#      cell_culture, well_n, interval, treatment, fccp, fccp_n2, sample_id
#Also edits the format of fibro_id
#Required columns in dt (as input) are plate_id, fibro_id and well1.
#This code assumes the Seahorse protocol was executed as published 
#in the accompanying paper (i.e., in-plate FCCP titration)
#Therefore, fibro_id should be written as "cell_line + " " + fccp + "/" fccp_n2"

#Determines interval for each measurement and adds this in a new column
create_column_interval <- function(dt) {
  x <- dt1$interval
  y <- dt1$time
  for (i in 1:length(dt1[[1]])){
    if (y[i]<5){
      x[i] <- "Int1"
    } else if (y[i]>4 & y[i]<8) {
      x[i] <- "Int2"
    } else if (y[i]>7 & y[i]<14) {
      x[i] <- "Int3"
    } else if (y[i]>13) {
      x[i] <- "Int4"
    } 
  }
  x
}

#Determines treatment for each measurement and adds this in a new column
create_column_treatment <- function(dt) {
  x <- dt1$treatment
  y <- dt1$interval
  for (i in 1:length(dt1[[1]])){
    if (y[i] == "Int1"){
      x[i] <- "basal"
    } else if (y[i]== "Int2") {
      x[i] <- "oligomycin"
    } else if (y[i]== "Int3") {
      x[i] <- "FCCP"
    } else if (y[i]== "Int4") {
      x[i] <- "Rot+AA"
    } 
  }
  x
}

#Add/edit the columns
dt1$ocr <- as.numeric(gsub(",","",dt1$ocr))
dt1$ecar <- as.numeric(gsub(",","",dt1$ecar))
dt1$ppr <- as.numeric(gsub(",","",dt1$ppr))
dt1$V1 <- NULL
dt1$fccp_n <- substr_right(dt1$fibro_id, 3, 3)
dt1$fccp_n2 <- substr_right(dt1$fibro_id, 1, 1)
dt1$fibro_id <- edit_id(dt1, control_ids)
dt1$cell_culture <- paste(dt1$plate_id, dt1$fibro_id, sep = "_")
dt1$well_n <- as.numeric(substring(dt1$well1, 2))
dt1$well <- paste0(substring(dt1$well1, 1, 1), dt1$well_n)
dt1$well1 <- NULL
dt1$interval <- create_column_interval(dt1)
dt1$treatment <- create_column_treatment(dt1)
dt1$sample_id <- paste(dt1$plate_id, dt1$well, sep = "_")
dt1$fccp_conc[between(dt1$time, 1,  7 )] <- paste0("FCCP0")
dt1$fccp_conc[between(dt1$time, 8,  10)] <- paste0("FCCP", dt1$fccp_n)[ between(dt1$time, 8, 10)]
dt1$fccp_conc[between(dt1$time, 11, 13)] <- paste0("FCCP", dt1$fccp_n2)[between(dt1$time, 11, 13)]
dt1$fccp_conc[between(dt1$time, 14, 16)] <- paste0("FCCP0")
dt1$fccp_n <- NULL
dt1$fccp_n2 <- NULL
dt1$fccp_conc <- as.factor(dt1$fccp_conc)
dt1$cc_int <- paste(dt1$cell_culture, dt1$interval, sep = "_")
dt1$cc_int_F <- paste(dt1$cell_culture, dt1$interval, 
                         dt1$fccp_conc, sep = "_")
dt1$id_int_F <- paste(dt1$sample_id, dt1$interval, 
                         dt1$fccp_conc, sep = "_")

