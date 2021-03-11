#-------------------bind all rate data-----------------------------
#requirements: tidyverse

#remember current working directory
curr_wd <- getwd()

#set the working directory to the folder containing the rate data.
setwd("rate_data")

#create a list of the files here
file_list <- list.files(path = ".")

#initiate a blank data frame
all_rate_data <- data.frame()

#read each file, add column headers, add column with of origin-file
for (i in 1:length(file_list)){
     temp_data <- data.frame()
     temp_data <- fread(file_list[i])
     colnames(temp_data) <- c("time", "well1", "fibro_id", 
                              "time_exact", "ocr", "ecar", "ppr")
     temp_data$plate_id <- sapply(substr_min_right(file_list[i], 4), 
                                  function(x){x[1]})
     all_rate_data <- rbind(all_rate_data, temp_data)
     rm(temp_data)
}

#set back to previous working directory
setwd(curr_wd)

#Write to a file
write.csv(all_rate_data, file = "all_rate_data.txt")

rm(i)
rm(file_list)
rm(curr_wd)
