#######   Seahorse-analyzeR   ########

##For more information, please visit Github: 
##For additional questions or remarks: 
##please contact Daan Viering at daan.viering@radboudumc.nl

#Set the right working directory (this command only works in Rstudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Source the necessary functions
source('configure.R')

#---------------------------read data---------------------------------
#If you have put all files containing the rate data in a folder 'rate_data',
#you can bind these together here
source('bind_all_rate_data.R') 
#This will also give the columns the proper column names
#The rate data can be obtained by going to the Data view in Wave 2.3.0,
#navigating to the Rate tab, selecting all data, and pasting it in a .txt file
#Make sure the decimal seperator is a '.' before binding the files together.

#Read the obtained file containing all data
dt_raw <- fread("all_rate_data.txt")


#---------------------------tidy data---------------------------------
#Removes wells that contain background measurements
#Also contains the option to transform or remove 'negative OCRs', 
#since these represent a theoretical impossibility in mammalian cells 
#and might interfere with downstream analyses
background <- c("A", "H")
dt1 <- remove_wells(dt_raw, method = "without transformation",
                       background = background)

#Define which cell lines are the control cell lines
control_ids <- c("c1", "c2", "c3")

#Add columns to dt (look in the file to see input specifications)
source('add_columns.R')


#--------------calculate and add group level data--------------------
#Calculate median of each interval for each sample
int_medians <- calculate_int_median(dt1)
dt2 <- left_join(dt1, int_medians, by = "id_int_F") #add it to dataframe

#Calculate average for each Interval (for each FCCP) for each cell culture 
#and add data to dataframe
int_F_averages <- average_cc_int_F(dt2)
dt3 <- left_join(dt2, int_F_averages, by = "cc_int_F")


#--------------keep only FCCP concentrations with max OCR----------------
#Calculate maximum OCR to find the 'optimal' FCCP concentration
ocr_maxes <- calculate_max_ocr(dt3)
dt4 <- left_join(dt3, ocr_maxes, by = "cc_int")
dt4 <- remove_other_fccp(dt4)


#-----------------quality control max_resp & spare_capacity------------------
#Find wells in which the normal mito-stress test response was observed
#i.e. OCR(Int3) > OCR(Int1) > OCR(Int2) > OCR(Int4)
#and remove the wells from the large datatable (dt5) if they failed.
dt_int_measures <- pivot_ints_to_cols(dt4)
dt_passed_int_measures <- find_stress_passes(dt_int_measures)
print(paste0(length(dt_passed_int_measures$sample_id),
             " wells removed because Int3 was NA, ",
             length(dt_int_measures$sample_id) - 2*length(dt_passed_int_measures$sample_id),
             " wells removed because they failed to follow stress test" ))

#dt7 <- add_cols_dt6(dt_passed_int_measures)
dt_without_fails <- remove_stress_fails(dt4, dt_passed_int_measures)
#dt10 <- left_join(dt9, dt7, by = "sample_id")


dt_with_max_resp1 <- add_bioenergetic_cols1(dt_passed_int_measures)
dt_with_extra_passes <- find_stress_passes_adj(dt_int_measures)
dt_with_other_bioenerg1 <- add_bioenergetic_cols2(dt_with_extra_passes)
#dt9 <- remove_stress_fails(dt4, dt_with_extra_passes)
dt10 <- left_join(dt_without_fails, dt_with_max_resp1, by = "sample_id")
#dt10 <- left_join(dt10, dt_with_other_bioenerg1, by = "sample_id")

#----------------summarize max_resp & spare_capacity------------------
dt11 <- summarize1(dt10)
dt12 <- summarize2(dt11)
dt13 <- summarize3(dt12)


#----------------------write results to csv file----------------------
write.table(dt13, "All-data_v2.csv", append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

