require(dplyr)
require(data.table)

# read file path
all_paths <-
  list.files(pattern = "*.csv",
             full.names = TRUE)

# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")

# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
# change column name
names(all_result)[9] <- "bacterial_strain"
all_result[] <- lapply(all_result, gsub, 
                       pattern='.csv', replacement='')


write.csv(all_result, "./all_defense_finder_systems.csv")
