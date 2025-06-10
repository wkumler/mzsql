
library(tidyverse)
library(RaMS)

# Write imzML data to Arrow for R ----
# 10.1007/s00216-012-5841-x
# Maybe also https://www.ebi.ac.uk/metabolights/editor/MTBLS313/descriptors?
library(Cardinal)
msi_files <- list.files("imzml_data", pattern = ".imzML", full.names = TRUE)[1:3]
all_msi <- pbapply::pblapply(msi_files, function(file_i){
  msi <- readImzML(file_i)
  
  coords <- coord(msi)
  mz_list <- as.list(mz(msi))
  int_list <- as.list(intensity(msi))
  
  result <- data.frame(
    filename = basename(file_i),
    scan = rep(seq_along(mz_list), lengths(mz_list)),
    x = rep(coords$x, lengths(mz_list)),
    y = rep(coords$y, lengths(mz_list)),
    mz = unlist(mz_list),
    int = unlist(int_list)
  )
  return(result)
})
v <- do.call("rbind", all_msi)
write_dataset(v, "imzml_data/pqds")

lapply(all_msi, function(df_i){
  df_i %>%
    arrange(desc(int)) %>%
    head()
})

# Write ms3 data to DuckDB for Python ----
# https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001661
# 10.1021/acs.analchem.0c04895
library(DBI)
ms3_files <- list.files("ms3_data", pattern = "mzML", full.names = TRUE)
new_grabMzmlMS3 <- function (xml_data, rtrange, file_metadata, incl_polarity) {
  ms3_xpath <- "//d1:spectrum[d1:cvParam[@name=\"ms level\" and @value=\"3\"]]"
  ms3_nodes <- xml2::xml_find_all(xml_data, ms3_xpath)
  rt_vals <- RaMS:::grabSpectraRt(ms3_nodes)
  fsnodes <- xml2::xml_find_all(ms3_nodes, "d1:scanList/d1:scan/d1:cvParam[@name=\"filter string\"]")
  shortstrings <- gsub("FTMS \\+ c NSI d Full ms3 |@hcd30.00.*", "", xml_attr(fsnodes, "value"))
  prepremz_vals <- as.numeric(gsub("@.*", "", shortstrings))
  premz_vals <- as.numeric(gsub(".*cid30.00 ", "", shortstrings))
  mz_vals <- RaMS:::grabSpectraMz(ms3_nodes, file_metadata)
  int_vals <- RaMS:::grabSpectraInt(ms3_nodes, file_metadata)
  data.table(rt = rep(rt_vals, lengths(mz_vals)), 
             prepremz = rep(prepremz_vals, lengths(mz_vals)), 
             premz = rep(premz_vals, lengths(mz_vals)), 
             fragmz = unlist(mz_vals), 
             int = as.numeric(unlist(int_vals)), 
             voltage = rep(rep(30, length(mz_vals)), lengths(mz_vals)))
}

new_grabMzmlData <- function (filename, grab_what, verbosity = 0, incl_polarity = FALSE, 
                              mz = NULL, ppm = NULL, rtrange = NULL, prefilter = -1) {
  if (verbosity > 1) {
    cat(paste0("\nReading file ", basename(filename), "... "))
    last_time <- Sys.time()
  }
  xml_data <- xml2::read_xml(filename)
  RaMS:::checkNamespace(xml_data)
  RaMS:::checkFileType(xml_data, "mzML")
  rtrange <- RaMS:::checkRTrange(rtrange)
  prefilter <- RaMS:::checkProvidedPrefilter(prefilter)
  output_data <- list()
  if ("everything" %in% grab_what) {
    extra_grabs <- setdiff(grab_what, "everything")
    if (any(c("MS1", "MS2", "BPC", "TIC", "metadata") %in% 
            extra_grabs) && verbosity > 0) {
      message(paste("Heads-up: grab_what = `everything` includes", 
                    "MS1, MS2, BPC, TIC, and meta data"))
      message("Ignoring duplicate specification")
    }
    grab_what <- unique(c("MS1", "MS2", "BPC", "TIC", "metadata", 
                          extra_grabs))
  }
  if (any(c("MS1", "MS2", "MS3", "DAD", "EIC", "EIC_MS2", "EIC_MS3", 
            "chroms") %in% grab_what)) {
    file_metadata <- RaMS:::grabMzmlEncodingData(xml_data)
  }
  if ("MS1" %in% grab_what) {
    if (verbosity > 1) 
      last_time <- RaMS:::timeReport(last_time, text = "Reading MS1 data...")
    output_data$MS1 <- RaMS:::grabMzmlMS1(xml_data = xml_data, rtrange = rtrange, 
                                          file_metadata = file_metadata, incl_polarity = incl_polarity, 
                                          prefilter = prefilter)
  }
  if ("MS2" %in% grab_what) {
    if (verbosity > 1) 
      last_time <- RaMS:::timeReport(last_time, text = "Reading MS2 data...")
    output_data$MS2 <- RaMS:::grabMzmlMS2(xml_data = xml_data, rtrange = rtrange, 
                                          incl_polarity = incl_polarity, file_metadata = file_metadata)
  }
  if ("MS3" %in% grab_what) {
    if (verbosity > 1) 
      last_time <- RaMS:::timeReport(last_time, text = "Reading MS3 data...")
    output_data$MS3 <- new_grabMzmlMS3(xml_data = xml_data, rtrange = rtrange, 
                                       incl_polarity = incl_polarity, file_metadata = file_metadata)
  }
  if ("BPC" %in% grab_what) {
    if (verbosity > 1) 
      last_time <- RaMS:::timeReport(last_time, text = "Reading BPC...")
    output_data$BPC <- RaMS:::grabMzmlBPC(xml_data = xml_data, rtrange = rtrange, 
                                          incl_polarity = incl_polarity)
  }
  if ("TIC" %in% grab_what) {
    if (verbosity > 1) 
      last_time <- RaMS:::timeReport(last_time, text = "Reading TIC...")
    output_data$TIC <- RaMS:::grabMzmlBPC(xml_data = xml_data, rtrange = rtrange, 
                                          incl_polarity = incl_polarity, TIC = TRUE)
  }
  if ("metadata" %in% grab_what) {
    if (verbosity > 1) {
      last_time <- RaMS:::timeReport(last_time, text = "Reading file metadata...")
    }
    output_data$metadata <- RaMS:::grabMzmlMetadata(xml_data = xml_data)
  }
  if (verbosity > 1) {
    time_total <- round(difftime(Sys.time(), last_time), 
                        digits = 2)
    cat(time_total, units(time_total), "\n")
  }
  output_data
}

ms3data <- new_grabMzmlData(ms3_files[1], grab_what = c("BPC", "MS1", "MS2", "MS3", "metadata"), verbosity = 2)

con <- dbConnect(duckdb::duckdb(), "ms3_data/ms3_data.duckdb")
dbWriteTable(con, "MS1", ms3data$MS1[rt%between%c(20, 25)], overwrite = TRUE)
dbWriteTable(con, "MS2", ms3data$MS2, overwrite = TRUE)
dbWriteTable(con, "MS3", ms3data$MS3, overwrite = TRUE)
dbWriteTable(con, "BPC", ms3data$BPC, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM MS1 ORDER BY int DESC LIMIT 10")
dbDisconnect(con)

con <- dbConnect(duckdb::duckdb(), "ms3_data/ms3_data.duckdb")
ms1_data <- dbGetQuery(con, "SELECT * FROM MS1 WHERE mz BETWEEN 282.1186 AND 282.1208")
ms2_data <- dbGetQuery(con, "SELECT * FROM MS2 WHERE premz BETWEEN 282.1186 AND 282.1208")
ms3_data <- dbGetQuery(con, "SELECT * FROM MS3 WHERE prepremz BETWEEN 282.1186 AND 282.1208 AND
                             premz BETWEEN 166.0718 AND 166.0731")
dbDisconnect(con)

# ggplot(ms1_data) + geom_line(aes(x=rt, y=int))
# ggplot(ms2_data) + geom_point(aes(x=fragmz, y=int))
# ggplot(ms3_data) + geom_point(aes(x=fragmz, y=int))



# Write chrom data to SQLite for Julia ----
ms_files <- "~/../Desktop/Will/PARAGON/mzMLs/pos" %>%
  list.files(pattern="Smp_Amm.*noIS", full.names = TRUE) %>%
  str_subset("Blk", negate = TRUE)
out_files <- paste0("chrom_data/", basename(ms_files))
minifyMSdata(ms_files, out_files, mz_include = c(147.0764, 147.0764+0.997035, 147.0764+0.997035*2,
                                                 118.0865, 118.0865+0.997035,
                                                 152.0567, 152.0567+0.997035*5), ppm = 20)

ms_files <- list.files("chrom_data", pattern = ".mzML", full.names = TRUE)
msdata <- grabMSdata(ms_files, grab_what = "MS1")
msdata$MS1$timepoint <- str_extract(msdata$MS1$filename, "T\\d+")

library(DBI)
con <- dbConnect(RSQLite::SQLite(), "chrom_data/chrom_data.sqlite")
dbExecute(con, "DROP TABLE IF EXISTS MS1")
dbWriteTable(con, "MS1", msdata$MS1, overwrite = TRUE)
dbDisconnect(con)

sapply(c(147.0764, 147.0764+0.997035, 147.0764+0.997035*2), function(mz_i){
  msdata$MS1[mz%between%pmppm(mz_i, 10)][rt%between%c(10, 12)] %>%
    qplotMS1data(color_col = "timepoint") %>%
    print()
})
