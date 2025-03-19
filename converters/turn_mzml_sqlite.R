
library(RaMS)
library(DBI)
library(RSQLite)
library(data.table)
library(xml2)

turn_mzml_sqlite <- function(ms_files, outfile, ordered=NULL){
  msdata <- grabMSdata(ms_files)
  
  simple_cols <- setdiff(names(msdata$metadata), "config_data")
  file_info <- msdata$metadata[,simple_cols,with = FALSE]
  
  ms_levels <- lapply(ms_files, grabAccessionData, "MS:1000511")
  ms_file_vec <- rep(basename(ms_files), sapply(ms_levels, nrow))
  ms_levels <- do.call(what = "rbind", ms_levels)
  colnames(ms_levels) <- c("name", "ms_level")
  ms_levels$filename <- ms_file_vec
  ms_levels$ms_level <- as.numeric(ms_levels$ms_level)
  ms_levels <- as.data.table(ms_levels)
  
  getScanNum <- function(filename){
    xml_data <- read_xml(filename)
    cont_nums <- xml_attr(xml_find_all(xml_data, "//d1:spectrum"), "id")
    scan_nums <- as.numeric(gsub(".*scan=", "", cont_nums))
    data.table(filename=basename(filename), scan_num=scan_nums)
  }
  scan_nums <- lapply(ms_files, getScanNum)
  scan_nums <- do.call(what = "rbind", scan_nums)
  scan_levels <- cbind(ms_levels, scan_nums[,filename:=NULL])
  scan_levels[,name:=NULL]
  ms1_nums <- scan_levels[ms_level==1,"scan_num"]
  tic <- cbind(ms1_nums, msdata$TIC)
  setnames(tic, c("scan_num", "rt", "TIC", "filename"))
  bpc <- cbind(ms1_nums, msdata$BPC)
  setnames(bpc, c("scan_num", "rt", "BPC", "filename"))
  scan_info <- merge(scan_levels, merge(tic, bpc), all.x = TRUE)
  
  conn <- dbConnect(RSQLite::SQLite(), outfile)
  dbWriteTable(conn, "file_info", file_info)
  dbWriteTable(conn, "scan_info", scan_info)
  dbWriteTable(conn, "MS1", msdata$MS1)
  dbWriteTable(conn, "MS2", msdata$MS2)
  
  if(!is.null(ordered)){
    if(ordered=="mz"){
      dbExecute(conn, "CREATE INDEX mz_idx ON MS1 (mz)")
      dbExecute(conn, "CREATE INDEX premz_idx ON MS2 (premz)")
    }
    if(ordered=="rt"){
      dbExecute(conn, "CREATE INDEX rt_idx ON scan_info (rt)")
      dbExecute(conn, "CREATE INDEX rt_idx ON MS1 (rt)")
      dbExecute(conn, "CREATE INDEX rt_idx ON MS2 (rt)")
    }
  }
  
  dbDisconnect(conn)
  return(outfile)
}

# Example usage:
ms_files <- list.files("demo_data", pattern = "mzML", full.names = TRUE)
turn_mzml_sqlite(ms_files, "demo_data/msdata.sqlite", ordered = "mz")

# Example chromatogram extraction
library(DBI)
pmppm <- function(mass, ppm)c(mass * (1 - ppm/1e+06), mass * (1 + ppm/1e+06))

conn <- dbConnect(RSQLite::SQLite(), "demo_data/msdata.sqlite")
chr_statement <- "SELECT * FROM MS1 WHERE mz BETWEEN ? AND ? ORDER BY filename, rt"
chrom <- dbGetQuery(conn, chr_statement, params=pmppm(118.0865, 10))
dbDisconnect(conn)
plot(chrom$rt, chrom$int, type = "l")
