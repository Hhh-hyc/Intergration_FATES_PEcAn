#-------------------------------------------------------------------------------
# Copyright (c) 2016 NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

# R Code to convert NetCDF CF met files into NetCDF FATES met files.

##' met2model wrapper for FATES
##' 
##' @title met2model for FATES
##' @export
##' @param in.path location on disk where inputs are stored
##' @param in.prefix prefix of input and output files
##' @param outfolder location on disk where outputs will be stored
##' @param start_date the start date of the data to be downloaded (will only use the year part of the date)
##' @param end_date the end date of the data to be downloaded (will only use the year part of the date)
##' @param lst timezone offset to GMT in hours
##' @param overwrite should existing files be overwritten
##' @param verbose should the function be very verbosefor(year in start_year:end_year)
##' @importFrom ncdf4 ncvar_get ncdim_def ncatt_get ncvar_put
met2model.FATES <- function(in.path, in.prefix, outfolder, start_date, end_date, lst = 0, lat, lon, 
                            overwrite = FALSE, verbose = FALSE, ...) {
  
  # General Structure- FATES Uses Netcdf so we need to rename vars, split files from years into months, and generate the header file
  # Get Met file from inpath.
  # Loop over years (Open nc.file,rename vars,change dimensions as needed,close/save .nc file)
  # close
  # defining temporal dimension needs to be figured out. If we configure FATES to use same tstep then we may not need to change dimensions  
  
  
  insert <- function(ncout, name, unit, data) {
    var   <- ncdf4::ncvar_def(name = name, units = unit, dim = dim, missval = -6999, verbose = verbose)
    ncout <- ncdf4::ncvar_add(nc = var, v = var, verbose = verbose) ##1. -6999?
    ncvar_put(nc = var, varid = name, vals = data)
    return(invisible(ncout))
  }
  sm <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365) * 86400  ## day of year thresholds
  
  ## Create output directory
  dir.create(outfolder)
  
  # Process start and end dates
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date   <- as.POSIXlt(end_date, tz = "UTC")
  start_year <- lubridate::year(start_date)
  end_year   <- lubridate::year(end_date)
  
  ## Build met
  for (year in start_year:end_year) {
    
    in.file <- file.path(in.path, paste(in.prefix, year, "nc", sep = "."))
    
    if (file.exists(in.file)) {
      
      ## Open netcdf file
      nc <- ncdf4::nc_open(in.file)
      
      ## extract variables. These need to be read in and converted to CLM names (all units are correct)
      time      <- ncvar_get(nc, "time")
      LATIXY    <- ncvar_get(nc, "latitude")
      LONGXY    <- ncvar_get(nc, "longitude")
      FLDS      <- ncvar_get(nc, "surface_downwelling_longwave_flux_in_air")  ## W/m2
      FSDS      <- ncvar_get(nc, "surface_downwelling_shortwave_flux_in_air")  ## W/m2
      PRECTmms  <- ncvar_get(nc, "precipitation_flux")  ## kg/m2/s -> mm/s (same val, diff name)
      PSRF      <- ncvar_get(nc, "air_pressure")  ## Pa
      QBOT      <- ncvar_get(nc, "specific_humidity")  ## g/g -> kg/kg
      TBOT      <- ncvar_get(nc, "air_temperature")  ## K
      WIND      <- sqrt(ncvar_get(nc, "eastward_wind") ^ 2 + ncvar_get(nc, "northward_wind") ^ 2)  ## m/s
      
      ## CREATE MONTHLY FILES
      for (mo in 1:12) {
        tsel <- which(time > sm[mo] & time <= sm[mo + 1]) #time? PEcAn-output consistent
        #define dim
        lat.dim  <- ncdim_def(name = "lat", units = "", vals = 1:1, create_dimvar = FALSE)
        lon.dim  <- ncdim_def(name = "lon", units = "", vals = 1:1, create_dimvar = FALSE)
        time.dim <- ncdim_def(name = "time", units = "seconds", vals = time,
                              create_dimvar = TRUE, unlim = TRUE)#left to CTSM automatically transfer
        scalar.dim <- ncdim_def(name='"scalar", units = "", vals = 1:1, creat_dimvar = FALSE) ## dimensions in forcing data
        dim      <- list(lat.dim, lon.dim, time.dim, scalar.dim)  ## Original question: docs say this should be time,lat,lon but get error writing unlimited first
        ## http://www.cesm.ucar.edu/models/cesm1.2/clm/models/lnd/clm/doc/UsersGuide/x12979.html
        
        # LATITUDE
        var <- ncdf4::ncvar_def(name = "LATIXY", units = "degree_north",
                         dim = list(lat.dim, lon.dim), missval = as.numeric(-9999)) ##2. dim consistent with .nc, missing value==NaNf?
        ncout <- ncdf4::nc_create(outfile, vars = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "LATIXY", vals = LATIXY) #same with FATES

        # LONGITUDE
        var <- ncdf4::ncvar_def(name = "LONGXY", units = "degree_east",
                         dim = list(lat.dim, lon.dim), missval = as.numeric(-9999))
        ncout <- ncdf4::ncvar_add(nc = ncout, v = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "LONGXY", vals = LONGXY)
        
        #TIME
        ncout <- insert(ncout, "time", "days", time) #3. neceaasry?
)
        # EDGEE
        var <- ncdf4::ncvar_def(name = "EDGEE", units = "degrees_east",
                         dim = list(scalar.dim, lat.dim, lon.dim), missval = as.numeric(-9999))
        ncout <- ncdf4::ncvar_add(nc = ncout, v = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "EDGEE", vals = LONGXY+0.005)

        # EDGES
        var <- ncdf4::ncvar_def(name = "EDGES", units = "degrees_north",
                         dim = list(scalar.dim, lat.dim, lon.dim), missval = as.numeric(-9999))
        ncout <- ncdf4::ncvar_add(nc = ncout, v = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "EDGES", vals = LONGXY-0.005)

        # EDGEN
        var <- ncdf4::ncvar_def(name = "EDGEN", units = "degrees_north",
                         dim = list(scalar.dim, lat.dim, lon.dim), missval = as.numeric(-9999))
        ncout <- ncdf4::ncvar_add(nc = ncout, v = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "EDGEN", vals = LONGXY+0.005)
        
        # EDGEW # *4  edge for resolution , edge-central 0.005, # PEcAn provide range of grid?
        var <- ncdf4::ncvar_def(name = "EDGEW", units = "degrees_east",
                         dim = list(scalar.dim, lat.dim, lon.dim), missval = as.numeric(-9999))
        ncout <- ncdf4::ncvar_add(nc = ncout, v = var, verbose = verbose)
        ncvar_put(nc = ncout, varid = "EDGEW", vals = LONGXY-0.005)
        ## saperately create files
        # precipitation
        outfile_prec <- file.path(outfolder, paste0("Prec", formatC(year, width = 4, flag = "0"), "-",
                                               formatC(mo, width = 2, flag = "0"), ".nc")) # 2. name beofre year? or CTSM_container_edit_forcing name with only year&month in user_datm_streams
        if (file.exists(outfile_prec) & overwrite == FALSE) {
          next
        }
        ncout_prec <- ncout 
        ## precipitation_flux
        ncout_prec <- insert(ncout_prec, "PRECTmms", "mm/s", PRECTmms)
        ncdf4::nc_close(ncout_prec)
       
        # solar
        outfile_slr <- file.path(outfolder, paste0("Slr", formatC(year, width = 4, flag = "0"), "-",
                                               formatC(mo, width = 2, flag = "0"), ".nc")) # 2. name beofre year? or CTSM_container_edit_forcing name with only year&month in user_datm_streams
        if (file.exists(outfile_Slr) & overwrite == FALSE) {
          next
        }
        ncout_slr <- ncout
        ## surface_downwelling_shortwave_flux_in_air
        ncout_slr <- insert(ncout_slr, "FSDS", "W m-2", FSDS)
        ncdf4::nc_close(ncout_slr)

        # temper
        outfile_tem <- file.path(outfolder, paste0("Tem", formatC(year, width = 4, flag = "0"), "-",
                                               formatC(mo, width = 2, flag = "0"), ".nc")) # 2. name beofre year? or CTSM_container_edit_forcing name with only year&month in user_datm_streams
        if (file.exists(outfile_Tem) & overwrite == FALSE) {
          next
        }
        ncout_tem <- ncout
        ## surface_downwelling_longwave_flux_in_air
        ncout_tem <- insert(ncout_tem, "FLDS", "W m-2", FLDS)
        ## air_pressure
        ncout_tem <- insert(ncout_tem, "PSRF", "Pa", PSRF)
        ## specific_humidity
        ncout_tem <- insert(ncout_tem, "QBOT", "kg/kg", QBOT)
        ## air_temperature
        ncout_tem <- insert(ncout_tem, "TBOT", "K", TBOT)
        ## eastward_wind & northward_wind
        ncout_tem <- insert(ncout_tem, "WIND", "m/s", WIND)
        ncdf4::nc_close(ncout_tem)
        }
        ncdf4::nc_close(nc)
    }  ## end input file exists
  }  ### end year loop over met files
  
  PEcAn.logger::logger.info("Done with met2model.FATES")
  
  return(data.frame(file = paste0(outfolder, "/"), 
                    host = c(PEcAn.remote::fqdn()), 
                    mimetype = c("application/x-netcdf"), 
                    formatname = c("CLM met"), 
                    startdate = c(start_date), 
                    enddate = c(end_date), 
                    dbfile.name = "", 
                    stringsAsFactors = FALSE))
} # met2model.FATES
