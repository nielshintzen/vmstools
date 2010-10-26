

 mergedTable2Fishframe <- function(general=general=list(output.path=
             file.path("H:","DIFRES","VMSanalysis","results_merged","DKWaters"),
                                          a.year=2009, a.country="DNK")){
                                          
  ve   <- mergedTable2FishframeVE  (general=general)
  vsl  <- mergedTable2FishframeVSL (general=general)

  # replace NA by " "
  #...
  
  # add fake columns to get the same ncol()
  #...


  # bind and order
  #(to get a VE line and then VSL lines, VE and then VSL lines, etc.)
  ff <- rbind(ve,vsl)
  library(doBy)
  ff <- orderBy(~c_square+LE_MET_level6+month+recordtype, data=ff)


  # save
  write.table(ff, file=file.path(general$output.path,
         paste("ff_", general$a.year, ".txt", sep='')),  dec=".", sep=";",
          quote=FALSE, row.names=FALSE)


  return(ff)
  }


# calls
# vsl <- mergedTable2FishframeVSL (general=list(output.path=file.path("C:","VMSanalysis","results_merged","DKWaters"),
#                                          a.year=2009, a.country="DNK") )
# ve <- mergedTable2FishframeVE (general=list(output.path=file.path("C:","VMSanalysis","results_merged","DKWaters"),
#                                          a.year=2009, a.country="DNK") )
# ff <- mergedTable2Fishframe (general=list(output.path=file.path("C:","VMSanalysis","results_merged","DKWaters"),
#                                          a.year=2009, a.country="DNK") )
