Level5=function(x){
  assemblage <- NULL
  if(length(x)==0){   assemblage <- c(assemblage,"NA")
  }else{
    for (i in seq(length(x))) {
      if (names(x[i]) %in% c('NEP','PRA','CRK','CRB','CRG','CRJ','CRA','CRQ','CRR','KCT','LBA','LBE','AES','PEN','PAN','CRU','ERS','CSH','CRE','KEF','SLO','CPR','LIO','SCR'))  assemblage <- c(assemblage,"CRU") else
      if (names(x[i]) %in% c('BLL','COD','HAD','RED','REG','REB','REN','HKS','HKR','PLA','WIT','YEL','GHL','HAL','FLW','FLS','FLD','FLX','FLY','ANG','SRA','TOM','ANT','CUN','USK','GRC','BLI','LIN','LUM',
              'KGF','PUF','ELZ','OPT','POC','RNG','RHG','SCU','CAT','CAA','CAS','GRO','BSS','GUG','LEM','LEZ','SRA','DAB','FLE','PLE','POK','POL','SAN','SOL','WHG','TUR','APO','ANF','ARU',
              'COE','DGS','DGX','LEF','LEM','MUX','MEG','SBR','SKA','SYC','SYT','SDV','GAG','AGN','JRS','RAJ','JOD','MGS','WEG','CTZ','GUR','SOS','MNZ','BIB','GUU','MGR','MGS','MUR','RJC',
              'RJM','MON','RJN','RJE'))  assemblage <- c(assemblage,"DEF") else # "Demersal_fish"  # FAO sp code
      if (names(x[i]) %in% c('ABR','HLT','KRS','MUS','KNK','KMS','SNE','MOL','CLR','CLH','CLQ','CLS','CLT','COC','SCB','SCC','ISC','SCA','SCX','OYA','WHX','PER','WHE','QSC','VSC','MYV','VEV','GKL','SCE','OYF','OYG'))  assemblage <- c(assemblage,"MOL") else # "Molluscs"  # FAO sp code
      if (names(x[i]) %in% c('ELE','STU'))  assemblage <- c(assemblage,"CAT") else # "Catadromus_species"   # FAO sp code
      if (names(x[i]) %in% c('ANE','HER','SPR','MAC','HKE','SAU','ANB','SSM','MHA','BUT','CAP','PIL','JAX','SHG','ARU','ARY','MAS','GAR','MUL','WHB','NOP','HOM'))  assemblage <- c(assemblage,"SPF") else # "Small_pelagic_fish"  # FAO sp code
      if (names(x[i]) %in% c('BLU','CVJ','FRI','KGM','SSM','SAI','WHM','BUM','SWO','ALB','BON','LTA','BET','BFT','SKJ','YFT','POR','DGS','GUQ','DGX','SRX','SKA','SRX','TUN','GUQ','SBG'))  assemblage <- c(assemblage,"LPF") else # "Large_pelagic_fish"  # FAO sp code
      if (names(x[i]) %in% c('CHR','SAL','APU','TRS','TRR','SLZ','TLV','PLN','VIV','LAR','STU','SME','TRO','SHG','SHD','PLN','LAU'))  assemblage <- c(assemblage,"ANA") else # "Anadromous_species"
      if (names(x[i]) %in% c('FBM','FPE','FPI','PLN','FBR','FCP','FTE','BRB','FBR','FBU','FCC','FCP','FPP','FRO','SAR','ELP','FFX'))  assemblage <- c(assemblage,"FWS") else # "Freshwater_species"
      if (names(x[i]) %in% c('SQL','SQI','SQU','OCC','CTC','SQR','IAX','SQZ','OMZ','CTL'))  assemblage <- c(assemblage,"CEP") else # "Cephalopods"
      if (names(x[i]) %in% c('HYD','ORY')) assemblage <- c(assemblage,"DWS") else  # "Deep_Water_Fish"
      if (names(x[i]) %in% c('MZZ','OTH','FRF','YYY','ZZZ'))  assemblage <- c(assemblage,"OTH") # "Other"  # FAO sp code
    }
  }
  assemblage <- do.call(paste,as.list(unique(assemblage)))
  return(assemblage)
}