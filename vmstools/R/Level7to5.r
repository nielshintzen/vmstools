# define assemblage of species according to DCR

Level7to5=function(x){
  if(x %in% c('NEP','PRA','CRK','CRB','CRG','CRJ','CRA','CRQ','CRR','KCT','LBA','LBE','AES','PEN','PAN','CRU','ERS','CSH','CRE','KEF','SLO','CPR','LIO','SCR')) x <- "CRU" # "Crustaceans" # FAO sp code
  
  if(x %in% c('BLL','COD','HAD','RED','REG','REB','REN','HKS','HKR','PLA','WIT','YEL','GHL','HAL','FLW','FLS','FLD','FLX','FLY','ANG','SRA','TOM','ANT','CUN','USK','GRC','BLI','LIN','LUM',
              'KGF','PUF','ELZ','OPT','POC','RNG','RHG','SCU','CAT','CAA','CAS','GRO','BSS','GUG','LEM','LEZ','SRA','DAB','FLE','PLE','POK','POL','SAN','SOL','WHG','TUR','APO','ANF','ARU',
              'COE','DGS','DGX','LEF','LEM','MUX','MEG','SBR','SKA','SYC','SYT','SDV','GAG','AGN','JRS','RAJ','JOD','MGS','WEG','CTZ','GUR','SOS','MNZ','BIB','GUU','MGR','MGS','MUR','RJC',
              'RJM','MON')) x <- "DEF" # "Demersal_fish"  # FAO sp code

  if(x %in% c('ABR','HLT','KRS','MUS','KNK','KMS','SNE','MOL','CLR','CLH','CLQ','CLS','CLT','COC','SCB','SCC','ISC','SCA','SCX','OYA','WHX','PER','WHE','QSC','VSC','MYV','VEV','GKL','SCE')) x <- "MOL" # "Molluscs"  # FAO sp code

  if(x %in% c('ELE','STU')) x <- "CAT" # "Catadromus_species"   # FAO sp code

  if (x %in% c('ANE','HER','SPR','MAC','HKE','SAU','ANB','SSM','MHA','BUT','CAP','PIL','JAX','SHG','ARU','ARY','MAS','GAR','MUL','WHB','NOP','HOM')) x <- "SPF" # "Small_pelagic_fish"  # FAO sp code

  if (x %in% c('BLU','CVJ','FRI','KGM','SSM','SAI','WHM','BUM','SWO','ALB','BON','LTA','BET','BFT','SKJ','YFT','POR','DGS','GUQ','DGX','SRX','SKA','SRX','TUN','GUQ','SBG')) x <- "LPF" # "Large_pelagic_fish"  # FAO sp code

  if (x %in% c('CHR','SAL','APU','TRS','TRR','SLZ','TLV','PLN','VIV','LAR','STU','SME','TRO','SHG','SHD','PLN','LAU')) x <- "ANA" # "Anadromous_species"

  if (x %in% c('FBM','FPE','FPI','PLN','FBR','FCP','FTE','BRB','FBR','FBU','FCC','FCP','FPP','FRO','SAR','ELP','FFX')) x <- "FWS" # "Freshwater_species"

  if(x %in% c('SQL','SQI','SQU','OCC','CTC','SQR','IAX','SQZ','OMZ')) x <- "CEP" # "Cephalopods"
  
  if(x %in% c('HYD','ORY')) x <- "DWS" # "Deep_Water_Fish"

  if(x %in% c('MZZ','OTH','FRF','YYY','ZZZ')) x <- "OTH" # "Other"  # FAO sp code
  
  return(x)
}
