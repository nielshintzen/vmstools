# Processing of industry data from questionnaire survey regarding vessel and gear specifications

# Introduction #

The R code to proceed the gear data obtained from the questionnaire survey is given as illustration. Equivalent coding is applied for OT, TBB, DS and DRB gear data (with some specifics for each of them). The final product is a table of parameters (estimates\_for\_gear\_param\_per\_metier.txt) defining the relationships to be used in the BENTHIS WP2 workflow.

This table will be disseminated by mails to partners so that they can incorporate it to the BENTHIS WP2 workflow.

Note that some BENTHIS metiers (e.g. OT\_CRU, etc.) are defined on the way which supposes that each EFLALO metier (e.g. OTB\_CRU\_70-89\_2\_35) should be linked to them in advance of using the table in the BENTHIS WP2 workflow. To proceed, it is basically one line of code that each partner should incorporate in their own workflow to convert the eflalo metier into the BENTHIS metier categorization. This could be for example levels(tacsatp$LE\_MET) <- c('OT\_CRU', 'OT\_DEM', 'OT\_MIX', etc.).


# Details #


Hereafter, the R code to produce the table of parameters to be incorporated into the BENTHIS WP2 workflow which will assign a gear width to each metier of each vessel recorded in the tacsatp merged data. The table of parameters (relationships to predict the gear width per BENTHIS metier from the kW or LOA) has been built from the following R code, after collating the outcomes of the analysis split by OT, TBB, DRB and DS for Otter trawl, beam trawl, dredge and seine respectively. Hereafter only the code for OT is providen.

```


# plot (with ggplot2)
library(ggplot2)

# paths
dataPath  <- file.path("C:", "BENTHIS", "data_gear_spec_questionnaire")
outPath   <- file.path("C:", "BENTHIS", "data_gear_spec_questionnaire")



##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-----------------OT------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------

# read file
ind_DK  <- read.table(file= file.path(dataPath, 'OT', 'OT_DK_11062014.csv'), sep=";", header=TRUE )
ind_IRE <- read.table(file= file.path(dataPath, 'OT', 'OT_IRE_11062014.csv'), sep=";", header=TRUE )
ind_SCO <- read.table(file= file.path(dataPath, 'OT', 'OT_SCO_11062014.csv'), sep=";", header=TRUE )
ind_SWE <- read.table(file= file.path(dataPath, 'OT', 'OT_SWE_11062014.csv'), sep=";", header=TRUE )
ind_NL  <- read.table(file= file.path(dataPath, 'OT', 'OT_NL_11062014.csv'), sep=";", header=TRUE )
ind_TUR <- read.table(file= file.path(dataPath, 'OT', 'OT_TUR_11062014.csv'), sep=";", header=TRUE )
ind_MED <- read.table(file= file.path(dataPath, 'OT', 'OT_MED_11062014.csv'), sep=";", header=TRUE ) # italian & spanish data
ind_NOR <- read.table(file= file.path(dataPath, 'OT', 'OT_NOR_11062014.csv'), sep=";", header=TRUE )
ind_BEL <- read.table(file= file.path(dataPath, 'OT', 'OT_BEL_11062014.csv'), sep=";", header=TRUE )
ind_FRA <- read.table(file= file.path(dataPath, 'OT', 'OT_FRA_11062014.csv'), sep=";", header=TRUE )
ind_GRE <- read.table(file= file.path(dataPath, 'OT', 'OT_GRE_11062014.csv'), sep=";", header=TRUE )


# small adjustements for not reliable records
ind_GRE [ind_GRE$Variable.name=="Vessel_LOA", "Value"] <- NA


# collate
cols <- c('Anonymous.vessel_ID','Variable.name','Variable','Value')
ind  <- rbind.data.frame (
cbind(country="DK", ind_DK [, cols]),
cbind(country="IRE", ind_IRE [, cols]),
cbind(country="SCO", ind_SCO [, cols]),
cbind(country="SWE", ind_SWE [, cols]),
cbind(country="NL", ind_NL [, cols]),
cbind(country="TUR", ind_TUR [, cols]),
cbind(country="MED", ind_MED [, cols]),
cbind(country="NOR", ind_NOR [, cols]),
cbind(country="BEL", ind_BEL [, cols]),
cbind(country="FRA", ind_FRA [, cols]),
cbind(country="GRE", ind_GRE [, cols])
)


# explore
head(ind [ind$Variable.name=="Trawl_model",])
levels(ind$Variable.name)
an <- function(x) as.numeric(as.character(x))

# check complete case for some variables (which need to be strictly informed together for a given vessel)
vids_Clp <- unique(as.character(ind [ind$Variable.name=="Clump_weight", "Anonymous.vessel_ID"]))
vids_DoW <- unique(as.character(ind [ind$Variable.name=="Door_weight", "Anonymous.vessel_ID"]))
vids_to_be_removed <- vids_Clp[!(vids_Clp %in% vids_DoW)]
vids_to_be_removed <- c(vids_to_be_removed, vids_DoW[!(vids_DoW %in% vids_Clp)])
#=> if any, then remove the concerned vessel.
ind <- ind [!ind$Anonymous.vessel_ID %in% vids_to_be_removed, ]




CT   <- ind [ind$Variable.name=="Consumption_trawling", "Value"]
CS   <- ind [ind$Variable.name=="Consumption_steaming", "Value"]
Str  <- ind [ind$Variable.name=="Speed_trawling", "Value"]
kW   <- ind [ind$Variable.name=="Vessel_kW", "Value"]
LOA  <- ind [ind$Variable.name=="Vessel_LOA", "Value"]
DoS  <- ind [ind$Variable.name=="Door_spread", "Value"]
DoW  <- ind [ind$Variable.name=="Door_weight", "Value"]
DoN  <- ind [ind$Variable.name=="Door_number", "Value"]
DoL  <- ind [ind$Variable.name=="Door_length", "Value"]
OcW  <- ind [ind$Variable.name=="Otherchain_weight", "Value"]
OcN  <- ind [ind$Variable.name=="Otherchain_number", "Value"]
TiW  <- ind [ind$Variable.name=="Ticklerchain_weight", "Value"]
TiN  <- ind [ind$Variable.name=="Ticklerchain_number", "Value"]
Clp  <- ind [ind$Variable.name=="Clump_weight", "Value"]
GrW  <- ind [ind$Variable.name=="Groundgear_weight", "Value"]
GrL  <- ind [ind$Variable.name=="Groundgear_length", "Value"]
sps  <- ind [ind$Variable.name=="Targetspecies_single", "Value"]
sp1  <- ind [ind$Variable.name=="Primarytarget_mixed", "Value"]
sp2  <- ind [ind$Variable.name=="Secondarytarget_mixed", "Value"]
sp3  <- ind [ind$Variable.name=="Thirdtarget_mixed", "Value"]
bt   <- ind [ind$Variable.name=="Bottom_type", "Value"]
nbt  <- ind [ind$Variable.name=="Trawl_number", "Value"]
sptr <- ind [ind$Variable.name=="Speed_trawling", "Value"]
spst <- ind [ind$Variable.name=="Speed_steaming", "Value"]
ctry <- ind [ind$Variable.name=="Door_spread", "country"]
vid  <- ind [ind$Variable.name=="Door_spread", "Anonymous.vessel_ID"]
area <- ind [ind$Variable.name=="Fishing_area", "Value"]
mesh <- ind [ind$Variable.name=="Codend_meshsize", "Value"]
swee <- ind [ind$Variable.name=="Sweerp_length", "Value"]
brdl <- ind [ind$Variable.name=="Bridle_length", "Value"]




# intermediate calculation
# total gear weight
dd    <- rbind(c(an(DoW) *  an(DoN) + an(GrW) * an(nbt)), an(Clp), c(an(OcW) * an(OcN)))
TotGW <- apply(dd,2, sum, na.rm=TRUE) # DoW and GrW need complete cases here.
TotGW[is.na(dd[1,])] <- NA

# intermediate calculation for getting a percentage per gear sub-component
brdl <- an(brdl)
brdl[is.na(an(brdl))] <-  mean(an(brdl), na.rm=T)  # ASSUMPTION TO FILL IN THE GAP


pw_swee_brdl        <- sin(13*pi/180)*(an(swee)+an(brdl))*2*an(nbt)   # 13 degrees
pw_grg_brdl         <- an(GrL)*0.4*an(nbt)
pw_doors_clump      <- an (DoL)*0.4*an(DoN)+(0.75*(an(nbt)-1))
pw_sum              <- pw_swee_brdl+ pw_grg_brdl + pw_doors_clump



# check - observations should be lying on the line
plot(an(DoS), an(pw_sum))
lines(1:300, 1:300)


# then, collate:
df1  <- cbind.data.frame(ctry, vid, CT, CS, sptr, spst, kW, LOA, DoS, DoL, DoW, DoN, GrW, nbt, Clp, OcW, OcN, TotGW, GrL, sps, sp1, sp2, sp3, bt, area, mesh, swee, brdl,
pw_swee_brdl, pw_grg_brdl, pw_doors_clump, pw_sum)

df1$pw_grg <- an(df1$pw_grg_brdl)
df1$pw_grg[is.na(df1$pw_grg)] <-  an(df1$DoS[is.na(df1$pw_grg_brdl)])-  an(df1$pw_swee_brdl[is.na(df1$pw_grg_brdl)]) - an(df1$pw_doors_clump[is.na(df1$pw_grg_brdl)])   # ASSUMPTION TO FILL IN THE GAP



# look at a potential proxy in case the door spread is not informed....
plot(an(df1$DoS),an(df1$swee), col=df1$ctry, pch=16)



# shortcut to retrieve a given Value
idx <- ind [ind$Variable.name=="Fishing_area", "Value"]==2.2
dd  <- ind [ind$Variable.name=="Fishing_area", ]
dd[idx,]

idx <- ind [ind$Variable.name=="Door_spread", "Value"]==0
dd  <- ind [ind$Variable.name=="Door_spread", ]
dd[idx,]


# refactor some variables (if needed)
df1$LOA_class     <- cut(an(df1$LOA), breaks=c(0,12,18,24,40,100))
df1$sptrawl       <- df1$sptr  # sauv trawling speed
df1$sptr          <- cut(an(df1$sptr), breaks= seq(0,10,0.75))     # trawling speed
df1$mesh_class    <- cut(an(df1$mesh), breaks= c(0,90,220))        # codend mesh size
df1$GrW_class     <- cut(an(df1$GrW), breaks= seq(0,1600, 200) )   # ground gear weight


# DCF metier coding
df1$metier <- as.character(df1$sps) # init
df1[df1$metier %in% c('NEP ','NEP','PRA','Nephrops','Nephrops trawl',
'TGS', 'ARA','DPS' ), 'metier'] <- 'OT_CRU'
df1[df1$metier %in% c('COD','PLE','SOL', 'LEM', 'WHG', 'WHI', 'POK',
'PDS','HAD','had','HKE','MON', 'MUT',
'NOP'), 'metier'] <- 'OT_DMF'
df1[df1$metier %in% c('SAN','SPR','CAP'), 'metier'] <- 'OT_SPF'
df1[df1$metier %in% c('NR', "NI",'0','','ni'), 'metier'] <- 'OT_MIX'
df1$metier <- as.factor(df1$metier)


# deal with mixed fisheries
df1$metier <- as.character(df1$metier) # init
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('NEP ','NEP',
'Nephrops','Nephrops trawl', 'PRA', 'CSH'
), 'metier'] <- 'OT_MIX_CRU_DMF'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('TGS') &  df1$sp2 %in% c('OCC'), 'metier'] <- 'OT_MIX_TGS_OCC'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('TGS') &  df1$sp2 %in% c('CTC'), 'metier'] <- 'OT_MIX_TGS_CTC'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('ARA'
), 'metier'] <- 'OT_MIX_ARA'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('DPS'
), 'metier'] <- 'OT_MIX_CRU'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('PLE','SOL', 'LEM',
'MON', 'MUT') & df1$sp1_type=="DMF_BEN", 'metier'] <- 'OT_MIX_DMF_BEN'
df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('COD','PLE','SOL', 'LEM',
'WHG', 'WHI', 'POK', 'PDS','HAD','had','HKE',
'MUT', 'NOP')& df1$sp1_type=="DMF_PEL", 'metier'] <- 'OT_MIX_DMF_PEL'

df1[df1$metier %in% c('OT_MIX') & df1$sp1 %in% c('SAN',
'SPR'), 'metier'] <- 'OT_MIX_SPF'
df1[df1$ctry %in% c('TUR'), 'metier'] <- 'OT_MIX'
df1$metier <- as.factor(df1$metier)



# nb of observations per country
table(df1$ctry)

# nb of observations per metier
table(df1$metier)

# nb of observations per country per metier
table(df1$ctry,df1$metier)

# nb of observations per country per metier
table(df1$ctry,df1$metier3)

##-------------------------------------------------
##-------------------------------------------------
# compute percentages per gear subcomponent
##-------------------------------------------------
##-------------------------------------------------

table(df1$ctry, is.na(df1$pw_sum))

table(df1$ctry, is.na(df1$pw_swee_brdl))

table(df1$ctry, is.na(df1$pw_grg))

table(df1$ctry, is.na(df1$pw_doors_clump))


table(df1$metier, is.na(df1$pw_sum))

plot(an(df1$swee),an(df1$brdl), col=df1$metier, pch=16)

table(df1$pw_grg <0)


# identify the non-single trawl
df1$metier_asterik <- as.character(df1$metier) # init
##idx          <-  an(df1$nbt)>=2
##idx[is.na(idx)] <- FALSE
##df1$metier_asterik[idx] <-paste(df1$metier_asterik[idx], "*", sep="")
df1$metier_asterik <- as.factor(df1$metier_asterik)




df_subcomponents <- df1[, c('metier_asterik','pw_swee_brdl','pw_grg', 'pw_doors_clump', 'DoS')]

df_subcomponents$pw_swee_brdl_percent   <- an(df_subcomponents$pw_swee_brdl)/an(df_subcomponents$DoS)*100
df_subcomponents$pw_grg_percent         <- an(df_subcomponents$pw_grg)/an(df_subcomponents$DoS)*100
df_subcomponents$pw_doors_clump_percent <- an(df_subcomponents$pw_doors_clump)/an(df_subcomponents$DoS)*100

df_subcomponents_cleaned <- df_subcomponents[
an(df_subcomponents$pw_doors_clump_percent) >0 &
an(df_subcomponents$pw_doors_clump_percent)<100 &
an(df_subcomponents$pw_swee_brdl_percent) >0 &
an(df_subcomponents$pw_swee_brdl_percent)<100 &
an(df_subcomponents$pw_grg_percent) >0 &
an(df_subcomponents$pw_grg_percent)<100 ,
]
df_subcomponents_cleaned <- df_subcomponents_cleaned[!is.na(df_subcomponents_cleaned[,1]),]
# => a lot of records are unfortunately discarded...but because eg Door spread estimate is smaller than a subcomponent which has no solution!
# note that the number of discarded records also depends very much on the chosen bridle angle (10, 13 or 15 degrees)

subcomponents_matrix <- cbind.data.frame(
av_pw_doors_clump_percent= tapply(df_subcomponents_cleaned$pw_doors_clump_percent, list(df_subcomponents_cleaned$metier), mean),
av_pw_swee_brdl_percent=   tapply(df_subcomponents_cleaned$pw_swee_brdl_percent, list(df_subcomponents_cleaned$metier), mean),
av_pw_grg_percent=         tapply(df_subcomponents_cleaned$pw_grg_percent, list(df_subcomponents_cleaned$metier), mean),
sd_pw_doors_clump_percent= tapply(df_subcomponents_cleaned$pw_doors_clump_percent, list(df_subcomponents_cleaned$metier), sd),
sd_pw_swee_brdl_percent=   tapply(df_subcomponents_cleaned$pw_swee_brdl_percent, list(df_subcomponents_cleaned$metier), sd),
sd_pw_grg_percent=         tapply(df_subcomponents_cleaned$pw_grg_percent, list(df_subcomponents_cleaned$metier), sd),
nb_observations=           tapply(df_subcomponents_cleaned$pw_doors_clump_percent, list(df_subcomponents_cleaned$metier), length)
)

# rescaling to 100%
a_sum <- apply(subcomponents_matrix[,c('av_pw_doors_clump_percent','av_pw_swee_brdl_percent', 'av_pw_grg_percent')], 1, sum)
subcomponents_matrix <- cbind.data.frame( subcomponents_matrix,
pathwidth_doors=   subcomponents_matrix[,'av_pw_doors_clump_percent'] /  a_sum*100,
pathwidth_sweep=   subcomponents_matrix[,'av_pw_swee_brdl_percent']   /  a_sum*100,
pathwidth_grounds= subcomponents_matrix[,'av_pw_grg_percent']         /  a_sum*100
)

# keep the metier of interest only
subcomponents_matrix <- subcomponents_matrix[c("OT_CRU", "OT_CRU*", "OT_DMF", "OT_DMF*", "OT_MIX", "OT_MIX*",
"OT_MIX_DMF_BEN", "OT_MIX_DMF_BEN*", "OT_MIX_DMF_PEL", "OT_MIX_DMF_PEL*", "OT_MIX_CRU", "OT_MIX_CRU_DMF", "OT_MIX_CRU_DMF*", "OT_SPF"),]

# export
write.table(subcomponents_matrix, file=file.path(outPath, "estimates_for_OT_subcomponents_matrix_6Jan15_13degrees.txt"), quote=FALSE)
# => for using this table in the workflow, partners should link each logbooks records to the metier categories found in that table...



##-------------------------------------------------
##-------------------------------------------------
## Door spread vs. kW
##-------------------------------------------------
##-------------------------------------------------

# look at the representativity....
df1$informedDoS <- ifelse(is.na(an(df1$DoS)) | an(df1$DoS)==0,0,1)
table(df1$ctry, df1$informedDoS)

# look at the representativity....
df1$informedarea <- ifelse(as.character(df1$area)=='0',0,1)
table(df1$ctry, df1$informedarea)


# look at the representativity....
table(df1$ctry, df1$metier, df1$informedDoS)

# filter out the 0
df1_DoS <- df1[!is.na(df1$DoS) & df1$DoS!=0,]


# look at the representativity....
df1$informedDoS <- ifelse(is.na(an(df1$DoS)) | an(df1$DoS)==0,0,1)
table(df1$ctry, df1$informedDoS)

# look at the representativity....
df1$informedarea <- ifelse(as.character(df1$area)=='0',0,1)
table(df1$ctry, df1$informedarea)


# look at the representativity....
table(df1$ctry, df1$metier, df1$informedDoS)

# filter out the 0
df1_DoS <- df1[!is.na(df1$DoS) & df1$DoS!=0,]


##-----------------------------------------
##-----------------------------------------
## Export a table of parameters
## for the Benthis vmstools R workflow
## after having made a choice of the most relevant
## categorisation
##-----------------------------------------
##-----------------------------------------

# get the coeffs for DoS~ LOA or kW
coeffs_DoS <- NULL
df1_DoS$DoS     <-  an(df1_DoS$DoS)
df1_DoS$LOA     <-  an(df1_DoS$LOA)
df1_DoS$kW      <-  an(df1_DoS$kW)
df1_DoS_c       <-  df1_DoS[!is.na(df1_DoS$LOA) & !is.na(df1_DoS$kW),]  # caution: keep complete cases for a true model comparison...

for (a_metier in c("OT_CRU", "OT_DMF", "OT_MIX",  "OT_MIX_DMF_BEN", "OT_MIX_DMF_PEL","OT_MIX_CRU", "OT_MIX_CRU_DMF","OT_SPF")){
a_nls_kW        <- nls(DoS~a*(kW^b), start=list(a=1,b=1),data=df1_DoS_c[df1_DoS_c$metier==a_metier,])
a_lm_kW         <- nls(DoS~a*kW+b, start=list(a=1,b=1),data=df1_DoS_c[df1_DoS_c$metier==a_metier,])
a_nls_LOA       <- nls(DoS~a*(LOA^b), start=list(a=50,b=1),data=df1_DoS_c[df1_DoS_c$metier==a_metier,])
a_lm_LOA        <- nls(DoS~a*LOA+b, start=list(a=50,b=1),data=df1_DoS_c[df1_DoS_c$metier==a_metier,])

#compare goodness of fit
residualSum    <- anova (a_nls_LOA, a_lm_LOA, a_nls_kW, a_lm_kW)
what_is_chosen <- c('a_nls_LOA','a_lm_LOA','a_nls_kW','a_lm_kW') [which.min(residualSum[,2])]
print(what_is_chosen)

nb_records      <- nrow(df1_DoS[df1_DoS$metier==a_metier,])

# then choose the best model....
# (and re-run on the full dataset)
if(what_is_chosen=="a_nls_LOA") {
a_nls_LOA      <- nls(DoS~a*(LOA^b), start=list(a=1,b=1), data=df1_DoS[df1_DoS$metier==a_metier,])  # redo with all the available data
coeffs_DoS     <- rbind.data.frame (coeffs_DoS,  cbind.data.frame(a_metier, param=c('a','b'), summary(a_nls_LOA)$coeff, equ="DoS=a*(LOA^b)", nb_records= nb_records))
}
if(what_is_chosen=="a_lm_LOA"){
a_lm_LOA        <- nls(DoS~a*LOA+b, start=list(a=50,b=1),data=df1_DoS[df1_DoS$metier==a_metier,])
coeffs_DoS      <- rbind.data.frame (coeffs_DoS,  cbind.data.frame(a_metier, param=c('a','b'), summary(a_lm_LOA)$coeff, equ="DoS=(a*LOA)+b", nb_records= nb_records))
}
if(what_is_chosen=="a_nls_kW"){
a_nls_kW        <- nls(DoS~a*(kW^b), start=list(a=1,b=1),data=df1_DoS[df1_DoS$metier==a_metier,])
#confint(a_nls_kW) in MASS
coeffs_DoS      <- rbind.data.frame (coeffs_DoS,  cbind.data.frame(a_metier, param=c('a','b'), summary(a_nls_kW)$coeff, equ="DoS=a*(kW^b)", nb_records= nb_records))
}
if(what_is_chosen=="a_lm_kW"){
a_lm_kW         <- nls(DoS~a*kW+b, start=list(a=1,b=1),data=df1_DoS[df1_DoS$metier==a_metier,])
coeffs_DoS      <- rbind.data.frame (coeffs_DoS,  cbind.data.frame(a_metier, param=c('a','b'), summary(a_lm_kW)$coeff, equ="DoS=(a*kW)+b", nb_records= nb_records))
}
}
rownames(coeffs_DoS) <- NULL


# export
write.table(coeffs_DoS, file=file.path(outPath, "estimates_for_OT_nls_DoS_vs_LOA_or_kW_per_metier_10Oct14.txt"))
# => for using this table in the workflow, partners should link each logbooks records to the metier categories found in that table...



a_nls_kW        <- nls(DoS~a*(kW^b), start=list(a=1,b=1),data=df1_DoS[df1_DoS$metier==a_metier,])
kW              <- seq(0, 10000, by=50)
an_equation     <- as.character(coeffs_DoS[coeffs_DoS$a_metier==a_metier,][1,'equ'])
a               <- coeffs_DoS[coeffs_DoS$a_metier==a_metier & coeffs_DoS$param=='a', 'Estimate']
b               <- coeffs_DoS[coeffs_DoS$a_metier==a_metier  & coeffs_DoS$param=='b', 'Estimate']
plot(df1_DoS[df1_DoS$metier==a_metier, "kW",], df1_DoS[df1_DoS$metier==a_metier, "DoS"], pch=16, col=df1_DoS[df1_DoS$metier==a_metier, "the_colors"], xlab="kW", ylab="Door spread (metre)", axes=FALSE, ylim=c(0,400))
lines(kW, eval(parse(text= an_equation)))
axis(1)
axis(2)

source("C:/BENTHIS/benthis_Rworkflows/predictNLS.R")
c.lim<- predictNLS (a_nls_kW, newdata= data.frame(kW=kW), level = 0.95, nsim = 10000)
lines(cbind(kW,c.lim[,'fit']), col="red", lty="dashed")
lines(cbind(kW,c.lim[,'2.5%']), col="blue", lty="dashed")
lines(cbind(kW,c.lim[,'97.5%']), col="green", lty="dashed")
polygon(x=c(kW, rev(kW)),y=c(c.lim[,'97.5%'], rev(c.lim[,'2.5%'])), border=grey(0.9), col=grey(0.9))


##-------------------------------------------------
##-------------------------------------------------
## do the correponding plot
##-------------------------------------------------
##-------------------------------------------------


tiff(filename = file.path(outPath, paste("plot_estimates_for_OT_nls_DoS_vs_LOA_or_kW_per_metier.tiff",sep="")),
width = 1600, height = 2000,   compression = c("lzw"),
units = "px", pointsize = 12,  res=300)   # high resolution plot
#windows(7,7)
par(mfrow=c(3,3))
par(oma=c(4,4,1,1))
par(mar=c(4,0,2,1))
df1_DoS$the_colors   <- df1_DoS$ctry
rgb.palette                 <- colorRampPalette(c("green", "blue", "red"),
space = "Lab")
levels(df1_DoS$the_colors)  <- rgb.palette (length(unique(df1_DoS$the_colors)))

levels(df1_DoS$the_colors)  <- c( "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",  "grey", "#FFFF99")[1:length(levels(df1_DoS$the_colors) )]

df1_DoS$the_colors          <- as.character(df1_DoS$the_colors)




count <-1
for(met in c("OT_CRU", "OT_DMF", "OT_MIX", "OT_MIX_DMF_BEN", "OT_MIX_DMF_PEL", "OT_MIX_CRU", "OT_MIX_CRU_DMF", "OT_SPF")){

an_equation <- as.character(coeffs_DoS[coeffs_DoS$a_metier==met,][1,'equ'])


if(length(grep('LOA' , an_equation))>0){
range_LOA <- range(df1_DoS[df1_DoS$metier==met, "LOA",], na.rm=TRUE)
LOA       <- seq(range_LOA[1]-10, range_LOA[2]+10, by=1)

# polygon
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']  + 1.96*coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='b', 'Std. Error']
up <- eval(parse(text= an_equation))
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']  - 1.96*coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='b', 'Std. Error']
down <- eval(parse(text= an_equation))


# median line
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']


plot(df1_DoS[df1_DoS$metier==met, "LOA",], df1_DoS[df1_DoS$metier==met, "DoS"], pch=16,
col=df1_DoS[df1_DoS$metier==met, "the_colors"], xlab="LOA (metre)", ylab="Door spread (metre)", axes=FALSE, ylim=c(0,300))
if(an_equation=="DoS=a*(LOA^b)"){
LOAs           <- LOA
a_nls_LOA      <- nls(DoS~a*(LOA^b), start=list(a=1,b=1), data=df1_DoS[df1_DoS$metier==met,])  # redo with all the available data
c.lim          <- predictNLS (a_nls_LOA, newdata= data.frame(LOA=LOAs), level = 0.95, nsim = 10000)
polygon(x=c(LOAs, rev(LOAs)),y=c(c.lim[,'97.5%'], rev(c.lim[,'2.5%'])), border=grey(0.9), col=grey(0.9))
#lines(LOAs, c.lim[,'median'], col=grey(0.5))
}else{
polygon(x=c(LOA,rev(LOA)), y=c(up, rev(down)), border=grey(0.9), col=grey(0.9))
}
points(df1_DoS[df1_DoS$metier==met, "LOA",], df1_DoS[df1_DoS$metier==met, "DoS"], pch=16,
col=df1_DoS[df1_DoS$metier==met, "the_colors"], xlab="LOA (metre)", ylab="Door spread (metre)", axes=FALSE, ylim=c(0,300))
lines(LOA, eval(parse(text= an_equation)))

}


if(length(grep('kW' , an_equation))>0){
range_kW <- range(df1_DoS[df1_DoS$metier==met, "kW",], na.rm=TRUE)
kW       <- seq(range_kW[1]-100, range_kW[2]+100, by=1)
kW       <- kW[kW>0]

# polygon
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']  +  1.96*coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='b', 'Std. Error']
up <- eval(parse(text= an_equation))
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']  -  1.96*coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='b', 'Std. Error']
down <- eval(parse(text= an_equation))


# median line
a           <- coeffs_DoS[coeffs_DoS$a_metier==met & coeffs_DoS$param=='a', 'Estimate']
b           <- coeffs_DoS[coeffs_DoS$a_metier==met  & coeffs_DoS$param=='b', 'Estimate']


plot(df1_DoS[df1_DoS$metier==met, "kW",], df1_DoS[df1_DoS$metier==met, "DoS"],
pch=16, col=df1_DoS[df1_DoS$metier==met, "the_colors"], xlab="kW", ylab="Door spread (metre)", axes=FALSE, ylim=c(0,300))
if(an_equation=="DoS=a*(kW^b)"){
kWs            <- kW
a_nls_kW       <- nls(DoS~a*(kW^b), start=list(a=1,b=1),data=df1_DoS[df1_DoS$metier==met,])
c.lim          <- predictNLS (a_nls_kW, newdata= data.frame(kW=kWs), level = 0.95, nsim = 10000)
polygon(x=c(kWs, rev(kWs)),y=c(c.lim[,'97.5%'], rev(c.lim[,'2.5%'])), border=grey(0.9), col=grey(0.9))
#lines(kW, c.lim[,'median'], col=grey(0.5))
#=> but we got back non-parametric curves!...hard to use in the work flow....
} else{
polygon(x=c(kW,rev(kW)), y=c(up, rev(down)), border=grey(0.9), col=grey(0.9))
}
points(df1_DoS[df1_DoS$metier==met, "kW",], df1_DoS[df1_DoS$metier==met, "DoS"], pch=16, col=df1_DoS[df1_DoS$metier==met, "the_colors"],
xlab="kW", ylab="Door spread (metre)", axes=FALSE, ylim=c(0,300))
lines(kW, eval(parse(text= an_equation)))

}

axis(1)
if(count==1 || count==4 || count==7 || count==10) axis(2, las=2)
box()
title(met, cex.main=0.85)
count <- count+1
}
plot(0,0,type="n", axes=FALSE, xlab="", ylab="")
legend("topright", legend=unique(df1_DoS$ctry), fill=unique(df1_DoS$the_colors), bty="n", cex=1.0, ncol=2)

mtext(side=2, text="Door spread (metre)", line=2.8, outer=TRUE)
#savePlot(file=file.path(outPath, "plot_estimates_for_OT_nls_DoS_vs_LOA_or_kW_per_metier.png"), type="png")

dev.off()






```