# Practical 3: Linking VMS and logbook data and exploring the benefits #

## Introduction ##
You have already come a long way by importing your VMS and logbook data into R, cleaning them to ensure that the basis of any analysis you do is solid. We can't stress it enough, as with normal statistics, that getting your data ready, and in the right format, has a crucial bearing for accurate results.

There is great potential utility in having these two datasets combined. Here we will show how it's done in VMStools, and how to make some first pictures and maps out of the results. You might be disappointed, however, when you see how much work it takes :-).

## Getting the data ready ##
If you still have an R version running with the results of the cleaning step, use those data here. If you had to start R anew, don't worry, we will show you how it works, and with your own data you can do it properly in one go (but note that you should only start to link VMS and logbook data after the cleaning process!).

## Linking VMS and logbook data ##
```

#If you save the tacsatClean and eflaloClean let's load them first
load("D:/output/tacsatClean.RData")
load("D:/output/eflaloClean.RData")

#We assume you still have got a good cleaned set of eflalo and tacsat available, if not, reload the data (data(tacsat);data(eflalo)), and continue without prior cleaning for example purposes
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)

head(tacsatp)
>         VE_COU VE_REF  SI_LATI  SI_LONG    SI_DATE  SI_TIME SI_SP SI_HE            SI_DATIM INTV FT_REF
662894 Atlantis     10 51.67319 3.176576 05/05/1800  4:54:00     2    43 1800-05-05 04:54:00  116 272058
662895 Atlantis     10 51.69432 3.229084 05/05/1800  6:48:00     2   294 1800-05-05 06:48:00  114 272058
662896 Atlantis     10 51.65607 3.160989 05/05/1800  8:44:00     6   144 1800-05-05 08:44:00  116 272058
662897 Atlantis     10 51.60704 3.083506 05/05/1800 10:40:00     7   248 1800-05-05 10:40:00  116 272058
662898 Atlantis     10 51.57692 2.987796 05/05/1800 14:30:00     6   228 1800-05-05 14:30:00  230 272058
662899 Atlantis     10 51.54655 2.876395 05/05/1800 16:24:00     7    50 1800-05-05 16:24:00  114 272058

print(length(which(tacsatp$FT_REF == 0))) #number of records not able to link
```

What we can observe is that another column has been added to the 'tacsatp' object; namely the 'FT\_REF' column. What happened 'under the hood' is that the VMS and logbook data have been linked, based on vessel name and year, and the time stamp, such that each VMS record or ping falls within the time frame of the departure and arrival of a logbook trip for that same vessel. Now we know which VMS records can be associated with which eflalo records. This now allows us to do some new types of analyses. Note that if no linking was possible, FT\_REF equals "0" and this happens in many instances, e.g. when vessels are in harbour.

```

#Copy the gear type from eflalo to tacsat using match.
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#Let's see if it worked correctly
subset(tacsatp,FT_REF == "326494")
>         VE_COU VE_REF  SI_LATI  SI_LONG    SI_DATE  SI_TIME SI_SP SI_HE            SI_DATIM INTV FT_REF LE_GEAR
1019702 Atlantis    102 51.63328 3.480378 20/12/1801  7:16:00   0.6   354 1801-12-20 07:16:00  114 326494     DRB
1019703 Atlantis    102 51.64547 3.462132 20/12/1801  9:12:00   0.0   258 1801-12-20 09:12:00  116 326494     DRB
1019704 Atlantis    102 51.63715 3.447098 20/12/1801 11:06:00   0.0   280 1801-12-20 11:06:00  114 326494     DRB
1019705 Atlantis    102 51.63734 3.451090 20/12/1801 13:02:00   0.2    86 1801-12-20 13:02:00  116 326494     DRB
1019706 Atlantis    102 51.63760 3.458493 20/12/1801 14:58:00   0.6   132 1801-12-20 14:58:00  116 326494     DRB
1019707 Atlantis    102 51.64877 3.528120 20/12/1801 16:52:00   6.8    72 1801-12-20 16:52:00  114 326494     DRB

subset(eflalo,FT_REF == "326494")
>      VE_REF VE_FLT   VE_COU VE_LEN VE_KW VE_TON FT_REF FT_DCOU FT_DHAR    FT_DDAT FT_DTIME FT_LCOU FT_LHAR    FT_LDAT FT_LTIME
48267    102    DRB Atlantis     29   127     NA 326494     nld      NL 20/12/1801 06:00:00     nld      NL 20/12/1801 18:00:00
```

### Exercise 1 ###
  1. Check if, in the example above, the merging went well. How do you check?
  1. What types of gears are available in the VMS dataset, how many of each?

## Exploring the benefits ##
Similar to the presentation, we will quickly look into four types of benefits
  1. Distinguish fishery types
  1. Define activity
  1. Fine scale effort patterns
  1. Distribution of catches / values

### Distinguish fishery types ###
Most often, fishery type is determined by the gear used (mesh size etc.) when fishing. As we have already assigned gear type from eflalo to tacsat, splitting the tacsat data up into bits per gear type is not too difficult.
```

#Splitting the tacsatp dataset into smaller chunks, a TBB and an OTB part
TBB <- subset(tacsatp,LE_GEAR == "TBB")
OTB <- subset(tacsatp,LE_GEAR == "OTB")

#Let's see what speed profiles are associated with these fishery types
par(mfrow=c(2,1))
hist(TBB$SI_SP,breaks=100,xlim=c(0,20),xlab="Knots",main="TBB")
hist(OTB$SI_SP,breaks=100,xlim=c(0,20),xlab="Knots",main="OTB")
```

Do not close this picture immediately, it comes in handy in the next example.

### Define activity ###
In many cases, activity can be defined based on speed. We've just made a picture of the speed profiles. Can you distinguish different behaviors by looking at these plots already? In general, there are 3 types of behavior to be distinguished (but with more information on the fishery, you could potentially indicate even more behaviours).
  1. Floating / being in harbour, associated with very low to zero speeds.
  1. Fishing, maybe towing a trawl, associated with intermediate speeds.
  1. Steaming, associated with high speeds usually going home.

If we now look back at the figure we just made, in the OTB panel, we can observe 2 'peaks'; one with a mean of around 3 knots ranging between 2 and 6-7 knots and one with a mean around 10 knots ranging between 7 and 13 knots. Everything below 1 knot can be assigned to 'floating' / being in harbour. All in all, this is a simple approach, but might still be very time consuming.

Therefore, based on these ideas, different functionalities have been developed in VMStools.
  1. segmentedTacsatSpeed, segmented regression on speed profile and automatic detection of fishing versus no-fishing
  1. activityTacsat, normal density curve fitting based on speed profile but informed by user decisions

We present short examples of both methods.

#### Segmented Tacsat Speed ####
Note that you can only run this example with 1 year of data!
```

#I want to keep track of the number of rows, so I add a unique number to each row
tacsatp$idx       <- 1:nrow(tacsatp)

require(segmented) #if not installed yet: install.packages("segmented",repos=getOption("repos"))
tacsatSegment     <- segmentedTacsatSpeed(tacsatp[which(tacsatp$LE_GEAR=="TBB" & year(tacsatp$SI_DATIM)==1800),],units="year",analyse.by="VE_REF",speed="instantanious")

#Merge the results of segmentTacsatSpeed with the original tacsatp dataset
tacsatp           <- merge(tacsatp,tacsatSegment[,c("idx","SI_STATE")],by=c("idx"),all.x=T)

#Lets look at the results, left = no-fishing, right = fishing
spds <- table(tacsatp$SI_SP,tacsatp$SI_STATE)
barplot(spds,names.arg=rep(dimnames(spds)[[1]],2),beside=T)
#-That looks a bit odd
tacsatp <- calculateSpeed(tacsatp,level="trip",fill.na=T)
par(mfrow=c(2,1))
hist(subset(tacsatp,SI_STATE=="nf")$SI_SPCA,breaks=100,main="Speeds associated with no fishing",xlab="Knots")
hist(subset(tacsatp,SI_STATE=="f" )$SI_SPCA,breaks=100,main="Speeds associated with fishing",xlab="Knots")


```

#### activityTacsat ####
This functionality works on more than 1 year and can be used to define activity of all similar gears or on a vessel level. In its most sophisticated form, it requires input from the user to determine peaks of fishing and non-fishing which will eventually result in better results.

```

#In order to get this working, I cannot have NA's as gear type, so we can assign them as 'NO_GEAR' by:
tacsatp$LE_GEAR[which(is.na(tacsatp$LE_GEAR)==T)] <- "NO_GEAR"

#Analyse a selection of the data, a window will pop-up!
require(tcltk) #- install.packages("tcltk",repos=getOption("repos"))
storeScheme   <- activityTacsatAnalyse(subset(tacsatp,LE_GEAR == "OTM" & year(SI_DATIM) == 1801),units="year",analyse.by="LE_GEAR")
```
Ok, have a look at the storeScheme first. The design is not so complicated, and especially if you have to analyse many vessels, it would pay off to create the 'storeScheme' yourself using e.g. 'expand.grid' functionality in R.
```

#If this gives an error, you need to install the package first, go to R -> packages -> install packages
require(mixtools) #- install.packages("mixtools",repos=getOption("repos"))

#After analysing the peaks, its time to let the magic fit normal density distributions to the data
res           <- activityTacsat(subset(tacsatp,LE_GEAR == "OTM" &
year(SI_DATIM) == 1801),units="year",analyse.by="LE_GEAR",storeScheme)

#Make sure the same selection is used when writing the results to the original tacsat dataset
tacsatSubset <- subset(tacsatp,LE_GEAR == "OTM" & year(SI_DATIM) == 1801 & is.na(SI_SP)==F)
tacsatSubset$SI_STATE <- na.omit(res)

#Convert character states into numeric
tacsatSubset$SI_STATE[tacsatSubset$SI_STATE=="h"] <- 1
tacsatSubset$SI_STATE[tacsatSubset$SI_STATE=="f"] <- 2
tacsatSubset$SI_STATE[tacsatSubset$SI_STATE=="s"] <- 3

#Lets look at the results
spds <- table(tacsatSubset$SI_SP,tacsatSubset$SI_STATE)
plot(y=spds[,1],x=an(dimnames(spds)[[1]]),type="h",ylim=range(spds),xlab="knots",ylab="Frequency")
lines(y=spds[,2],x=(an(dimnames(spds)[[1]])+0.1),col=2,type="h")
lines(y=spds[,3],x=(an(dimnames(spds)[[1]])+0.2),col=3,type="h")
legend("topright",legend=c("No fishing","Fishing","Steaming"),col=1:3,box.lty=0,lwd=2,lty=1)
```

From both these examples you can tell already that it isn't easy to determine activity, and it will most certainly take you a considerable amount of time to do this properly. However, it is very important as the results of this analysis are most often used as input to follow-up analyses!

### Fine scale effort patterns ###
With the gear type available in tacsat, we can now compare the effort available in the VMS data, and compare that with the effort data in the logbooks. Therefore let's first calculate the effort according to the VMS data.

```

#Start fresh
data(tacsat); data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#Sort the data and calculate the effort
tacsatp <- sortTacsat(tacsatp) #Make sure it is ordered correctly
tacsatp <- effort(tacsatp, level = "trip", unit = "hours",fill.na=T)
eflalop <- effort(eflalo,level="trip",unit="hours",fill.na=T)
eflalop$LE_EFF_LOG <- an(eflalop$LE_EFF_LOG)

#Let's quickly introduce a plotting routine, you'll hear more about that later on in the course
require(RColorBrewer)

#First the tacsat effort
plotTools(subset(tacsatp,LE_GEAR == "OTM" & year(SI_DATIM) == 1801),level="gridcell",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=F,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm="LE_EFF_VMS"))

#Now the eflalo effort
x11()
plotTools(subset(eflalop,LE_GEAR == "OTM" & year(SI_DATIM) == 1801),level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=F,color=NULL,control.eflalo=list(clm="LE_EFF_LOG"))
```

Quite a difference, not! Having VMS and logbook data available allows you to dispatch catches, values or effort at a much finer scale. Whether this is as good as, and as believable as normal logbook data is up to you to decide; the tools are available however.

### Distribution of catches / values ###
Similar to the analyses above, we can make a picture of, for example, plaice landings on the spatial scale of VMS. Note here, however, that we are about to dispatch landings at a much finer spatial scale than at which the original data are stored (ie. the logbooks). We might hope that this is a sensible interpretation of the two data sources available, exploiting each to permit the finest possible spatial resolution. It is up to you, however, to decide if this is a wise analysis to perform!
Make sure that the column 'LE\_CDAT' in the eflalo dataset contains real dates. Also, make sure you have at least 2 columns starting with either LE\_KG_or LE\_EURO_.

```

#This code will dispatch KG and EURO on the ping scale to VMS.  We need, however, to add activity to tacsatp (weve done this before, but let's take a quick and dirty approach here)
idx <- which(tacsatp$SI_SP >=2 & tacsatp$SI_SP <= 8)
tacsatp$SI_STATE <- 0
tacsatp$SI_STATE[idx] <- 1

#Crucial in this example below is that the elements 'kg in eflalo', 'kg in merged tacsat' and 'kg removed from eflalo' present you with identical values
# you might encounter some warnings or even errors, that is not a problem.

tacsatEflalo      <- splitAmongPings(tacsat=subset(tacsatp,LE_GEAR == "TBB" & year(SI_DATIM)==1801),eflalo=subset(eflalo,LE_GEAR == "TBB"),
variable="all",level="day",conserve=T)

#Quickly add a 'year' column to eflalo
eflalo$LE_CYEAR   <- year(as.POSIXct(eflalo$LE_CDAT, tz = "GMT", format = "%d/%m/%Y"))

#If 'conserve' was turned on in the function above, you would expect that tacsatEflalo would hold the same amount of KG and EURO as in the eflalo dataset
sum(tacsatEflalo[,kgeur(colnames(tacsatEflalo))])
sum(subset(eflalo,LE_GEAR == "TBB" & eflalo$LE_CYEAR==1801)[,kgeur(colnames(eflalo))])

#Lets make some plots again
plotTools(tacsatEflalo,level="gridcell",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=F,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm="LE_KG_PLE"))
eflalo$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,eflalo$FT_DTIME),format="%d/%m/%Y %H:%M")

x11()
plotTools(subset(eflalo,LE_GEAR == "TBB" & year(FT_DDATIM) == 1801),level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=F,color=NULL,control.eflalo=list(clm="LE_KG_PLE"))
```

In this practical you've come across a lot of example code where each of the functions also has many options. Once you get into these analyses yourselves back home, it won't hurt to read the manuals that go with these functions.

### Exercise 2 ###
  1. In the last example, how many kilo difference is there between all plaice landings in the tacsatEflalo set and in the eflalo dataset?
  1. Where might this difference come from?