# Practical 1: Getting the VMS and logbook data into R #

## Introduction ##
Welcome to the first practical of the ICES VMStools course. Within this series of practicals you will learn how to deal with VMS and logbook data to produce a variety of end results, ranging from the calculation of indicators of seabed disturbance to linking your data with, for example, sea temperature datasets.

The practicals you will work on will consist of some guidance text, some code that can be copied directly into R, and some examples to test your skills. In all instances, the R code necessary to do the job will be included and hence allows you, after the course, to immediately get started with your own data!

### It is recommended to carefully read through the text that goes with each of the practicals ###

## Getting VMStools installed ##
Before doing anything, make sure you have installed a good editor where you can 'store' your code. Examples are Tinn-R or R-studio. However, programs such as notepad work as well (but do not have function and bracket color recognition). It is also important that you have a good R version running. Make sure that it is [R2](https://code.google.com/p/vmstools/source/detail?r=2).14.x or higher (But not [R3](https://code.google.com/p/vmstools/source/detail?r=3).0.x yet)! (How can I check...)
```

R.Version()$version.string```
If this says anything below [R2](https://code.google.com/p/vmstools/source/detail?r=2).14.x, please upgrade. If this says [R3](https://code.google.com/p/vmstools/source/detail?r=3).0.x, you will run into problems too. Please use up to [R2](https://code.google.com/p/vmstools/source/detail?r=2).15.3 for this course. VMStools is still in development and R is as well. Over the past months considerable changes have taken place which do no not allow backward compatibility of VMStools with older R versions. VMStools will become available for [R3](https://code.google.com/p/vmstools/source/detail?r=3).0.x soon.

Now start a version of your favorite editor and start R. VMStools is an add-on package to R, and depends on several other packages. Hence, before we can run VMStools, we need to install these other packages too. This is how you do that:

**Note that all code inside a grey box starting with '>' refers to output of R. If it has no '>' in front, it refers to code of R that can directly be copied into R (but preference is given to go through the code line by line, to look at intermediate results etc).**

```

vmstoolsPackages <-
c("cluster","data.table","doBy","maps","mapdata","maptools","PBSmapping","sp")

for(i in vmstoolsPackages) try(install.packages(pkgs=i,repos=getOption("repos")))
```

A pop-up box will appear allowing you to select a country, and place, to download the packages from; choose one nearby (e.g. Denmark).

Now it is time to install the VMStools package, you have to download it first however from this website: [VMStools google repository](http://code.google.com/p/vmstools/downloads/list). Take either the vmstools\_0.66.zip version if you work under Windows or Mac and the vmstools\_0.66.tar.gz version if you work under Linux. Download the file to your local disk. Next you must install the VMStools package itself. Do this manually by clicking in R the 'Packages' tab and select the option 'Install packages from local zip files' and navigate to where the file was downloaded.  Alternatively it can be done, as follows, from the R command line:

```

path2VMStoolsDownload <- "C:/downloads/" #this folder is an example folder
install.packages(paste(path2VMStoolsDownload,"vmstools_0.66.zip",sep=""))
```

If all that went well, you have now successfully installed VMStools. Note that you only need to do this this once, and then all the packages will be available to you every time you start R.

## Getting started ##
Once more, open your editor and R version (if it is still running, that is fine as well). To load the VMStools package, you need to type the following at the 'R command lin'e or in the 'R interpreter':

```

library(vmstools)
```

All functions available in VMStools are loaded now and available to you. Note that a lot of messages are printed to R. You will need to make sure they do not contain errors, although a few warning messages are usually fine.
Getting your own data into R might be difficult right from the start, hence, we've included some example VMS and logbook data within the package itself.  These have been checked and cleaned allowing you to experiment with the package right away. You do, however, need to load these datasets first with the simple command below.

```

data(tacsat); data(eflalo)

#Let's have a look at what they look like
str(tacsat)

>'data.frame':   97015 obs. of  8 variables:
$ VE_COU : chr  "Atlantis" "Atlantis" "Atlantis" "Atlantis" ...
$ VE_REF : chr  "10" "10" "10" "10" ...
$ SI_LATI: num  51.4 51.4 51.4 51.4 51.4 ...
$ SI_LONG: num  3.58 3.58 3.6 3.59 3.59 ...
$ SI_DATE: chr  "04/05/1800" "04/05/1800" "04/05/1800" "04/05/1800" ...
$ SI_TIME: chr  "13:32:00" "15:28:00" "17:22:00" "19:18:00" ...
$ SI_SP  : num  0 0 0 0 0 0 0 12 2 2 ...
$ SI_HE  : num  25 25 25 25 25 25 111 292 43 294 ...

#- We can get a summary too
summarizeTacsat(tacsat)
```
Apparently, tacsat has class 'data.frame', has 97015 observations (rows), and 8 variables (columns). Each column name is given on a new row in the example above. The names might be a bit strange, but follow the agreed format as you can download from [Exchange\_EFLALO2\_v2-1.doc](http://code.google.com/p/vmstools/downloads/list).

### Exercise 1 ###
  1. Go through the columns of the dataset and try to understand what they mean.
  1. What are the different classes the tacsat data contains, and what do they mean?
  1. Which other common classes have you come across in R?

Have a quick look at the first 5 lines of the tacsat (VMS) data using the R function ```
head```.
```

head(tacsat)

>         VE_COU VE_REF  SI_LATI  SI_LONG    SI_DATE  SI_TIME SI_SP SI_HE
658237 Atlantis     10 51.44772 3.583731 04/05/1800 13:32:00     0    25
658238 Atlantis     10 51.44067 3.583847 04/05/1800 15:28:00     0    25
662888 Atlantis     10 51.44074 3.595529 04/05/1800 17:22:00     0    25
662889 Atlantis     10 51.44315 3.586862 04/05/1800 19:18:00     0    25
662890 Atlantis     10 51.44351 3.586948 04/05/1800 21:12:00     0    25
662891 Atlantis     10 51.44811 3.593958 04/05/1800 23:08:00     0    25
```
The vessel country is Atlantis, and the first vessel listed has 'name' 10. It's speed is zero and heading 25 degrees. It must be an old vessel as it was in this position in the 1800's!! Obviously, VMS wasn't invented by that time and this dataset is for example purposes only. However, it is real VMS data from a Dutch fishery, but we cannot give you more details because of confidentiality issues.

Let's do the same thing for the eflalo (logbook) data.

```

str(eflalo)
str(eflalo[,1:20])

>'data.frame':   4539 obs. of  20 variables:
$ VE_REF  : chr  "238" "731" "742" "742" ...
$ VE_FLT  : chr  "TBB" "TBB" "TBB" "TBB" ...
$ VE_COU  : chr  "Atlantis" "Atlantis" "Atlantis" "Atlantis" ...
$ VE_LEN  : num  24 23 22.1 22.1 23 ...
$ VE_KW   : num  221 221 221 221 221 ...
$ VE_TON  : num  NA NA NA NA NA NA NA NA NA NA ...
$ FT_REF  : chr  "271253" "271268" "271270" "271270" ...
$ FT_DCOU : chr  "nld" "nld" "nld" "nld" ...
$ FT_DHAR : chr  "SCH" "SL" "SCH" "SCH" ...
$ FT_DDAT : chr  "01/05/1800" "01/05/1800" "01/05/1800" "01/05/1800" ...
$ FT_DTIME: chr  "08:00:00" "15:00:00" "04:00:00" "04:00:00" ...
$ FT_LCOU : chr  "nld" "nld" "nld" "nld" ...
$ FT_LHAR : chr  "SL" "SL" "SL" "SL" ...
$ FT_LDAT : chr  "01/05/1800" "02/05/1800" "02/05/1800" "02/05/1800" ...
$ FT_LTIME: chr  "23:00:00" "07:00:00" "04:00:00" "04:00:00" ...
$ LE_ID   : chr  "271253-TBB-32F3" "271268-TBB-33F2" "271270-TBB-32F3" "271270-TBB-33F4" ...
$ LE_CDAT : chr  "01/05/1800" "01/05/1800" "01/05/1800" "01/05/1800" ...
$ LE_STIME: chr  NA NA NA NA ...
$ LE_ETIME: chr  NA NA NA NA ...
$ LE_SLAT : num  NA NA NA NA NA NA NA NA NA NA ...

#- Get a summary of eflalo
summarizeEflalo(eflalo)
```

This dataset contains a staggering amount of columns, 189! and only 4539 rows (observations). The data are structured in such a way that each species caught get's its own column. Hence, it might be more difficult to get the right information on your screen by just typing

```

head(eflalo)```

And therefore you might want to switch to

```

head(eflalo[,1:20])```

Again, there is a lot of information available. It starts with vessel 238, which apparently is a beamtrawler (VE\_FLT = TBB), has a vessel length of 24 meters and a 221kW engine. From column 32 onwards, the weight (kgs) of fish caught and their associated cash value is given (LE\_KG\_spec and LE\_EURO\_spec).

```

head(eflalo[,32:52])
```

### Exercise 2 ###
  1. Go through the different columns and try to understand what they mean.
  1. Two other classes are given in the eflalo dataset, what are they and do they make sense?
  1. How many kilos of plaice (PLE) has been registered in this dataset?

If you brought your own data, we can now try to get this into R, if not, we will just save the built-in datasets to disk and read them in.

#### Your own data ####
```

path2OwnData <- "C:/" #Fill out your own path here
dataName <- "tacsat2011.csv" #The tacsat file name, needs to be in csv format!
tacsat <- readTacsat(paste(path2OwnData,dataName,sep=""))
```

#### Example data ####
```

dir.create("C:/tmpTacsat")
write.table(tacsat, file = "C:/tmpTacsat/tacsat.csv", quote = TRUE, sep = ",",
eol = "\n", na = "NA", dec = ".", row.names = F,col.names = TRUE)

#Read in tacsat file
tacsat <- readTacsat("C:/tmpTacsat/tacsat.csv")
str(tacsat) #How does that look?

#Or we can use the build-in read.csv function from R
tacsat <- read.csv("C:/tmpTacsat/tacsat.csv")
str(tacsat) #How does that look?

tacsat <- formatTacsat(tacsat)
str(tacsat) #How about now?

```

Obviously the same thing works for eflalo. However, in most instances, the format of the data (i.e. the classes of each of the columns) is not correct. This will cause problems along the way, and therefore they must be put into the correct format using the VMStoosl R functions as follows:

```

tacsat <- formatTacsat(tacsat)
eflalo <- formatEflalo(eflalo)
```

### Exercise 3 ###
  1. Read in your own tacsat and eflalo data (or use the example data), check the formats, are they all correct?
  1. What happens if you run the formatTacsat and formatEflalo functions?
  1. Do the new formats make more sense?

## Lots of data ##
VMS and logbook data can be very large in size (e.g. multiple years, lot's of vessels etc) and reading them all in at once might pose problems. It isn't a problem to do the whole job on, for example, a year-by-year basis using a 'loop' and combine the datasets.

```

#First split the tacsat data into a 1800 and 1801 dataset
# To get a simple tacsat date timestap, use the sortTacsat function
tacsat <- sortTacsat(tacsat)
tacsat$YEAR <- year(tacsat$SI_DATIM)

#Split the data
tacsat1800 <- subset(tacsat,YEAR == 1800)
tacsat1801 <- subset(tacsat,YEAR == 1801)

#Check if both datasets together hold the same number of records as the original tacsat dataset
nrow(tacsat1800) + nrow(tacsat1801) #96882
nrow(tacsat) #97051!!! So apparently, some date + time's are not converted properly, up to you to fix this!

#Write files to disk
write.table(tacsat1800, file = "C:/tmpTacsat/tacsat1800.csv", quote = TRUE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = TRUE)

write.table(tacsat1801, file = "C:/tmpTacsat/tacsat1801.csv", quote = TRUE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = TRUE)

#Read it in again
years <- 1800:1801
path  <- "C:/tmpTacsat/"
totalTacsat <- numeric()
for(iYears in years){
tacsatYr <- readTacsat(paste(path,"tacsat",iYears,".csv",sep=""))
totalTacsat <- rbindTacsat(totalTacsat,tacsatYr)
}
head(totalTacsat) #Note the first rows of NA's, needs to be removed
totalTacsat <- totalTacsat[-1,]
nrow(totalTacsat) #96882 again!
```

### Exercise 4 ###
  1. Split the example dataset into a set with vessel names smaller than 1000 and bigger than thousand and save these to disk.
  1. Read the separate files in and combine them into one 'totalTacsat' file.
  1. How many rows of data do you have, and how many columns? Does that match up with the number of rows and columns mentioned earlier in this practical?