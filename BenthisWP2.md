# Benthis WP2 workflow example #

Welcome to the Benthis workflow example wiki. The page is hosted on the VMStools page
as for this analyses you will be required to install the VMStools software for R.
Please note that at this stage (May 2013), VMStools will only operate on R-2.15.x versions.

## Objectives ##
  * To map habitat types and sea bed impact from fishing activities in EU waters to identify where fisheries potentially compromises seafloor integrity and conflicts of ecosystem services
  * To develop and implement new methodology, combining VMS, logbook and industry data, to assess actual seabed impact from large scale fishing activities on an appropriate spatial and temporal scale

## workflow ##

Please find below an example workflow to treat VMS and Logbook in a systematic manner to create pressure maps (swept area and severity). The workflow should ultimately produce a dataset with VMS positions coupled to an estimate of swept area per VMS ping for each country, further combined across countries. Swept areas are further broken-down to each of the sub-components of the gear and these sub-components are assigned a specific severity of impact.

### Installing VMStools ###
  * The VMStools software is available [here](http://code.google.com/p/vmstools/downloads/list)
  * The procedure for installing the VMStools library under R is described [here](http://code.google.com/p/vmstools/wiki/Introduction?tm=6)

### Loading the required libraries into the R console ###
We load three different libraries, to start with the VMStools one, followed by two packages that help us to create maps. If you have not installed them, please use `install.packages("maps",repos=getOption("repos"))` and select the correct packages.
```

rm(list=ls())
library(vmstools)
library(maps)
library(mapdata)
memory.size(4000)
```

In general, workflows can be easily adapted to each users' computer by setting working directories first. You may need to adjust these to direct to the <br> folders you would like to store the results in.<br>
<br>
<h3>setting working directories</h3>

<pre><code><br>
codePath  &lt;- "D:/Benthis/Code/R/"<br>
dataPath  &lt;- "D:/Benthis/Data/"<br>
pricePath &lt;- "D:/Benthis/Data/Price/"<br>
outPath   &lt;- "D:/Benthis/Output/"<br>
polPath   &lt;- "D:/Benthis/Data/Polygons/"<br>
</code></pre>

<h3>loading VMS and Logbook data</h3>

The VMS data we use in this analyses must be formatted to the 'tacsat' style. The logbook data must be formatted in the 'eflalo' style. You can download a description of these formats <a href='http://code.google.com/p/vmstools/downloads/list'>here</a> and if you want another example, in the VMStools package an example tacsat and eflalo dataset are embedded. Load the VMStools library and call these by typing <code>data(tacsat)</code> or <code>data(eflalo)</code>.<br>
If you have the data in R format already, you can read them in as below.<br>
Make sure that they are in the right column format too!<br>
<pre><code><br>
load(file.path(dataPath,"tacsat.RData")); # get the tacsat object<br>
load(file.path(dataPath,"eflalo.RData")); # get the eflalo object<br>
tacsat &lt;- formatTacsat(tacsat) # format each of the columns to the specified class<br>
eflalo &lt;- formatEflalo(eflalo) # format each of the columns to the specified class<br>
</code></pre>

If they are in a .csv format, you could use the following lines.<br>
<br>
<pre><code><br>
tacsat &lt;- readTacsat(file.path(dataPath,"tacsat.csv"))<br>
eflalo &lt;- readEflalo(file.path(dataPath,"eflalo.csv"))<br>
</code></pre>

Make sure that the fishing activities / metier (the 'LE_MET' column in eflalo) is present and complete.<br>
<br>
<h3>Country specific variables</h3>

<pre><code><br>
VMS_ping_rate_in_hour &lt;- 2<br>
</code></pre>

<h3>load VMStools additional datasets</h3>
<pre><code><br>
data(euharbours)<br>
data(ICESareas)<br>
data(europa)<br>
</code></pre>

Let's take a look at e.g. the ICES areas and harbours dataset.<br>
<br>
<pre><code><br>
map("worldHires",xlim=c(-15,15),ylim=c(40,65),fill=T,col="#7FBC41",oma=rep(0,4),mar=rep(2.5,4))<br>
plot(ICESareas,add=T,col="lightblue")<br>
points(harbours$lon,harbours$lat,pch=16,col="blue",cex=0.4)<br>
map.axes()<br>
</code></pre>
It should look something like this:<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/harbours.png' />


<h3>Cleaning the tacsat dataset</h3>
In general, the tacsat dataset contains many errors, such as duplicates, points that are on land or not even on the globe. These need to be filtered out before we execute the analyses. Most of these 'problems' have been identified by different authors of peer-reviewed papers, and we treat the most common ones below.<br>
<br>
<br>
<ul><li>It is a good idea to keep track of those records you have removed<br>
<pre><code><br>
remrecsTacsat     &lt;- matrix(NA,nrow=6,ncol=2,dimnames= list(c("total","duplicates","notPossible",<br>
"pseudoDuplicates","harbour","land"),<br>
c("rows","percentage")))<br>
remrecsTacsat["total",] &lt;- c(nrow(tacsat),"100%")<br>
</code></pre></li></ul>

<ul><li>Remove duplicate records<br>
<pre><code><br>
tacsat$SI_DATIM &lt;- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "),<br>
tz="GMT", format="%d/%m/%Y  %H:%M")<br>
uniqueTacsat    &lt;- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,tacsat$SI_DATIM)<br>
tacsat          &lt;- tacsat[!duplicated(uniqueTacsat),]<br>
remrecsTacsat["duplicates",] &lt;- c(nrow(tacsat),100+round((nrow(tacsat) -<br>
an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove points that cannot be possible<br>
<pre><code><br>
spThres         &lt;- 20   #Maximum speed threshold in analyses in nm<br>
idx           &lt;- which(abs(tacsat$SI_LATI) &gt; 90 | abs(tacsat$SI_LONG) &gt; 180)<br>
idx           &lt;- unique(c(idx,which(tacsat$SI_HE &lt; 0 | tacsat$SI_HE &gt; 360)))<br>
idx           &lt;- unique(c(idx,which(tacsat$SI_SP &gt; spThres)))<br>
if(length(idx)&gt;0) tacsat          &lt;- tacsat[-idx,]<br>
remrecsTacsat["notPossible",] &lt;- c(nrow(tacsat),100+round((nrow(tacsat) -<br>
an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove points which are pseudo duplicates as they have an interval rate < x minutes<br>
<pre><code><br>
intThres        &lt;- 5    # Minimum difference in time interval in minutes to prevent pseudo duplicates<br>
tacsat          &lt;- sortTacsat(tacsat)<br>
tacsatp         &lt;- intervalTacsat(tacsat,level="vessel",fill.na=T)<br>
tacsat          &lt;- tacsatp[which(tacsatp$INTV &gt; intThres | is.na(tacsatp$INTV)==T),-grep("INTV",colnames(tacsatp))]<br>
remrecsTacsat["pseudoDuplicates",] &lt;- c(nrow(tacsat),100+round((nrow(tacsat) -<br>
an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove points in harbour<br>
<pre><code><br>
idx             &lt;- pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,harbours)<br>
pih             &lt;- tacsat[which(idx == 1),]<br>
save(pih,file=paste(outPath,"pointInHarbour.RData",sep=""))<br>
tacsat          &lt;- tacsat[which(idx == 0),]<br>
remrecsTacsat["harbour",] &lt;- c(nrow(tacsat),100+round((nrow(tacsat) -<br>
an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove points on land<br>
<pre><code><br>
pols            &lt;- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),<br>
function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,SI_LATI=subset(europa,SID==x)$Y)}))<br>
idx             &lt;- pointOnLand(tacsat,pols);<br>
pol             &lt;- tacsat[which(idx == 1),]<br>
save(pol,file=paste(outPath,"pointOnLand.RData",sep=""))<br>
tacsat          &lt;- tacsat[which(idx == 0),]<br>
remrecsTacsat["land",] &lt;- c(nrow(tacsat),100+round((nrow(tacsat) -<br>
an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))<br>
</code></pre>
If you would plot those points on land, it would look something like this:<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/pointsOnLand.png' /></li></ul>


<ul><li>Save the remrecsTacsat file<br>
<pre><code><br>
save(remrecsTacsat,file=file.path(outPath,"remrecsTacsat.RData"))<br>
</code></pre>
The result, based on the VMStools tacsat dataset, looks like:<br>
<table><thead><th> <b>Category</b> </th><th> <b>Records left</b> </th><th> <b>Percentage of original</b> </th></thead><tbody>
<tr><td> total           </td><td> 97015               </td><td> 100%                          </td></tr>
<tr><td> duplicates      </td><td> 97015               </td><td> 100%                          </td></tr>
<tr><td> not possible    </td><td> 97001               </td><td> 99.99%                        </td></tr>
<tr><td> pseudo duplicates </td><td> 84774               </td><td> 87.38%                        </td></tr>
<tr><td> in harbour      </td><td> 70948               </td><td> 73.13%                        </td></tr>
<tr><td> on land         </td><td> 70063               </td><td> 72.22%                        </td></tr></li></ul></tbody></table>

<ul><li>Save the cleaned tacsat file<br>
<pre><code><br>
save(tacsat,file=file.path(outPath,"cleanTacsat.RData"))<br>
</code></pre></li></ul>


<h3>cleaning eflalo</h3>
Very similar to cleaning the tacsat dataset, we need to clean the eflalo dataset too.<br>
<br>
<ul><li>Keep track of removed points<br>
<pre><code><br>
remrecsEflalo     &lt;- matrix(NA,nrow=5,ncol=2,dimnames=list(c("total","duplicated","impossible time",<br>
"before 1st Jan","departArrival"),<br>
c("rows","percentage")))<br>
remrecsEflalo["total",] &lt;- c(nrow(eflalo),"100%")<br>
</code></pre></li></ul>

<ul><li>Remove non-unique trip numbers<br>
<pre><code><br>
eflalo            &lt;- eflalo[!duplicated(paste(eflalo$LE_ID,eflalo$LE_CDAT,sep="-")),]<br>
remrecsEflalo["duplicated",] &lt;- c(nrow(eflalo),100+round((nrow(eflalo) -<br>
an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove impossible time stamp records<br>
<pre><code><br>
eflalo$FT_DDATIM  &lt;- as.POSIXct(paste(eflalo$FT_DDAT,eflalo$FT_DTIME, sep = " "),<br>
tz = "GMT", format = "%d/%m/%Y  %H:%M")<br>
eflalo$FT_LDATIM  &lt;- as.POSIXct(paste(eflalo$FT_LDAT,eflalo$FT_LTIME, sep = " "),<br>
tz = "GMT", format = "%d/%m/%Y  %H:%M")<br>
<br>
eflalo            &lt;- eflalo[!(is.na(eflalo$FT_DDATIM) |is.na(eflalo$FT_LDATIM)),]<br>
remrecsEflalo["impossible time",] &lt;- c(nrow(eflalo),100+round((nrow(eflalo) -<br>
an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove trip starting before 1st Jan<br>
<pre><code><br>
year              &lt;- min(year(eflalo$FT_DDATIM))<br>
eflalo            &lt;- eflalo[eflalo$FT_DDATIM&gt;=strptime(paste(year,"-01-01 00:00:00",sep=''),<br>
"%Y-%m-%d %H:%M:%S"),]<br>
remrecsEflalo["before 1st Jan",] &lt;- c(nrow(eflalo),100+round((nrow(eflalo) -<br>
an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Remove records with arrival date before departure date<br>
<pre><code><br>
eflalop           &lt;- eflalo<br>
eflalop$FT_DDATIM &lt;- as.POSIXct(paste(eflalo$FT_DDAT,  eflalo$FT_DTIME,   sep=" "),<br>
tz="GMT", format="%d/%m/%Y  %H:%M")<br>
eflalop$FT_LDATIM &lt;- as.POSIXct(paste(eflalo$FT_LDAT,  eflalo$FT_LTIME,   sep=" "),<br>
tz="GMT", format="%d/%m/%Y  %H:%M")<br>
idx               &lt;- which(eflalop$FT_LDATIM &gt;= eflalop$FT_DDATIM)<br>
eflalo            &lt;- eflalo[idx,]<br>
remrecsEflalo["departArrival",] &lt;- c(nrow(eflalo),100+round((nrow(eflalo) -<br>
an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))<br>
</code></pre></li></ul>

<ul><li>Save the remrecsEflalo file<br>
<pre><code><br>
save(remrecsEflalo,file=file.path(outPath,"remrecsEflalo.RData"))<br>
</code></pre></li></ul>

The result, based on the VMStools eflalo dataset, looks like:<br>
<table><thead><th> <b>Category</b> </th><th> <b>Records left</b> </th><th> <b>Percentage of original</b> </th></thead><tbody>
<tr><td> total           </td><td> 4539                </td><td> 100%                          </td></tr>
<tr><td> duplicated      </td><td> 4446                </td><td> 97.95%                        </td></tr>
<tr><td> impossible time </td><td> 4446                </td><td> 97.95%                        </td></tr>
<tr><td> before 1st of Jan </td><td> 4446                </td><td> 97.95%                        </td></tr>
<tr><td> departure before arrival </td><td> 4446                </td><td> 97.95%                        </td></tr></tbody></table>

<ul><li>Save the cleaned eflalo file<br>
<pre><code><br>
save(eflalo,file=file.path(outPath,"cleanEflalo.RData"))<br>
</code></pre></li></ul>


As we now have cleaned both the tacsat and eflalo datasets, we can start to explore<br>
the potential of combining these datasets. Both datasets have two aspects in common:<br>
vessel name or identifier and a date / time stamp. On this basis, both datasets<br>
can be merged together. One tries to find those tacsat records that fall within<br>
the trip departure and arrival date registered in the logbooks.<br>
<br>
In VMStools, a special function is designed to do this automatically.<br>
<br>
<h3>merge eflalo and tacsat</h3>

The code itself is simple. Note that the trip identifier (FT_REF) in the eflalo<br>
dataset is assigned to the tacsat file if a match is found. If not, the FT_REF<br>
column in the tacsat dataset will state '0'.<br>
<pre><code><br>
tacsatp           &lt;- mergeEflalo2Tacsat(eflalo,tacsat)<br>
</code></pre>

As we've now established some links, we can assign more attributes from the<br>
eflalo data to the tacsat data, such as gear characteristics (LE_GEAR), vessel<br>
length (VE_LEN) or fishery type (LE_MET).<br>
<br>
<ul><li>Link gear characteristics to tacsat</li></ul>

<pre><code><br>
tacsatp$LE_GEAR   &lt;- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]<br>
tacsatp$VE_LEN    &lt;- eflalo$VE_LEN[ match(tacsatp$FT_REF,eflalo$FT_REF)]<br>
tacsatp$LE_MET    &lt;- eflalo$LE_MET[ match(tacsatp$FT_REF,eflalo$FT_REF)]<br>
save(tacsatp,   file=paste(outPath,"tacsatMerged.RData",   sep=""))<br>
</code></pre>
In some cases, no match could be found. The result will look something like this:<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/missingMatch.png' />

<ul><li>read in look-up table (industry data from interviews / experimental studies)<br>
To calculate bottom impact in a later stage, we need the width of the gear. It is not<br>
always easy to derive these widths though.<br>
The gear width for a (beam) trawl is straightforward, the gear width for a seine<br>
could be the length of the net and the gear width of a gillnets or other could be<br>
the number of nets released per operation times the unit net length. However<br>
, this needs to be confirmed by the industry and potentially making use of a functional<br>
relationships between vessel HP or length to gear size.</li></ul>

A table specifying all these gear - width combinations should be put together and<br>
read into the workflow.<br>
<br>
This is making use of a table of parameters specifying the relationships linking gear widths (door spread, beam width, dredge width, or length of the seine rope) to vessel specifications (KW or LOA).<br>
<br>
Please remember that before making this step you first need to link  your metier definition in tacsatp$LE_MET to the BENTHIS metier classification included in this table. This can be done by renaming by using the very useful 'levels()<-' R function.<br>
<br>
<pre><code><br>
# a fake example for illustration<br>
dd&lt;- data.frame(LE_MET=c("OTB_DEM_&gt;120mm", "OTM_SPF_&gt;16mm"), data=c(1,1))<br>
print(levels(dd$LE_MET))<br>
levels(dd$LE_MET) &lt;- c('OT_DMF', 'OT_SPF') # renaming the levels<br>
print(dd)<br>
#=&gt; in your true case, replace 'dd' by 'tacsatp' and adapt for renaming your list of metiers.<br>
</code></pre>

The table of parameters is the outcome of the <a href='https://code.google.com/p/vmstools/wiki/Industry_data'>questionnaire survey analysis</a> and will be circulated by mail to partners.<br>
<br>
Only the length of the rope is predicted for seiners at this stage of the work flow. The swept area for the seiners (computed from the length of the rope) is addressed further below in the workflow.<br>
<br>
<pre><code><br>
<br>
# import the param table obtained from the industry_data R analyses<br>
gear_param_per_metier       &lt;- read.table(file=file.path(outPath, "estimates_for_gear_param_per_metier.txt"))<br>
<br>
GearWidth                   &lt;- tacsatp[!duplicated(data.frame(tacsatp$VE_REF,tacsatp$LE_MET)), ]<br>
GearWidth                   &lt;- GearWidth[,c('VE_REF','LE_MET','VE_KW', 'VE_LEN') ]<br>
GearWidth$GEAR_WIDTH        &lt;- NA<br>
GearWidth$GEAR_WIDTH_LOWER  &lt;- NA<br>
GearWidth$GEAR_WIDTH_UPPER  &lt;- NA<br>
for (i in 1:nrow(GearWidth)) { # brute force...<br>
<br>
kW      &lt;- GearWidth$VE_KW[i]<br>
<br>
LOA     &lt;- GearWidth$VE_LEN[i]<br>
<br>
this    &lt;- gear_param_per_metier[gear_param_per_metier$a_metier==GearWidth$LE_MET[i],]<br>
a &lt;- NULL ; b &lt;- NULL<br>
a       &lt;- this[this$param=='a', 'Estimate']<br>
b       &lt;- this[this$param=='b', 'Estimate']<br>
GearWidth[i,"GEAR_WIDTH"]  &lt;-   eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km<br>
a       &lt;- this[this$param=='a', 'Estimate'] +2*this[this$param=='a', 'Std..Error']<br>
b       &lt;- this[this$param=='b', 'Estimate'] +2*this[this$param=='b', 'Std..Error']<br>
GearWidth[i,"GEAR_WIDTH_UPPER"]  &lt;-  eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km<br>
a       &lt;- this[this$param=='a', 'Estimate'] -2*this[this$param=='a', 'Std..Error']<br>
b       &lt;- this[this$param=='b', 'Estimate'] -2*this[this$param=='b', 'Std..Error']<br>
GearWidth[i,"GEAR_WIDTH_LOWER"]  &lt;-  eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km<br>
}<br>
save(GearWidth, file=file.path(dataPath,"gearWidth.RData"))<br>
load(file.path(dataPath, "gearWidth.RData"))<br>
tacsatp                     &lt;- merge(tacsatp, GearWidth,by=c("VE_REF","LE_MET"),<br>
all.x=T,all.y=F)<br>
save(tacsatp,   file=file.path(outPath,"tacsatMergedWidth.RData"))<br>
<br>
</code></pre>

<h3>Separate tacsat that is not merged</h3>
<pre><code><br>
tacsatpmin        &lt;- subset(tacsatp,FT_REF == 0)<br>
save(tacsatpmin, file=paste(outPath,"tacsatNotMerged.RData",sep=""))<br>
</code></pre>

<h3>define (fishing) activity</h3>
A common way to define activity of a vessel is to look at its speed profiles.<br>
In general, you can identify for most bottom trawling gear 3 peaks. One peak<br>
would be situated right at the zero-knots boundary, one peak somewhere between<br>
4-8 knots and a third one between 8 and 12 knots. For example, see the speed profile<br>
below (top: instantaneous speeds, below: calculated speeds):<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/speedProfiles.png' />


The following step identifies the fishing activities from the steaming phases<br>
by analysing the vessel-specific speed profile and knowing the type of gear used from the logbooks.<br>
<br>
<br>
First remove points with NA's in them in critial places<br>
<pre><code><br>
idx               &lt;- which(is.na(tacsatp$VE_REF) == T   | is.na(tacsatp$SI_LONG) == T | is.na(tacsatp$SI_LATI) == T |<br>
is.na(tacsatp$SI_DATIM) == T |  is.na(tacsatp$SI_SP) == T)<br>
if(length(idx)&gt;0) tacsatp         &lt;- tacsatp[-idx,]<br>
</code></pre>

The identification of these activities is sensitive to starting conditions on where<br>
the peaks are approximately. Hence, a pre-analyses needs to take place to roughly<br>
estimate where each of the three peaks is situated. You can either spot 2 or 3 peaks<br>
(one around zero knots, around 4-8 knots and around 8-12 knots or one around zero<br>
nots and one around 6-12 knots). To reliably estimate these peaks, the speed profile<br>
is mirrored in the y-axis and therefore you can either spot 3 or 5 peaks. The<br>
following routine will ask for the number of peaks and where they roughly are<br>
situated.<br>
<pre><code><br>
storeScheme       &lt;- activityTacsatAnalyse(tacsatp, units = "year", analyse.by = "LE_GEAR",identify="means")<br>
storeScheme       &lt;- storeScheme[-which(is.na(storeScheme$analyse.by)==T),]<br>
save(storeScheme, file=paste(outPath,"storeScheme.RData",sep=""))<br>
</code></pre>
In case you have many vessels, you could create the storeScheme object yourself outside<br>
of this function.<br>
<br>
It is in many cases easier to perform the calculation on a year-by-year basis.<br>
<pre><code><br>
storeScheme       &lt;- storeScheme[storeScheme$years==1800,]<br>
storeScheme$years &lt;- as.numeric(as.character(storeScheme$years))<br>
tacsatp$year      &lt;- format(tacsatp$SI_DATIM, "%Y")<br>
tacsatp           &lt;- tacsatp[tacsatp$year %in% 1800,]<br>
</code></pre>

The following routine is the function that actually identifies the peaks and the<br>
result is returned in the 'activity' object.<br>
<pre><code><br>
activity          &lt;- activityTacsat(tacsatp,units="year",analyse.by="LE_GEAR",storeScheme,<br>
plot=TRUE,level="all",sigma=1)<br>
tacsatp$SI_STATE  &lt;- activity<br>
</code></pre>

In most cases however, you have still a few records that have no activity associated<br>
as yet. You need to 'fix' this manually. An example is given below.<br>
<br>
<pre><code><br>
idx               &lt;- which(is.na(tacsatp$SI_STATE))<br>
tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] &gt;= 1.5 &amp;<br>
tacsatp$SI_SP[idx] &lt;= 7.5)]] &lt;- 'f'<br>
tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] &lt;  1.5)]] &lt;- 'h'<br>
tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] &gt;  7.5)]] &lt;- 's'<br>
</code></pre>

Make sure you save the results, else you might need to do the analyse all over!<br>
<pre><code><br>
save(storeScheme, file=file.path(outPath,"storeScheme.RData"))<br>
save(tacsatp,     file=file.path(outPath,"tacsatActivity.RData"))<br>
</code></pre>

<h3>Alternative activity identification: (Based on Fock 2008)</h3>

We can read in a table that specifies speed ranges for vessels fishing for different<br>
gears. This table can either be constructed on the basis of published work by<br>
Fock 2008 (Fisheries in the context of marine spatial planning: Defining principal<br>
areas for fisheries in the German EEZ) or on a dataset created yourself.<br>
<br>
<pre><code><br>
gearSpeeds            &lt;- read.csv(file.path(dataPath,"gearSpeeds.csv"))<br>
colnames(gearSpeeds)  &lt;- c("LE_GEAR","SI_SPmin","SI_SPmax")<br>
#- In this example, gearSpeeds is a table with 3 columns. One column = gear code,<br>
#   second column is the minimum speed associatd with fishing and the third<br>
#   column is the maximum speed associated with fishing<br>
<br>
tacsatp$SI_STATE &lt;- 's'<br>
for(iGear in sort(gearSpeeds$LE_GEAR)){<br>
idxGearSpeed  &lt;- subset(gearSpeeds,LE_GEAR == iGear)<br>
idxtacsat     &lt;- which(tacsatp$LE_GEAR == iGear &amp;<br>
tacsatp$SI_SP &gt;= gearSpeeds$SI_SPmin &amp;<br>
tascatp$SI_SP &lt;= gearSpeeds$SI_SPmax)<br>
tacsatp$SI_STATE[idx] &lt;- 'f'<br>
}<br>
</code></pre>


<h3>Labelling each haul</h3>
Assign an identifier (HL_ID) to each of the fishing sequences. This can especially<br>
be useful when you want to count the number of hauls in an aggregated grid.<br>
<br>
<pre><code><br>
tacsatp   &lt;- labellingHauls(tacsatp)<br>
</code></pre>

<h3>Interpolate tacsat data (or not)</h3>
In some occasions, and especially for towed trawl gear, it can be useful to artificially<br>
increase the density of VMS points. Through the use of an interpolation algorithm<br>
this is possible. Based on speed and heading a spline will be fitted between two<br>
consecutive fishing pings and the result can be converted into the tacsat format<br>
again, now containing an additional set of intermediate pings.<br>
<br>
<pre><code><br>
load(file=paste(outPath,"tacsatActivity.RData",   sep=""))<br>
<br>
#- Select towed gears only<br>
towed_gears       &lt;- c('OTB', 'TBB', 'PTB', 'PTM', 'DRB')  # TO DO: list to be checked<br>
tacsatp           &lt;- orderBy(~VE_REF+SI_DATIM,data=tacsatp)<br>
</code></pre>

As these interpolations can become rather large, it is sensible to loop over the<br>
gears and store the results in between. Therefore, we first create an output folder<br>
by typing <code>dir.create(file.path(outPath, "interpolated"))</code>.<br>
<br>
Interpolation should only occur on the fishing points, so we need to get rid to<br>
all non-fishing points. To be conservative, we start each 'haul' with a steaming<br>
point and finish with a steaming point too.<br>
<br>
<pre><code><br>
tacsatp$SI_STATE_num &lt;- NA<br>
tacsatp$SI_STATE_num[which(tacsatp$SI_STATE=="h")] &lt;- 1<br>
tacsatp$SI_STATE_num[tacsatp$SI_STATE=="f"] &lt;- 2<br>
tacsatp$SI_STATE_num[tacsatp$SI_STATE=="s"] &lt;- 3<br>
is_transition     &lt;- c(0,diff(tacsatp$SI_STATE_num))<br>
is_transition2    &lt;- c(diff(tacsatp$SI_STATE_num), 0)<br>
tacsatp           &lt;- tacsatp[ !is.na(tacsatp$SI_STATE_num) &amp; (tacsatp$SI_STATE_num ==2 |<br>
is_transition!=0 | is_transition2!=0),]<br>
tacsatp           &lt;- tacsatp[,-grep("SI_STATE_num",colnames(tacsatp))]<br>
tacsatp$SI_STATE  &lt;- "f"<br>
</code></pre>

Now that we remain with fishing points only, we can perform the interpolation<br>
<br>
<pre><code><br>
for(iGr in towed_gears){<br>
tacsatpGear        &lt;- tacsatp[!is.na(tacsatp$LE_GEAR) &amp; tacsatp$LE_GEAR==iGr,]<br>
<br>
for(iVE_REF in sort(unique(tacsatpGear$VE_REF))){<br>
tacsatpGearVEREF &lt;- tacsatpGear[tacsatpGear$VE_REF %in% iVE_REF,]<br>
if(nrow(tacsatpGearVEREF)&gt;3) {<br>
<br>
#- Interpolate according to the cubic-hermite spline interpolation<br>
interpolationcHs &lt;- interpolateTacsat(tacsatpGearVEREF,<br>
interval=VMS_ping_rate_in_hour*60,<br>
margin=10, # i.e. will make disconnected interpolations if interval out of the 110 130 min range<br>
res=100,<br>
method="cHs",<br>
params=list(fm=0.2,distscale=20,sigline=0.2,st=c(2,6)),<br>
headingAdjustment=0,<br>
fast=FALSE)<br>
<br>
#- Convert the interpolation to tacsat style data<br>
tacsatIntGearVEREF &lt;- interpolation2Tacsat(interpolationcHs, tacsatpGearVEREF)<br>
<br>
#- Correct for HL_ID skipping (not in use if labellingHauls() called before)<br>
#tacsatIntGearVEREF$HL_ID[which(diff(tacsatIntGearVEREF$HL_ID)&lt;0)] &lt;- tacsatIntGearVEREF$HL_ID[(which(diff(tacsatIntGearVEREF$HL_ID)&lt;0)-1)]<br>
<br>
#  the swept area (note that could work oustide the loop area as well....)<br>
tacsatIntGearVEREF$SWEPT_AREA_KM2 &lt;- NA<br>
a_dist &lt;- distance(c(tacsatIntGearVEREF$SI_LONG[-1],0),  c(tacsatIntGearVEREF$SI_LATI[-1],0),<br>
tacsatIntGearVEREF$SI_LONG, tacsatIntGearVEREF$SI_LATI)<br>
tacsatIntGearVEREF$SWEPT_AREA_KM2 &lt;- a_dist * tacsatIntGearVEREF$GEAR_WIDTH<br>
tacsatIntGearVEREF$SWEPT_AREA_KM2_LOWER &lt;- a_dist * tacsatIntGearVEREF$GEAR_WIDTH_LOWER<br>
tacsatIntGearVEREF$SWEPT_AREA_KM2_UPPER &lt;- a_dist * tacsatIntGearVEREF$GEAR_WIDTH_UPPER<br>
<br>
<br>
<br>
<br>
save(tacsatIntGearVEREF, file=file.path(outPath, "interpolated",<br>
paste("tacsatInt",iVE_REF, "_", iGr, ".RData", sep="")))<br>
}<br>
}<br>
}<br>
</code></pre>

This is how an interpolation could look like, including the sampled points from it.<br>
The blue points are the original VMS pings while the black solid line represents<br>
the connection of all interpolated points (100). The 8 intermediate red points are<br>
sampled from this set of 100 (representing equal distance between points).<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/interpolation.png' />

Additionally, we can calculate the swept area as we know its distance travelled<br>
and the width of the gear.<br>
<br>
<pre><code><br>
tacsatIntGearVEREF$SWEPT_AREA_KM2 &lt;- NA<br>
tacsatIntGearVEREF$SWEPT_AREA_KM2 &lt;- distance(c(tacsatIntGearVEREF$SI_LONG[-1],0),  c(tacsatIntGearVEREF$SI_LATI[-1],0),<br>
tacsatIntGearVEREF$SI_LONG,           tacsatIntGearVEREF$SI_LATI) *<br>
tacsatIntGearVEREF$GEAR_WIDTH<br>
</code></pre>

<h3>Buffer area</h3>

For those vessels where we cannot interpolate (because it doesn't make sense,<br>
we can apply a certain buffer area, representing some sort of circle around a<br>
VMS ping. To derive the swept area from that circle we apply simple mathematics.<br>
<br>
This step is concerning Danish and Scottish seines where the swept area is deduced from the total length of the rope (stored in the (confusing for this case) LE_GEAR_WIDTH field).<br>
<br>
Note that all passive gears such as gillnets, pots, etc. are not considered in the workflow because considered irrelevant for bottom impact assessment in WP2.<br>
<br>
We do that here for those gear types that were not in our towed_gear list.<br>
<br>
<pre><code><br>
all_gears            &lt;- sort(unique(tacsatp$LE_GEAR))<br>
seine_gears          &lt;- c('SDN', 'SSC')  # Danish and Scottish seiners (all other gears should not be included because irrelevant for bottom impact).<br>
for(iGr in seine_gears){<br>
tacsatpGear        &lt;- tacsatp[!is.na(tacsatp$LE_GEAR) &amp; tacsatp$LE_GEAR==iGr,]<br>
<br>
for(iVE_REF in sort(unique(tacsatpGear$VE_REF)){<br>
tacsatpGearVEREF &lt;- tacsatpGear[tacsatpGear$VE_REF %in% iVE_REF,]<br>
tacsatpGearVEREF &lt;- tacsatpGearVEREF[tacsatpGearVEREF$SI_STATE=='f',] # keep fishing pings only<br>
<br>
tacsatpGearVEREF$SWEPT_AREA_KM2         &lt;- pi*(tacsatpGearVEREF$GEAR_WIDTH/(2*pi))^2<br>
tacsatpGearVEREF$SWEPT_AREA_KM2_LOWER   &lt;- pi*(tacsatpGearVEREF$GEAR_WIDTH_LOWER/(2*pi))^2<br>
tacsatpGearVEREF$SWEPT_AREA_KM2_UPPER   &lt;- pi*(tacsatpGearVEREF$GEAR_WIDTH_UPPER/(2*pi))^2<br>
<br>
haul_duration                           &lt;- 3 # assumption of a mean duration based from questionnaires to seiners<br>
tacsatpGearVEREF$SWEPT_AREA_KM2         &lt;- tacsatpGearVEREF$SWEPT_AREA_KM2 * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.<br>
tacsatpGearVEREF$SWEPT_AREA_KM2_LOWER         &lt;- tacsatpGearVEREF$SWEPT_AREA_KM2_LOWER * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.<br>
tacsatpGearVEREF$SWEPT_AREA_KM2_UPPER         &lt;- tacsatpGearVEREF$SWEPT_AREA_KM2_UPPER * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.<br>
idx                                     &lt;- grep('SSC', as.character(tacsatpGearVEREF$LE_GEAR))<br>
tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2'] &lt;- tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2'] *1.5 # ad hoc correction to account for the SSC specificities<br>
tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2_LOWER'] &lt;- tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2_LOWER'] *1.5 # ad hoc correction to account for the SSC specificities<br>
tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2_UPPER'] &lt;- tacsatpGearVEREF[idx, 'SWEPT_AREA_KM2_UPPER'] *1.5 # ad hoc correction to account for the SSC specificities<br>
<br>
tacsatIntGearVEREF &lt;- tacsatpGearVEREF<br>
<br>
<br>
save(tacsatIntGearVEREF, file=file.path(outPath, "interpolated",<br>
paste("tacsatSweptArea_",iVE_REF, "_", iGr, ".RData", sep="")))<br>
}<br>
}<br>
</code></pre>

<h3>Create dataset incl. swept area and effort</h3>

As we've now calculated the swept area for both the interpolated gear types and<br>
the remaining ones via the buffer method, we can now combine all of them.<br>
<br>
<pre><code><br>
allInts &lt;- list()<br>
counter &lt;- 0<br>
for(iGr in all_gears){<br>
for(iVE_REF in sort(unique(tacsatpVE_REF))){<br>
<br>
dat &lt;- try(get(load(file=file.path(outPath, "interpolated",<br>
paste("tacsatSweptArea_",iVE_REF, "_", iGr, ".RData", sep="")))))<br>
#- Assign effort to dataset<br>
tacsatp$effort_KWdays    &lt;- intervalTacsat(tacsatp,level="trip",fill.na=T)$INTV * as.numeric(as.character(tacsatp$VE_KW))<br>
<br>
if(class(dat) != "try-error"){<br>
counter &lt;- counter + 1<br>
allInts[[counter]] &lt;- dat<br>
}<br>
}<br>
}<br>
tacsatp &lt;- do.call(rbind,allInts)<br>
save(tacsatp,file=file.path(outPath,"tacsatSweptArea.RData"))<br>
</code></pre>

<h3>Link habitat map to vms</h3>

For some regions, habitat maps are available as shapefiles. If these datasets<br>
exist, we can read them into R and overlay the VMS points with the habitat map.<br>
Through this overlap we can identify a habitat type with VMS ping.<br>
<br>
First we load the habitat map and load the tacsat data with swept area indication.<br>
Thereafter, we convert the tacsat data into a 'coordinate' dataset which we then<br>
can overlay with the habitat map. The result of that exercise is a vertor which<br>
indicates which habitat can be associated with each VMS ping. Thereafter we assign<br>
this habitat type to the tacsat file.<br>
<br>
<b>Via shapefile</b>

<pre><code><br>
library(maptools)<br>
<br>
# load a habitat map shape file (Baltic)<br>
habitat_map           &lt;- readShapePoly(file.path(polPath,"sediment_lat_long"),<br>
proj4string=CRS("+proj=longlat +ellps=WGS84"))<br>
# load a habitat map shape file (North Sea)<br>
habitat_map           &lt;- readShapePoly(file.path(polPath,"ModelledSeabedHabitats"),<br>
proj4string=CRS("+proj=longlat +ellps=WGS84"))<br>
<br>
# get 'tacsatp' with all data<br>
load(file.path(outPath, "tacsatSweptArea.RData"))<br>
<br>
# Turn the habitat map into spatial polygons<br>
sp &lt;- SpatialPolygons(habitat_map@polygons)<br>
proj4string(sp) &lt;-  CRS("+proj=longlat +ellps=WGS84")<br>
<br>
#Turn the VMS point into a spatial points object<br>
spo             &lt;- SpatialPoints(coordinates(data.frame(SI_LONG=tacsatp$SI_LONG,<br>
SI_LATI=tacsatp$SI_LATI)))<br>
proj4string(spo) &lt;-  CRS("+proj=longlat +ellps=WGS84")<br>
<br>
# Use the magic 'over' function to match habitat type with VMS position<br>
idx &lt;- over(spo,sp)<br>
<br>
#- Baltic<br>
tacsatp$SUBSTRATE &lt;- habitat_map$BAL_CODE[idx]<br>
<br>
#- North Sea<br>
tacsatp$SUBSTRATE &lt;- habitat_map$substrate[idx]<br>
<br>
# Make a plot of the result (baltic)<br>
plot(habitat_map, xlim=c(11,14), ylim=c(55,56))<br>
axis(1) ; axis(2, las=2) ; box()<br>
points(tacsatp[, c("SI_LONG","SI_LATI")], col=tacsatp$SUBSTRATE, pch=".")<br>
<br>
# Make a plot of the results (North Sea)<br>
plot(1,1,col="white",xlim=c(-4,10), ylim=c(50,60),xlab="Longitude",ylab="Latitude")<br>
plot(habitat_map, add=T,border="grey")<br>
box()<br>
points(tacsatp[, c("SI_LONG","SI_LATI")], col=tacsatp$SUBSTRATE, pch=".")<br>
<br>
# Save the plot<br>
savePlot(filename=file.path(outPath, "VMSpingsAttachedToSedimentMap.jpeg"), type="jpeg")<br>
</code></pre>
The result might look something like this:<br>
<img src='http://vmstools.googlecode.com/svn/wiki/figures/habitatMapNorthSea.png' />

Here we've used shapefiles of habitat type and linked VMS points, but the habitat<br>
maps might be available in raster types only (already gridded files such as a TIF)<br>
<br>
<b>Via raster</b>

<pre><code><br>
<br>
sh_coastlines            &lt;- readShapePoly(file.path(polPath,"francois_EU"))<br>
<br>
## use point-raster overlay.......<br>
library(raster)<br>
landscapes       &lt;- raster(file.path(polPath, "landscapes.tif"))    # probably need an update of rgdal here....<br>
newproj          &lt;- "+proj=longlat +datum=WGS84"<br>
landscapes_proj  &lt;- projectRaster(landscapes, crs=newproj)<br>
<br>
save(landscapes_proj, file=file.path(polPath, "landscapes_proj.RData"))<br>
<br>
load(file.path(polPath, "landscapes_proj.RData"))<br>
<br>
load(file.path(outPath, "tacsatSweptArea.RData")) # get 'tacsatp' with all data<br>
#....or load only one instance eg load("C:\\merging\\BENTHIS\\outputs\\interpolated\\tacsatSweptArea_DNK000005269_OTB.RData"))<br>
#tacsatp &lt;- tacsatInt_gr_vid<br>
<br>
<br>
coord &lt;- cbind(x=anf(tacsatp$SI_LONG), y=anf(tacsatp$SI_LATI))<br>
<br>
dd &lt;- extract (landscapes_proj, coord[,1:2]) # get the landscape on the coord points!<br>
<br>
coord &lt;- cbind(coord,  landscapes_code=cut(dd, breaks=c(0,100,200,300,400,500,600)))<br>
<br>
tacsatp &lt;- cbind(tacsatp, landscapes_code= coord[,'landscapes_code'])<br>
<br>
# plot and save...<br>
plot(landscapes_proj, xlim=c(10,14), ylim=c(54.5,56.5))<br>
plot(sh_coastlines,  xlim=c(10,14), ylim=c(54.5,56.5), add=TRUE)  # better for plotting the western baltic sea coastline!<br>
points(coord[,"x"], coord[,"y"], col=coord[,"landscapes_code"], pch=".", cex=1)<br>
<br>
# save<br>
savePlot(filename=file.path(outPath, "VMSpingsAttachedToLandscapeMap.jpeg"), type="jpeg")<br>
<br>
</code></pre>




<ul><li>severity of activity: read in a look-up table (from interviews / experimental studies)</li></ul>

This should make use of the info we know about sub-component of the gears. Look at the Table 6.1 in the Deliverable 1.1b, the four gear types are described and broken down into individual components. The gear-seabed contact by component is conceptualized in “gear footprints” of each gear. These footprints form the basis of an industry questionnaire designed to deliver information of the dimensions of the individual gear components for those gears currently in use in European and Turkish fisheries.<br>
<br>
This table is about to be produced, but for the time being hereafter a fake example below:<br>
<br>
<br>
<pre><code><br>
## by gear/metier by habitat<br>
## (i.e. a mutiplying factor in relative terms)<br>
<br>
# create a fake input file to show the required format<br>
dd &lt;- tacsatp[!duplicated(data.frame(tacsatp$VE_REF,tacsatp$LE_GEAR, tacsatp$LE_MET)), ]<br>
fake_gear_metier_habitat_severity_table &lt;- dd[,c('VE_REF', 'LE_GEAR', 'LE_MET')]<br>
fake_gear_metier_habitat_severity_table &lt;- cbind(fake_gear_metier_habitat_severity_table, HAB_SEVERITY=1)<br>
gear_metier_habitat_severity_table  &lt;- fake_gear_metier_habitat_severity_table<br>
gear_metier_habitat_severity_table  &lt;-   gear_metier_habitat_severity_table[complete.cases( gear_metier_habitat_severity_table),]<br>
save(gear_metier_habitat_severity_table,   file=paste(dataPath,"gear_metier_habitat_severity_table.RData",   sep=""))<br>
<br>
# load a table for HAB_SEVERITY per vid LE_REF per gr LE_GEAR per met LE_MET<br>
load(file.path(dataPath, "gear_metier_habitat_severity_table.RData"))<br>
tacsatp &lt;- merge(tacsatp, gear_metier_habitat_severity_table)<br>
save(tacsatp,   file=paste(outPath,"tacsatMergedHabSeverity.RData",   sep=""))<br>
</code></pre>


<ul><li>pressure: weigh the swept area by the severity</li></ul>

<pre><code><br>
tacsatp$pressure &lt;- tacsatp$HAB_SEVERITY * tacsatp$SWEPT_AREA_KM2<br>
</code></pre>



<ul><li>define grid<br>
<ul><li>mid-point polygon distance distribution (to help defining a suitable grid resolution)<br>
<pre><code><br>
</code></pre></li></ul></li></ul>








<ul><li>random process of fishing activity (to help defining a suitable grid resolution)</li></ul>

<pre><code><br>
</code></pre>






<ul><li>maps + sum / average different grids<br>
</li></ul><blockquote>using a quick gridding code at various resolution.</blockquote>

<ul><li>For example, grid the swept area, or the number of hauls, or the fishing pressure, etc.<br>
<pre><code><br>
##-----------------------------------------------------------------------------<br>
#<br>
<br>
## GRIDDING (IN DECIMAL DEGREES OR IN UTM COORD)<br>
# using a quick gridding code at various resolution.<br>
<br>
# For example, grid the swept area, or the number of hauls, or the fishing pressure, etc.<br>
sh1 &lt;- readShapePoly(file.path(polPath,"francois_EU")) # coastline<br>
<br>
<br>
<br>
### the swept area-------------------------------------------<br>
<br>
## user selection here----<br>
what                 &lt;- "SWEPT_AREA_KM2"<br>
#what                &lt;- "HL_ID"<br>
#what                &lt;- "effort_days"<br>
#what                &lt;- "effort_KWdays"<br>
#what                &lt;- "pressure"<br>
#a_func              &lt;- function(x) {unique(length(x))} # for HL_ID, to be tested.<br>
a_func               &lt;- "sum"<br>
is_utm               &lt;- FALSE<br>
all_gears            &lt;- sort(unique(tacsatp$LE_GEAR))<br>
towed_gears          &lt;- c('OTB', 'TBB', 'PTB', 'PTM', 'DRB')  # TO DO: list to be checked<br>
passive_gears        &lt;- all_gears[!all_gears %in% towed_gears]<br>
we &lt;- 10; ea &lt;- 13; no &lt;- 59; so &lt;- 55;<br>
##------------------------<br>
<br>
# subset for relevant fisheries<br>
this            &lt;- tacsatp [tacsatp$LE_GEAR %in% towed_gears, ]<br>
<br>
# restrict the study area<br>
# (it is likely that the same bounding box should be used when stacking different layers e.g. from different countries)<br>
this &lt;- this[this$SI_LONG&gt;we &amp; this$SI_LONG&lt;ea &amp; this$SI_LATI&gt;so &amp; this$SI_LATI&lt;no,]<br>
<br>
# grid the data (in decimal or in UTM)<br>
if(is_utm){<br>
dx &lt;- 0.0002 # 5km x 5km<br>
# convert to UTM<br>
library(sp)<br>
library(rgdal)<br>
SP &lt;- SpatialPoints(cbind(as.numeric(as.character(this$SI_LONG)), as.numeric(as.character(this$SI_LATI))),<br>
proj4string=CRS("+proj=longlat +datum=WGS84"))<br>
this &lt;- cbind(this,<br>
spTransform(SP, CRS("+proj=utm  +ellps=intl +zone=32 +towgs84=-84,-107,-120,0,0,0,0,0")))    # convert to UTM<br>
this            &lt;- this [, c('SI_LONG', 'SI_LATI', 'SI_DATE', 'coords.x1', 'coords.x2', what)]<br>
this$round_long &lt;- round(as.numeric(as.character(this$coords.x1))*dx)<br>
this$round_lat  &lt;- round(as.numeric(as.character(this$coords.x2))*dx)<br>
this            &lt;- this[, !colnames(this) %in% c('coords.x1', 'coords.x2')]<br>
this$cell       &lt;- paste("C_",this$round_long,"_", this$round_lat, sep='')<br>
this$xs         &lt;- (this$round_long/(dx))<br>
this$ys         &lt;- (this$round_lat/(dx))<br>
<br>
}  else {<br>
dx &lt;- 20 # 0.05 degree<br>
this &lt;- this [, c('SI_LONG', 'SI_LATI', 'SI_DATE', what)]<br>
this$round_long &lt;- round(as.numeric(as.character(this$SI_LONG))*dx*2) # 0.1<br>
this$round_lat  &lt;- round(as.numeric(as.character(this$SI_LATI))*dx)   # 0.05<br>
this$cell       &lt;- paste("C_",this$round_long,"_", this$round_lat, sep='')<br>
this$xs         &lt;- (this$round_long/(dx*2))<br>
this$ys         &lt;- (this$round_lat/(dx))<br>
}<br>
# if the coordinates in decimal then dx=20 corresponds to grid resolution of 0.05 degrees<br>
# i.e. a 3� angle = 3nm in latitude but vary in longitude (note that a finer grid will be produced if a higher value for dx is put here)<br>
<br>
colnames(this) &lt;- c('x', 'y', 'date', 'what', 'round_long', 'round_lat', 'cell', 'xs', 'ys')<br>
<br>
<br>
# retrieve the geo resolution in degree, for info<br>
#long &lt;- seq(1,15,by=0.01)<br>
#res_long &lt;- diff( long [1+which(diff(round(long*dx)/dx/2)!=0)] )<br>
#res_lat &lt;- diff( long [1+which(diff(round(long*dx/2)/dx)!=0)] )<br>
#print(res_long) ; print(res_lat)<br>
<br>
<br>
# a quick gridding method...<br>
background &lt;- expand.grid(<br>
x=0,<br>
y=0,<br>
date=0,<br>
what=0,<br>
round_long=seq(range(this$round_long)[1], range(this$round_long)[2], by=1),<br>
round_lat=seq(range(this$round_lat)[1], range(this$round_lat)[2], by=1),<br>
cell=0,<br>
xs=0,<br>
ys=0<br>
)<br>
this &lt;- rbind(this, background)<br>
the_points &lt;- tapply(this$what,<br>
list(this$round_lat, this$round_long), a_func)<br>
<br>
xs &lt;- (as.numeric(as.character(colnames(the_points)))/(dx*2))<br>
ys &lt;- (as.numeric(as.character(rownames(the_points)))/(dx))<br>
<br>
the_breaks &lt;-  c(0, (1:12)^3.5 ) # to be decided...<br>
graphics:::image(<br>
x=xs,<br>
y=ys,<br>
z= t(the_points)  ,<br>
breaks=c(the_breaks),<br>
col = terrain.colors(length(the_breaks)-1),<br>
useRaster=FALSE,<br>
xlab="",<br>
ylab="",<br>
axes=FALSE,<br>
xlim=range(xs), ylim=range(ys),<br>
add=FALSE<br>
)<br>
title("")<br>
<br>
# land<br>
sh1 &lt;- readShapePoly(file.path(polPath,"francois_EU"),  proj4string=CRS("+proj=longlat +datum=WGS84"))<br>
if(is_utm) sh1 &lt;- spTransform(sh1, CRS("+proj=utm  +ellps=intl +zone=32 +towgs84=-84,-107,-120,0,0,0,0,0"))<br>
plot(sh1, add=TRUE, col=grey(0.7))<br>
<br>
legend("topright", fill=terrain.colors(length(the_breaks)-1),<br>
legend=round(the_breaks[-1],1), bty="n", cex=0.8, ncol=2, title="")<br>
box()<br>
axis(1)<br>
axis(2, las=2)<br>
<br>
if(is_utm){<br>
mtext(side=1, "Eastings", cex=1, adj=0.5, line=2)<br>
mtext(side=2, "Northings", cex=1, adj=0.5, line=2)<br>
} else{<br>
mtext(side=1, "Longitude", cex=1, adj=0.5, line=2)<br>
mtext(side=2, "Latitude", cex=1, adj=0.5, line=2)<br>
#points (tacsatp [,c('SI_LONG', 'SI_LATI')], pch=".", col="white")<br>
}<br>
<br>
<br>
# save<br>
savePlot(filename=file.path(outPath, "GriddedSweepAreaExample.jpeg"), type="jpeg")<br>
<br>
<br>
# export the quantity per cell and date<br>
library(data.table)<br>
DT                &lt;-  data.table(this)<br>
qu                &lt;-  quote(list(sum(what)))<br>
quantity_per_cell &lt;- DT[,eval(qu), by=list(cell, xs,ys)]<br>
quantity_per_date &lt;- DT[,eval(qu), by=list(date, xs,ys)]<br>
quantity_per_cell_date &lt;- DT[,eval(qu), by=list(cell,date, xs,ys)]<br>
quantity_per_cell_date &lt;- as.data.frame(quantity_per_cell_date)<br>
colnames(quantity_per_cell_date)   &lt;- c('cell','date', 'xs','ys', 'quantity')<br>
rm(DT) ; gc(reset=TRUE)<br>
save(quantity_per_cell_date, res_long, res_lat,  we, ea, no, so, file=file.path(outPath,"quantity_per_cell_date.RData") )<br>
<br>
<br>
# export the cumul per cell and date<br>
quantity_per_cell_date &lt;- orderBy(~date, data=quantity_per_cell_date)<br>
quantity_cumul_per_cell_date &lt;- do.call("rbind", lapply(<br>
split(quantity_per_cell_date, f=quantity_per_cell_date$cell),<br>
function(x){<br>
x$quantity &lt;- cumsum(x$quantity)<br>
x<br>
})  )<br>
# check the cumul on a given cell<br>
# head(quantity_cumul_per_cell_date[quantity_cumul_per_cell_date$cell=="C_197_2722",])<br>
save(quantity_cumul_per_cell_date, res_long, res_lat,  we, ea, no, so, file=file.path(outPath,"quantity_cumul_per_cell_date.RData") )<br>
<br>
<br>
<br>
<br>
<br>
</code></pre></li></ul>


<blockquote>[[Image:<a href='http://code.google.com/p/vmstools/source/browse/wiki/GriddedSweepAreaExample.jpeg|thumb|center|upright=2.0|alt=Gridded'>http://code.google.com/p/vmstools/source/browse/wiki/GriddedSweepAreaExample.jpeg|thumb|center|upright=2.0|alt=Gridded</a> Swept Area.<br>
Gridded Swept Area example in km2 for one vessel in the Baltic Sea. | Original VMS pings are overlaid]]</blockquote>


or using createGrid()<br>
<pre><code><br>
</code></pre>

...and stack some layers making use of mosaic() from raster layers.<br>
<br>
<pre><code><br>
# from ascgridfile eg. returned by VMSGridCreate()<br>
# TO DO:<br>
dd &lt;-<br>
nn &lt;-<br>
fun1  &lt;-  function(x)  { x[ is.na(x)]  &lt;- 0;  return(x)  }            # remove NA such zero<br>
dd &lt;-  calc(dd,  fun1)<br>
nn &lt;-  calc(nn,  fun1)<br>
a_stack &lt;-mosaic(dd, nn, na.rm=TRUE, fun = sum, keepres=TRUE)<br>
</code></pre>


<ul><li>non-vms maps. For non VMS-equipped vessels, the mapping of the fishing effort<br>
</li></ul><blockquote>can be generated at the ICES rectangle resolution as a proxy of the fishing pressure on the benthic habitats</blockquote>

<pre><code><br>
# TO DO from e.g. plotTools()<br>
#....<br>
# TO DO: stack several layers with various geographical resolutions<br>
</code></pre>





<ul><li>testing</li></ul>

<pre><code><br>
</code></pre>


