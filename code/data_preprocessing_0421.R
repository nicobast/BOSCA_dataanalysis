#COMMENT: Preprocessing has to be done with R-version < 4!

rm(list=ls())

#### 0.1. Version Control ####
#17.12.2020
#-left out validity code >= 2 exclusion in preprocessing
#-left in data exclusion for gazepoints >1 & <0 (not on screen)
#-rpd defined by PD / tonic PD across PLR conditions
#-plotted rpd over task duration per condition (black/white) and group (ASD/TD)
#-start values of amplitude and latency codes changed to tonic PD across conditions in PLR task
#05.01.2021
#-changed data frame for analyses so that trials can be differentiated (not only conditions)
#-added variables: id, condition, trial number
#-added LMMs for analyses
#08.01.2021
#-NA exclusion now only after pd.colscreen creation
#-added trial number variable without gaps that exist due to NA exclusion
#-deletet "old" trial number calculation
#11.01.2021
#-loading of demographics
#-included age variable and merged to PLR_measures
#12.01.2021
#-included sex variable
#-adapted LMMs for analyses

#---------------------------------------------------------------------------------------------------#

#COMMENT: Preprocessing has to be done with R-version < 4!

#### 0.2. Load Packages ####
require(ggplot2)
require(zoo)
require(reshape2)
require(lme4)
require(lmerTest)
require(MatchIt)
require(emmeans)

#-----------------------------------------------------------------------------------------------------#

#### 1. Load Data ####
#datapath<-"D:/Polzer/Daten/"
datapath<-'~/PowerFolders/data_AFFIP/data' #preprocessed on LINUX machine

#read data from datapath and store in according objects
data.files<-list.files(path=datapath,full.names=T)
data.time<-data.files[grep('_timestamps',data.files)]
data.events<-data.files[grep('_event',data.files)]
data.files<-data.files[grep('_gazedata',data.files)]

#reduce datafile to relevant datasets
data.files<-data.files[-grep('RECOVER',data.files)] #exclude recover files
data.files<-data.files[c(grep('t2_gazedata.csv',data.files),grep('K_gazedata.csv',data.files))] #exclude other test points
data.time<-data.time[-grep('RECOVER',data.time)] #exclude recover files
data.time<-data.time[c(grep('t2_timestamps.csv',data.time),grep('K_timestamps.csv',data.time))] #exclude other test points
data.events<-data.events[-grep('RECOVER',data.events)] #exclude recover files
data.events<-data.events[c(grep('t2_event.csv',data.events),grep('K_event.csv',data.events))] #exclude other test points
# -> all three file lists (gazedata, time, events) should have the same length

#read all gaze data into list class
df.list.data<-lapply(data.files,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
df.list.time<-lapply(data.time,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
df.list.events<-lapply(data.events,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)

#relabel data
name.labels<-c('eyepos.X_L','eyepos.Y_L','eyepos.Z_L','releyepos.X_L','releyepos.Y_L','eyepos.Z_L','gazepos2D.X_L','gazepos2D.Y_L','gazepos3D.X_L','gazepos3D.Y_L','gazepos3D.Z_L','pupildil_L','validcode_L','eyepos.X_R','eyepos.Y_R','eyepos.Z_R','releyepos.X_R','releyepos.Y_R','eyepos.Z_R','gazepos2D.X_R','gazepos2D.Y_R','gazepos3D.X_R','gazepos3D.Y_R','gazepos3D.Z_R','pupildil_R','validcode_R')
for(i in 1:length(df.list.data)){names(df.list.data[[i]])<-name.labels}
for(i in 1:length(df.list.time)){names(df.list.time[[i]])<-'timestamp'}
for(i in 1:length(df.list.events)){names(df.list.events[[i]])<-c('event','ev.ts')}

#add timestamp to gaze data - only possible for same length files
df.list<-mapply(data.frame,df.list.time,df.list.data,SIMPLIFY = F)

#Add event data to gaze data
##Function to map events to according timestamps:
eventfunc<-function(x,y,z,i){
  event<-which(x[,1]>=y[i,2])
  z[event]<-y[i,1]
  return(z)}
##apply this function to all list elements of df.list (for loop in for loop)
df<-list()
for(i in 1:length(df.list.data)){
  eventlog<-rep(NA,nrow(df.list.events[[i]]))
  for(j in 1:nrow(df.list.events[[i]])){eventlog<-eventfunc(df.list[[i]],df.list.events[[i]],eventlog,j)}
  eventlog<-as.factor(eventlog) #changes the values
  levels.eventlog<-levels(df.list.events[[i]]$event)[as.numeric(levels(eventlog))] #retrieve those levels that are in eventlog
  levels(eventlog)<-levels.eventlog #apply retrieved events
  df[[i]]<-cbind(eventlog,df.list[[i]])
}


#add names to list
id.names<-substr(data.files,nchar(datapath)+2,nchar(datapath)+6) #08.07.19: now independent of relative to datapath
names(df)<-id.names

#remove lists: clear ram
rm(df.list.data, df.list.events, df.list.time)


#---------------------------------------------------------------------------------------------------#

#### 2. Data Preprocessing ####

#create index variable (specific to eventlog: INDEX.EVENT, TS.EVENT) #
event.change<-lapply(df,function(x){which(diff(as.numeric(x$eventlog))!=0)}) #index event change per participant
length.datasets<-sapply(df,nrow)
event.change<-mapply(function(x,y){x<-c(0,x,y)},event.change,length.datasets) #add start and end event to index
event.change<-lapply(event.change,diff) #create a difference value - LEONIE: problem - difference row wise?
#event.change<-diff(event.change) # LEONIE eingefÃ¼gt
#event.change<-split(event.change, rep(1:ncol(event.change), each = nrow(event.change))) # LEONIE eingefÃ¼gt: turn matrix into list (one column = one list element)
index.event<-lapply(event.change,function(x){rep(seq_along(x),times=x)})
ts.event<-lapply(index.event,function(x){as.numeric(do.call(c,by(x,x,seq_along, simplify=F)))})

df<-mapply(function(x,y,z){data.frame(x,y,z)},df,index.event,ts.event,SIMPLIFY = F)
#look at first eventlog change to check if it worked
#df[["046_t"]][280:284,c("eventlog", "y", "z")]
#y=index.event, z=ts.event - variables are renamed below (after turning df into data frame)


###PD Preprocessing (Kret, 2018)###
func_pd_preprocess<-function(x){

  #define variables
  Left_Diameter<-x$pupildil_L
  Right_Diameter<-x$pupildil_R
  RemoteTime<-x$timestamp

  #constant for MAD caluclation
  constant<-3 ##--> if change speed is higher than constant * median change --> values are excluded
  #constant<-3 #default value

  # STEP 1 - exclude invalid data ####
  pl <- ifelse((Left_Diameter<2|Left_Diameter>8), NA, Left_Diameter)
  pr <- ifelse((Right_Diameter<2|Right_Diameter>8), NA, Right_Diameter)
  #table(is.na(pl))
  #table(is.na(pr))
  # STEP 2 - filtering ####
  ## A) normalized dilation speed, take into account time jumps with Remotetimestamps: ####
  #maximum change in pd compared to last and next pd measurement
  #Left
  pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
  pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
  pl.speed1<-c(NA,pl.speed1)
  pl.speed2<-c(rev(pl.speed2),NA)
  pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
  rm(pl.speed1,pl.speed2)
  #Right
  pr.speed1<-diff(pr)/diff(RemoteTime)
  pr.speed2<-diff(rev(pr))/diff(rev(RemoteTime))
  pr.speed1<-c(NA,pr.speed1)
  pr.speed2<-c(rev(pr.speed2),NA)
  pr.speed<-pmax(pr.speed1,pr.speed2,na.rm=T)
  rm(pr.speed1,pr.speed2)
  #median absolute deviation -SPEED
  #constant<-3
  pl.speed.med<-median(pl.speed,na.rm=T)
  pl.mad<-median(abs(pl.speed-pl.speed.med),na.rm = T)
  pl.treshold.speed<-pl.speed.med+constant*pl.mad #treshold.speed units are mm/microsecond
  #plot(abs(pl.speed))+abline(h=pl.treshold.speed)
  pr.speed.med<-median(pr.speed,na.rm=T)
  pr.mad<-median(abs(pr.speed-pr.speed.med),na.rm = T)
  pr.treshold.speed<-pr.speed.med+constant*pr.mad #treshold.speed units are mm/microsecond
  #plot(abs(pr.speed))+abline(h=pr.treshold.speed)
  #correct pupil dilation for speed outliers
  pl<-ifelse(abs(pl.speed)>pl.treshold.speed,NA,pl)
  pr<-ifelse(abs(pr.speed)>pr.treshold.speed,NA,pr)
  ## B) delete data around blinks - not applied ####
  ## C) normalized dilation size - median absolute deviation -SIZE ####
  #applies a two pass approach
  #first pass: exclude deviation from trend line derived from all samples
  #second pass: exclude deviation from trend line derived from samples passing first pass
  #-_> reintroduction of sample that might have been falsely excluded due to outliers
  #estimate smooth size based on sampling rate
  smooth.length<-150 #measured in ms
  #take sampling rate into account (300 vs. 120):
  #smooth.size<-round(smooth.length/mean(diff(RemoteTime)/1000)) #timestamp resolution in microseconds
  smooth.size<-round(smooth.length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
  is.even<-function(x){x%%2==0}
  smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)
  #Left
  pl.smooth<-na.approx(pl,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - FIRST pass
  pl.pass1<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass1<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  #Left
  pl.smooth<-na.approx(pl.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - SECOND pass
  pl.pass2<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass2<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  pl<-pl.pass2
  pr<-pr.pass2

  ## D) sparsity filter - not applied ####
  # STEP 3 - processing valid samples  ####
  #take offset between left and right into account
  pd.offset<-pl-pr
  pd.offset<-na.approx(pd.offset,rule=2)
  #mean pupil dilation across both eyes
  pl <- ifelse(is.na(pl)==FALSE, pl, pr+pd.offset)
  pr <- ifelse(is.na(pr)==FALSE, pr, pl-pd.offset)
  pd <- (pl+pr)/2
  # end of function --> return ####
  #detach(x)
  return(pd)
}

pd.list<-lapply(df, func_pd_preprocess)


#create long format variables: id, group, ts, pd (preprocessed)
id<-rep(names(pd.list),as.numeric(lapply(pd.list,length)))
group<-ifelse(grepl('_K',id),'TD','ASD') #define group
pd<-as.numeric(do.call(c,pd.list))
ts<-lapply(df,function(x){x<-x$timestamp})
ts<-lapply(ts,function(x){abs(head(x,n=1)-x)/1000000})
ts<-as.numeric(do.call(c,ts))

#create main data frame
df<-do.call(rbind,df)
df<-data.frame(ts,id,group,pd,df)

##rename event variables
names(df)[names(df) == "y"]<-'index.event'
names(df)[names(df) == "z"]<-'ts.event'


# EXTRACT TOTAL EVENT DURATION (samples_per_event) #
samples_per_event<-as.numeric(with(df,by(ts.event,droplevels(interaction(eventlog,index.event)),max,na.rm=T))) #maximum ts.event per distinct event per participant
event_label<-as.factor(with(df,by(eventlog,droplevels(interaction(eventlog,index.event)),head,n=1))) #extract eventlog label
levels(event_label)<-levels(df$eventlog) #copy the levels labels
samples_per_event<-by(samples_per_event,event_label,median,na.rm=T) #median samples per eventlog
names(samples_per_event)<-levels(event_label) #copy the levels labels

#control for GAZE BEHAVIOR on screen
##TO DO: compare data quantity with and without excluded data points
attach(df)
#exclude implausible values
xl <- ifelse((gazepos2D.X_L<0|gazepos2D.X_L>1), NA, gazepos2D.X_L)
xr <- ifelse((gazepos2D.X_R<0|gazepos2D.X_R>1), NA, gazepos2D.X_R)
yl <- ifelse((gazepos2D.Y_L<0|gazepos2D.Y_L>1), NA, gazepos2D.Y_L)
yr <- ifelse((gazepos2D.Y_R<0|gazepos2D.Y_R>1), NA, gazepos2D.Y_R)
#take offset between left and right into account
x.offset<-xl-xr
x.offset<-na.approx(x.offset,rule=2)
y.offset<-yl-yr
y.offset<-na.approx(y.offset,rule=2)
#mean gaze across both eyes
xl <- ifelse(is.na(xl)==FALSE, xl, xr+x.offset)
xr <- ifelse(is.na(xr)==FALSE, xr, xl-x.offset)
yl <- ifelse(is.na(yl)==FALSE, yl, yr+y.offset)
yr <- ifelse(is.na(yr)==FALSE, yr, yl-y.offset)
gazepos.x<-(xl+xr)/2
gazepos.y<-(yl+yr)/2
#add the gaze position data to df
df<-data.frame(df,gazepos.x,gazepos.y)
#exclude data, for which gaze position or PD = NA (TO DO: check if there are PD values for which gaze position = NA)
df<-df[!is.na(df$gazepos.x) & !is.na(df$gazepos.y),]
df<-df[!is.na(df$pd),]

### SAVE PREPROCESSED DATA ####

object.size(df)

save(df,file='~/PowerFolders/data_AFFIP/data_preprocessed_300421.Rdata')


