# HEADER ####

## Script purpose: PREPROCESS JOINT ATTENTION DATA - BOSCA PROJECT
##
##
## Author: Nico Bast
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
## Email: nico.bast@kgu.de
##
##  INFO:
##
##  - extracted stimulus onset info from splitting videos intro frames with custom python script
##    e.g.: CMD line = C:\Users\Nico\PowerFolders\project_video_salience\code\python_code_salience_extraction\split_video_to_images.py "neutral_rabbit_left" "C:\Users\Nico\PowerFolders\ETbat\Jointattention\stimuli\blockB"
##    --> tested with four videos: length ~ 550 frames = 11 s, face onset ~ 64-72 frames ~ 1.35 sec, cueing onset ~ 165-180 frames ~  3.46 sec, stimulus_shaking ~ 150 frames = 3 sec
##
##  - SYSTEM REQUIREMENTS: at least 64GB of RAM to read data of 314 sessions (13.03.2024)
##  - Line 1-170: is independent of the task and describes general preprocessing of BOSCA data
##  - Line 1100: demographic data is preprocessed by script sample_overview_Jan2024.Rmd - this data is used here

# SETUP ####

#check for OS --> define home path (script independent of OS)
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

#project path
project_path<-'/PowerFolders/project_BOSCA_battery'

#data path
data_path<-'/PowerFolders/data_AFFIP/data'

#user-defined functions
fun_rename<-function(x,variable_position,new_name){
  names(x)[variable_position]<-new_name
  return(x)}


# Load Packages ####
require(zoo) #na.approx
require(data.table) #fread uses parallelization and thus mucd faster than read.csv

#analysis, visualization
require(ggplot2)
require(lme4) #linear mixed models
require(emmeans) #linear mixed models
require(lmerTest) #linear mixed models
require(ggpubr) #significance bars in figures
require(pbapply) #lapply with progress bar

#--> loaded packages
sessionInfo()

# Load Data ####
start_time <- Sys.time()

datapath<-paste0(home_path,data_path) #preprocessed on LINUX machine

#datapath<-"G:/BACKUP_Polzer_ETBattery/data_AFFIP/data"

#read data from datapath and store in according objects
data.files<-list.files(path=datapath,full.names=T)
data.time<-data.files[grep('_timestamps',data.files)]
data.events<-data.files[grep('_event',data.files)]
data.files<-data.files[grep('_gazedata',data.files)]

#reduce datafile to relevant datasets
data.files<-data.files[-grep('RECOVER',data.files)] #exclude recover files
data.time<-data.time[-grep('RECOVER',data.time)] #exclude recover files
data.events<-data.events[-grep('RECOVER',data.events)] #exclude recover files
data.files<-data.files[-grep('test',data.files)] #exclude recover files
data.time<-data.time[-grep('test',data.time)] #exclude recover files
data.events<-data.events[-grep('test',data.events)] #exclude recover files

#remove double entry of 094_t6
#data.events<-data.events[-134]
data.events<-data.events[data.events!="C:/Users/Nico/PowerFolders/data_AFFIP/data/094_t6__event.csv"]


##remove participants that have not all data (events,gazedata,time)
id_names_time<-substr(data.time,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.time)-nchar("_timestamps.csv"))
id_names_files<-substr(data.files,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.files)-nchar("_gazedata.csv"))
id_names_events<-substr(data.events,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.events)-nchar("_event.csv"))

id_names_files[!(id_names_files %in% id_names_time)]
id_names_events[!(id_names_events %in% id_names_time)]
##--> empty

data.files<-data.files[data.files!="C:/Users/Nico/PowerFolders/data_AFFIP/data/066_fu2_gazedata.csv"]
# data.events<-data.events[data.events!="C:/Users/Nico/PowerFolders/data_AFFIP/data/066_fu2_event.csv"]

#sorting - that files names match (k_fu introduced different sorting in data.files versus data.events)
###--> sort according to participant id
order_data<-order(substr(data.files,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.files)-nchar("_gazedata.csv")))
order_events<-order(substr(data.events,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.events)-nchar("_event.csv")))
order_time<-order(substr(data.time,nchar("C:/Users/Nico/PowerFolders/data_AFFIP/data/")+1,nchar(data.time)-nchar("_timestamps.csv")))

data.files<-data.files[order_data]
data.events<-data.events[order_events]
data.time<-data.time[order_time]

#do files match? --> visual inspection
cbind(substr(data.files,44,57),substr(data.events,44,54),substr(data.time,44,54))

#READ
#read csv is very slow
df.list.data<-list(0)
for(i in 1:length(data.files)){
  df.list.data[[i]]<-fread(data.files[i])
  print(paste0('read: ',i))
}
###16 gb fails after 118 entries --> requires at least 32GB (June 2022)
###--> fine for 64gb of memory (March 2024 --> 314 data entries)
### for loop is slower, but does not kill memory

#df.list.data<-lapply(data.files,read.csv)
df.list.time<-pblapply(data.time,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
df.list.events<-pblapply(data.events,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
###--> may takes some time

end_time <- Sys.time()
end_time - start_time

## --> RELABEL + CONCATENATE data ####
start_time <- Sys.time()

#relabel data
name.labels<-c('eyepos.X_L','eyepos.Y_L','eyepos.Z_L','releyepos.X_L','releyepos.Y_L','eyepos.Z_L','gazepos2D.X_L','gazepos2D.Y_L','gazepos3D.X_L','gazepos3D.Y_L','gazepos3D.Z_L','pupildil_L','validcode_L','eyepos.X_R','eyepos.Y_R','eyepos.Z_R','releyepos.X_R','releyepos.Y_R','eyepos.Z_R','gazepos2D.X_R','gazepos2D.Y_R','gazepos3D.X_R','gazepos3D.Y_R','gazepos3D.Z_R','pupildil_R','validcode_R')
for(i in 1:length(df.list.data)){names(df.list.data[[i]])<-name.labels}
for(i in 1:length(df.list.time)){names(df.list.time[[i]])<-'timestamp'}
for(i in 1:length(df.list.events)){names(df.list.events[[i]])<-c('event','ev.ts')}

#add timestamp to gaze data - only possible for same length files
df.list<-pbmapply(data.frame,df.list.time,df.list.data,SIMPLIFY = F)

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
  print(paste0('matched: ',i))
}

#add names to list
#id.names<-substr(data.files,nchar(datapath)+2,nchar(datapath)+7) #08.07.19: now independent of relative to datapath
id.names<-substr(data.files,nchar(datapath)+2,nchar(data.files)-nchar('_gazedata.csv')) #06.03.23: amended with FUs

names(df)<-id.names

## --> SAVE MERGED DATA (TEMP) ####

#remove lists: clear ram
rm(df.list.data, df.list.events, df.list.time)

#save(df,file="F:/temp_data_AFFIP/all_data_merged_040722.Rdata")
#save(df,file="C:/Users/nico/Desktop/BOSCA_joint_attention_all_data_merged_060323.Rdata")
#save(df,file="C:/Users/nico/Desktop/BOSCA_joint_attention_all_data_merged_170423.Rdata")
save(df,file="C:/Users/nico/Desktop/BOSCA_all_data_merged_130324.Rdata")
###-->save takes very long

#--> now found on NAS

end_time <- Sys.time()
end_time - start_time

# --> ##load("F:\temp_data_AFFIP\all_data_merged_040722.Rdata")

#---------------------------------------------------------------------------------------#

## ----------------------- #####

## CREATE TIMESTAMP VARIABLE ####
###--> before task data is selected

#create long format variable: ts (in seconds format)
ts<-lapply(df,function(x){x<-x$timestamp})
ts<-lapply(ts,function(x){abs(head(x,n=1)-x)/1000000})

#add ts and change name of the variable
df<-pbmapply(cbind,df,ts,SIMPLIFY = FALSE) #simplify needs to be in capital letters
df<-pblapply(df,fun_rename,variable_position=29,new_name='ts')

## --> SELECT TASK DATA: JOINT ATTENTION ####

  #diagnostics - find right eventlogs
  table(df[[100]]['eventlog'])
  table(df[[220]]['eventlog'])
  list_eventlog<-lapply(df,function(x){x['eventlog']})
  df_eventlog<-dplyr::bind_rows(list_eventlog)
  table(droplevels(df_eventlog$eventlog[grepl('jointattention_',df_eventlog$eventlog)]))
  table(droplevels(df_eventlog$eventlog))

list_jointatt<-pblapply(df,function(x){x[grepl('jointattention_',x$eventlog),]})

#remove useless data
list_jointatt<-pblapply(list_jointatt,function(x){x[!grepl('jointattention_blockA',x$eventlog),]})
list_jointatt<-pblapply(list_jointatt,function(x){x[!grepl('jointattention_blockB',x$eventlog),]})

## --> CREATE VARIABLES: trial, ts.trial  ####

#FUNCTION: identify trials --> function that based on change to identify trials and add trial sequence variable
fun_define_trials<-function(block_data){

  #create index variable - indicates when event changes #
  event.change<-pblapply(block_data,function(x){which(diff(as.numeric(as.factor(x$eventlog)))!=0)}) #index event change per participant
  length.datasets<-sapply(block_data,nrow) #length of each data set (per participant)
  event.change<-pbmapply(function(x,y){x<-c(0,x,y)},event.change,length.datasets) #add start and end event to index
  event.change<-pblapply(event.change,diff) #create a difference value

  #create index for these new trials
  index.trial<-pblapply(event.change,function(x){rep(seq_along(x),times=x)})

  #sequence over each new trial
  ts.event<-pblapply(index.trial,function(x){as.numeric(do.call(c,by(x,x,seq_along, simplify=F)))}) #tts over each trial

  #add index and sequence data to list
  block_data<-pbmapply(function(x,y,z){data.frame(x,y,z)},block_data,index.trial,ts.event,SIMPLIFY = F)
  block_data<-pblapply(block_data,fun_rename,variable_position=30,new_name='index_trial')
  block_data<-pblapply(block_data,fun_rename,variable_position=31,new_name='ts_event')
  return(block_data)
}

#apply functions that add trial and trial sequence data
list_jointatt<-fun_define_trials(list_jointatt)

#drop participants that have more than 16 trials (recording error for two participants)
list_jointatt<-pblapply(list_jointatt,function(x){x[x$index_trial<=16,]})

#---------------------------------------------------------------------------------------------------#
# EYE TRACKING DATA PREPROCESSING ####

rm(ts)

#-- control screen attention
fun_screen_att<-function(x){

  attach(x)
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

  #remove outside screen
  gazepos.x<-ifelse(gazepos.x>1 | gazepos.x<0,NA,gazepos.x)
  gazepos.y<-ifelse(gazepos.y>1 | gazepos.y<0,NA,gazepos.y)

  #estimate center deviation
  center_deviation<-sqrt((gazepos.x-0.5)^2 + (gazepos.y-0.5)^2)

  x[,'gazepos.x']<-gazepos.x
  x[,'gazepos.y']<-gazepos.y
  x[,'center_dev']<-center_deviation

  # #exclude data, for which gaze position or PD = NA (TO DO: check if there are PD values for which gaze position = NA)
  # df<-df[!is.na(df$gazepos.x) & !is.na(df$gazepos.y),]
  detach(x)

  return(x)
}
list_jointatt<-pblapply(list_jointatt,fun_screen_att)

## -- retrieve screen distance ####
fun_screen_dist<-function(x){

  dist_L<-x$eyepos.Z_L #in mm from tracker
  dist_R<-x$eyepos.Z_R #in mm from tracker

  #exclude implausible values (smaller 500mm and larger 800mm is outside track box)
  dist_L<-ifelse(dist_L > 800 | dist_L < 500, NA, dist_L)
  dist_R<-ifelse(dist_R > 800 | dist_R < 500, NA, dist_R)

  #take offset between left and right into account
  offset<-dist_L-dist_R
  offset<-na.approx(offset,rule=2)

  #mean gaze across both eyes
  dist_L <- ifelse(is.na(dist_L)==FALSE, dist_L, dist_R+offset)
  dist_R <- ifelse(is.na(dist_R)==FALSE, dist_R, dist_R-offset)

  screen_dist<-(dist_L+dist_R)/2
  x[,'screen_dist']<-screen_dist

  return(x)

}
list_jointatt<-pblapply(list_jointatt,fun_screen_dist)

## -- drop unecessary data + particiapnts without data --> subsequent preprocessing requires far less RAM####
fun_required_necessary_data<-function(x){

  #drop raw eye tracking data
  x<-x[,!(grepl('gazepos2D',names(x)) | grepl('gazepos3D',names(x)) | grepl('eyepos',names(x)))]
  return(x)

}
list_jointatt<-pblapply(list_jointatt,fun_required_necessary_data)

## -- remove participants without data
dropped_participants<-which(sapply(list_jointatt,function(x){nrow(x)==0}))
list_jointatt<-list_jointatt[-dropped_participants]



hist(sapply(list_jointatt,nrow))
hist(sapply(list_jointatt,function(x){length(unique(x$index_trial))}))

paste('date:',paste(Sys.Date()),',data of participants: n=',length(list_jointatt))

## ----------> GAZE PREPROCESSING (Nyström, 2010) ---------------------- ####

start_time <- Sys.time()

#variables for gaze preprocessing
screen_width<-510 #mm on a 23 inch screen wiht 16:9 aspect ratio (Tobii TX 300 screen)
screen_height<-290 #mm on a 23 inch screen wiht 16:9 aspect ratio (Tobii TX 300 screen)
degrees_by_radian<-180/pi #fixed conversion facor
velocity_cutoff<-1000 #visual degress per second
acceleration_cutoff<-100000 #visual degress per second
initial_velocity_cutoff<-200
median_samples_trial<-3300

#--> Savitksy Golay filter of length 15 ~ 50ms --> see for coefficients: http://www.statistics4u.info/fundstat_eng/cc_savgol_coeff.html
#filter_sg15<-c(-78,-13,42,87,122,147,162,167,162,147,122,87,42,-13,-78)/1105
filter_sg21<-c(-171,-76,9,84,149,204,249,284,309,324,329,324,309,284,249,204,149,84,9,-76,-171)/3059

#required functions for gaze preprocessing:

#A. - blink identification function - within 75ms to 250ms interval
#--> INFO: identifies blinks and NAs 8 samples before and after it (~25ms)
fun_blink_cor <- function(signal,lower_threshold=23,upper_threshold=75,samples_before=8,samples_after=8) {
  #change NA to 999 for rle()-function
  findna <- ifelse(is.na(signal),999,signal)
  #find blinks:
  #output of rle(): how many times values (NA) are repeated
  repets <- rle(findna)
  #stretch to length of PD vector for indexing
  repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
  #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
  #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
  #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
  repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
  #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
  repets[findna!=999 & repets==1] <- 0
  #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
  changes <- c(diff(repets),0)
  #define start (interval before blink/missing data)
  changes.start<-which(changes==1) #where NA-sequence starts
  #gives out row numbers of NA (blink) and previous 8 frames
  start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
  repets[start.seq]<-1
  #define end (interval after blink/missing data)
  changes.end<-which(changes==-1)+1 #where NA.sequence ends
  #gives out row numbers of NA (blink) and subsequent 8 frames
  end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
  repets[end.seq]<-1
  #replace PD data in blink interval (start to end) with NA
  signal[repets==1]<-NA
  return(signal)
}

#B. - return velocity
fun_return_speed<-function(a,b,time){
  time<-unlist(time)
  gaze_diff_x<-diff(a)
  gaze_diff_y<-diff(b)
  gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
  gaze_speed<-gaze_diff/diff(time)
  gaze_speed<-c(NA,gaze_speed)
  return(gaze_speed)
}

#C. - return acceleration
fun_return_accel<-function(x,time){
  time<-unlist(time)
  diff_speed<-diff(x)
  gaze_accel<-diff_speed/diff(time)
  gaze_accel<-c(NA,gaze_accel)
  return(gaze_accel)
}

#D. - function data driven velocity threshold to identify saccades --> see Nyström et al., 2010
fun_velocity_threshold<-function(x){

  cutoff<-initial_velocity_cutoff
  diff_cutoff<-cutoff
  k<-x

  while(diff_cutoff>1){

    k<-k[k<cutoff]
    mean_speed<-median(k,na.rm=T)
    sd_speed<-mean(k,na.rm=T)

    cutoff_new<-mean_speed+(3*sd_speed)
    if(is.na(cutoff_new)){
      cutoff<-NA
      break}

    diff_cutoff<-cutoff-cutoff_new
    if(cutoff_new>cutoff){cutoff_new<-cutoff}
    cutoff<-cutoff_new

  }

  # #testing
  # return(cutoff)

  saccade_peak<-ifelse(x>cutoff,T,F)
  return(saccade_peak)
  #also return cutoffs in a list

}

#idnetififes most often values in a sequence
fun_most_values <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  return(ux[tab == max(tab)])

}

#function gaze preprocessing (requires functions above)
#--> inspired by Nyström et al. 2010 - denoising with Savitsky-Golay filter and adaptive velocity threshold filter
fun_gaze_preprocess<-function(x){

  #testing
  #x<-list_jointatt[[2]]

  #preprocess based on per trial level
  x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
  split_by_trial<-split(x,as.factor(x$index_trial))

  #split_by_trial<-split(x,as.factor(x$index_trial))
  timestamp<-sapply(split_by_trial,function(x){x<-x['ts']})
  screen_distance<-sapply(split_by_trial,function(x){x<-x['screen_dist']})

  #1. blink correction --> as noise
  gazepos_x<-sapply(split_by_trial,function(x){fun_blink_cor(x$gazepos.x)})
  gazepos_y<-sapply(split_by_trial,function(x){fun_blink_cor(x$gazepos.y)})

  #3. drop trials with less than 50% of data
  gazepos_x<-sapply(gazepos_x,function(x){if(sum(is.na(x))>0.5*median_samples_trial){x<-as.numeric(rep(NA,length(x)))}else{return(x)}})
  gazepos_y<-sapply(gazepos_y,function(x){if(sum(is.na(x))>0.5*median_samples_trial){x<-as.numeric(rep(NA,length(x)))}else{return(x)}})

  #2. convert relative gaze to degrees visual angle --> degrees from point of origin
  gazepos_x_deg<-mapply(function(x,y){x<-x*atan(screen_width/y)*degrees_by_radian},x=gazepos_x,y=screen_distance,SIMPLIFY = F)
  gazepos_y_deg<-mapply(function(x,y){x<-x*atan(screen_height/y)*degrees_by_radian},x=gazepos_y,y=screen_distance,SIMPLIFY = F)

  #put together
  unsplitting_factor<-as.factor(x$index_trial)
  unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
  gazepos_x<-unsplit(gazepos_x,f=unsplitting_factor)
  gazepos_y<-unsplit(gazepos_y,f=unsplitting_factor)
  gazepos_x_deg<-unsplit(gazepos_x_deg,f=unsplitting_factor)
  gazepos_y_deg<-unsplit(gazepos_y_deg,f=unsplitting_factor)

  x<-data.frame(unsplit_by_trials,gazepos_x,gazepos_y,gazepos_x_deg,gazepos_y_deg)
  x$index_trial[x$index_trial==0]<-NA
  return(x)

}

list_jointatt<-pblapply(list_jointatt,fun_gaze_preprocess)

#saccade identification based on velocity
fun_saccade_ident<-function(x){

  #x<-fun_gaze_preprocess(data_block1[[100]])

  #preprocess based on per trial level
  x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
  split_by_trial<-split(x,as.factor(x$index_trial))
  timestamp<-sapply(split_by_trial,function(x){x<-x['ts']})
  gazepos_x_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_x_deg']})
  gazepos_y_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_y_deg']})


  # A. VELOCITY: FILTERING, DENOISING, PEAK IDENTIFICATION

  #3. return velocity and acceleration variable --> in degrees per second
  gaze_speed<-mapply(fun_return_speed,a=gazepos_x_deg,b=gazepos_y_deg,time=timestamp,SIMPLIFY = F)
  gaze_accel<-mapply(fun_return_accel,x=gaze_speed,time=timestamp,SIMPLIFY = F)

  #4. remove biologically implausible acceleration and velocity values
  gaze_speed<-sapply(gaze_speed,function(x){ifelse(abs(x)>velocity_cutoff,NA,x)})
  gaze_accel<-sapply(gaze_accel,function(x){ifelse(abs(x)>acceleration_cutoff,NA,x)})

  #5. Denoising by Savitzky-Golay filter (length 21)
  gaze_speed<-lapply(gaze_speed,function(x){
    ifelse(length(x)<length(filter_sg21),
           k<-x,
           k<-as.numeric(stats::filter(x,filter_sg21)))
    return(k)})

  #--> SOME FORM OF DENOISING HAS TO BE APPLIED
  #literature of signal denoise
  #use a Savitzky Golay filter with length 15 ~ 50ms --> has been proven to be superior, see below


  #filter/denoise benchmarking #
  # sg7<-c(-2,3,6,7,6,3,-2)/21 #Savitsky Golay filter of length 7
  # sg5<-c(-3,12,17,12,-3)/35 #Savitsky Golay filter of length 5
  # sg21<-c(-171,-76,9,84,149,204,249,284,309,324,329,324,309,284,249,204,149,84,9,-76,-171)/3059
  # sg15<-c(-78,-13,42,87,122,147,162,167,162,147,122,87,42,-13,-78)/1105
  # ##--> see for coefficients: http://www.statistics4u.info/fundstat_eng/cc_savgol_coeff.html
  #
  # x<-gaze_speed[[15]]
  # par(mfrow=c(3,2))
  # plot(x,main='data')
  # plot(as.numeric(filter(x,sg5)),main='S-G5 filter')
  # plot(as.numeric(filter(x,sg15)),main='S-G15 filter')
  # plot(as.numeric(filter(x,sg21)),main='S-G21 filter')
  # #plot((na.approx(x,maxgap = interpolate_cutoff)),main='linear interpol')
  # plot(pracma::savgol(na.approx(x),fl=21),main='S-G long')
  # plot(spline(x,n=length(x)),main='spline')
  # par(mfrow=c(1,1))
  # ###--> compare different denoise / filtering methods
  # ##--> Savitzky-Golay filter seems to be the winner - use a length of 21

  #6.data-driven velocity threshold and identify peaks (saccades) based on it
  velocity_peak<-lapply(gaze_speed,fun_velocity_threshold)

  # # #  # #testing
  # plot(gaze_speed[[10]],col=as.factor(unlist(velocity_peak[[10]])))
  # plot(unlist(gaze_speed),col=as.factor(unlist(velocity_peak)))

  #B.: DISPERSION

  # fun_dispersion_filter<-function(a,b){
  #
  #   gaze_diff_x<-diff(a)
  #   gaze_diff_y<-diff(b)
  #   gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
  #   gaze_diff<-c(NA,gaze_diff)
  # }
  #
  # gaze_diff<-mapply(fun_dispersion_filter,a=gazepos_x,b=gazepos_y,SIMPLIFY=F)

  #C.: DURATION

  #modus in a window of ten samples
  velocity_continues<-frollapply(velocity_peak,n=10,fun_most_values,align='center')
  velocity_continues<-lapply(velocity_continues,as.logical)

  #D.: DEFINE SACCADE
  saccade<-mapply(function(x,y){ifelse(x & y,T,F)},x=velocity_peak,y=velocity_continues,SIMPLIFY=F)

  #return trial-wise list to data.frame
  unsplitting_factor<-as.factor(x$index_trial)
  unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
  gaze_speed<-unsplit(gaze_speed,f=unsplitting_factor)
  gaze_accel<-unsplit(gaze_accel,f=unsplitting_factor)
  velocity_peak<-unsplit(velocity_peak,f=unsplitting_factor)
  velocity_continues<-unsplit(velocity_continues,f=unsplitting_factor)
  saccade<-unsplit(saccade,f=unsplitting_factor)


  x<-data.frame(unsplit_by_trials,gaze_speed,gaze_accel,velocity_peak,velocity_continues,saccade)
  x$index_trial[x$index_trial==0]<-NA
  return(x)
}

list_jointatt<-pblapply(list_jointatt,fun_saccade_ident)

#fixation identification
fun_fixation_ident<-function(x,degree_fixation_cutoff=1,duration_fixation_cutoff=30){

  #testing
  # x<-data_block1[[100]]
  # x<-fun_gaze_preprocess(x)
  # x<-fun_saccade_ident(x)

  x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
  split_by_trial<-split(x,as.factor(x$index_trial))
  gazepos_x_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_x_deg']})
  gazepos_y_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_y_deg']})


  #gaze difference from sampel to sample
  fun_gaze_diff<-function(a,b){

    gaze_diff_x<-diff(a)
    gaze_diff_y<-diff(b)
    gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
    gaze_diff<-c(NA,gaze_diff)

  }

  gaze_diff<-mapply(fun_gaze_diff,a=gazepos_x_deg,b=gazepos_y_deg,SIMPLIFY=FALSE)

  #identify significant movement in subsequent or preceding samples

  fun_movement_ident<-function(x){

    no_movement_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(mean(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='left')
    no_movement_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(mean(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='right')
    #no_movement_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(all(x<degree_fixation_cutoff),T,F)},align='left')
    #no_movement_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(all(x<degree_fixation_cutoff),T,F)},align='right')
    no_movement<-ifelse(no_movement_next_samples | no_movement_last_samples,T,F)

  }

  no_movement<-sapply(gaze_diff,fun_movement_ident)

  #identify significant drifts in subsequent or preceding samples

  fun_gaze_diff_abs<-function(x){

    gaze_diff_abs<-diff(x)
    gaze_diff_abs<-c(NA,gaze_diff_abs)

  }

  gaze_diff_abs_x<-sapply(gazepos_x_deg,fun_gaze_diff_abs)
  gaze_diff_abs_y<-sapply(gazepos_y_deg,fun_gaze_diff_abs)

  fun_drift_ident<-function(x){

    no_drift_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(sum(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='left')
    no_drift_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(sum(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='right')
    no_drift<-ifelse(no_drift_next_samples | no_drift_last_samples,T,F)

  }

  no_drift_x<-sapply(gaze_diff_abs_x,fun_drift_ident)
  no_drift_y<-sapply(gaze_diff_abs_y,fun_drift_ident)


  #--> define a fixation sample if no movement or drift in the last or next 100ms ~ 30 samples
  fixation<-mapply(function(x,y,z){ifelse(x & y & z,T,F)},x=no_movement,y=no_drift_x,z=no_drift_y)

  #add to data
  unsplitting_factor<-as.factor(x$index_trial)
  unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
  fixation<-unsplit(fixation,f=unsplitting_factor)
  x<-data.frame(unsplit_by_trials,fixation)
  x$index_trial[x$index_trial==0]<-NA
  return(x)

}

list_jointatt<-pblapply(list_jointatt,fun_fixation_ident)

## ----------> PUPIL Preprocessing (Kret, 2018) -------------- ####

fun_blink_cor <- function(signal,lower_threshold=23,upper_threshold=75,samples_before=8,samples_after=8) {
  #change NA to 999 for rle()-function
  findna <- ifelse(is.na(signal),999,signal)
  #find blinks:
  #output of rle(): how many times values (NA) are repeated
  repets <- rle(findna)
  #stretch to length of PD vector for indexing
  repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
  #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
  #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
  #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
  repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
  #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
  repets[findna!=999 & repets==1] <- 0
  #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
  changes <- c(diff(repets),0)
  #define start (interval before blink/missing data)
  changes.start<-which(changes==1) #where NA-sequence starts
  #gives out row numbers of NA (blink) and previous 8 frames
  start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
  repets[start.seq]<-1
  #define end (interval after blink/missing data)
  changes.end<-which(changes==-1)+1 #where NA.sequence ends
  #gives out row numbers of NA (blink) and subsequent 8 frames
  end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
  repets[end.seq]<-1
  #replace PD data in blink interval (start to end) with NA
  signal[repets==1]<-NA
  return(signal)
}

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
  #gaps=missing data sections > 75ms; Leonie: also <=250ms, otherwise not likely to be a blink
  #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
  pl<-fun_blink_cor(pl)
  pr<-fun_blink_cor(pr)

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

  #interpolation of NA (for <=300ms)
  pl<-na.approx(pl, na.rm=F, maxgap=90, rule=2)
  pr<-na.approx(pr, na.rm=F, maxgap=90, rule=2)

  pd <- (pl+pr)/2
  # end of function --> return ####
  #detach(x)

  x[,'pd']<-pd
  return(x)
}

# pd.list<-lapply(df, func_pd_preprocess)
list_jointatt<-pblapply(list_jointatt, func_pd_preprocess)

### - select relevant data ####

fun_select_data<-function(x){
  x<-x[,names(x) %in% c('eventlog','timestamp','ts','index_trial',
                        'ts_event','gazepos.x','gazepos.y',
                        'pd','center_deviation','screen_dist','gazepos_x','gazepos_y','gaze_speed','gaze_accel',
                        'velocity_peak','velocity_continues','saccade','fixation')]
  return(x)
}

list_jointatt<-pblapply(list_jointatt,fun_select_data)

### - create ID variable - required before concatenating blocks ####
fun_retrieve_id<-function(x,id){
  k<-nrow(x)
  id<-rep(id,k)
  x[,'id']<-id
  return(x)
}

list_jointatt<-mapply(fun_retrieve_id,x=list_jointatt,id=names(list_jointatt),SIMPLIFY=F)

### - calculate baseline pupil size - corrected pupil size ####
fun_baseline<-function(x){

  split_by_trial<-split(x,as.factor(x$index_trial))
  baseline_data<-sapply(split_by_trial,function(x){mean(x$pd[x$ts_event<=150],na.rm=T)}) #select between event
  #rpd<-unlist(mapply(function(x,y){x$pd-y},x=split_by_trial,y=baseline_data)) #correct for baseline by subtraction
  rpd<-unsplit(mapply(function(x,y){x$pd-y},x=split_by_trial,y=baseline_data),f=as.factor(x$index_trial)) #correct for baseline by subtraction
  baseline_pd<-rep(baseline_data,sapply(split_by_trial,nrow)) #calculate baseline pd

  x[,'rpd']<-rpd
  x[,'baseline_pd']<-baseline_pd
  return(x)

}

list_jointatt<-pblapply(list_jointatt,fun_baseline)

###------ save preprocessed data (single sample level) #####
df_jointatt<-dplyr::bind_rows(list_jointatt)

save(df_jointatt,file=paste0(home_path,project_path,"/data/all_data_preprocessed_jointatt_130324.Rdata"))

#### limit data to first 11 seconds ####

#load(paste0(home_path,project_path,"/data/all_data_preprocessed_jointatt_060922.Rdata"))
#load(paste0(home_path,project_path,"/data/all_data_preprocessed_jointatt_150922.Rdata"))

hist(df_jointatt$ts_event,40)
df_jointatt<-df_jointatt[df_jointatt$ts_event<=3300,]

### - CREATE PIC, GROUP, TIMEPOINT, CONDITION, STIMULUS, POSITION, VARIABLE ####

#create group variable
df_jointatt$group<-ifelse(grepl('_K',df_jointatt$id),'TD','ASD')

#create timepoint variable
df_jointatt$timepoint<-ifelse(grepl('_t2',df_jointatt$id) | grepl('_T2',df_jointatt$id),'T2',
                             ifelse(grepl('_t4',df_jointatt$id) | grepl('_T4',df_jointatt$id),'T4',
                                    ifelse(grepl('_t6',df_jointatt$id) | grepl('_T6',df_jointatt$id),'T6',
                                           ifelse(grepl('_fu2',df_jointatt$id) | grepl('_FU2',df_jointatt$id),'FU2',
                                                         ifelse(grepl('_fu3',df_jointatt$id) | grepl('_FU3',df_jointatt$id),'FU3',
                                                                       ifelse(grepl('K_fu',df_jointatt$id) | grepl('K_FU',df_jointatt$id) | grepl('k_FU',df_jointatt$id) | grepl('k_fu',df_jointatt$id),'K_FU','K'))))))

table(df_jointatt$timepoint)

#create individual id variable
df_jointatt$pic<-substr(df_jointatt$id,1,3)

#create CONDITION
condition<-substr(df_jointatt$eventlog,16,20)
condition<-ifelse(condition=='inten','intense',
                  ifelse(condition=='mild_','mild',
                         ifelse(condition=='neutr','neutral','point')))
df_jointatt$condition<-condition
rm(condition)

#create STIMULUS
label_length<-nchar(as.character(df_jointatt$eventlog))
stimulus<-substr(df_jointatt$eventlog,label_length-16,label_length-9)
stimulus<-ifelse(grepl('flower',stimulus),'flower',
                 ifelse(grepl('ball',stimulus),'ball',
                        ifelse(grepl('truck',stimulus),'truck',
                               ifelse(grepl('rab',stimulus),'rabbit',NA))))
df_jointatt$stimulus<-stimulus
rm(stimulus)

#create POSITION
position<-substr(df_jointatt$eventlog,label_length-8,label_length-4)
table(position)
df_jointatt$position<-ifelse(position=='right','target_right','target_left')
rm(position)

### - CREATE GAZE BEHAVIOR VARIABLES ####

#gaze_fixcross, gaze_face, gaze_stimulus_left, gaze_stimulus_right

fun_ident_hits<-function(df,position,time_window){

  xstart<-position[1]
  xend<-position[2]
  ystart<-position[3]
  yend<-position[4]

  time_start<-time_window[1]
  time_end<-time_window[2]

  gaze_x<-df$gazepos_x
  gaze_y<-df$gazepos_y
  time<-df$ts_event

  #identify hits - compare gaze position to target position + matching time
  hit<-rep(NA,nrow(df))
  hit<-ifelse((gaze_x > xstart) &
                (gaze_x < xend) &
                (gaze_y > ystart) &
                (gaze_y < yend),T,F)

  hit<-ifelse((time >= time_start) & (time <= time_end),hit,NA)

  return(hit)

}

fixation_cross_position<-c(0.4,0.6,0.025,0.325) # in relative position - 0.2 * 0.3
fixation_cross_timing<-c(1,450) #in frames (first 1.5 seconds)
df_jointatt$gaze_fixcross<-fun_ident_hits(df=df_jointatt,
                                          position=fixation_cross_position,
                                          time_window=fixation_cross_timing)

head_position<-c(0.4,0.6,0.025,0.325) # in relative position - 0.2 * 0.3
head_timing<-c(451,3300) #in frames (last 9 seconds)
df_jointatt$gaze_head<-fun_ident_hits(df=df_jointatt,
                                          position=head_position,
                                          time_window=head_timing)

head_position<-c(0.4,0.6,0.025,0.325) # in relative position - 0.2 * 0.3
head_timing<-c(0,3300) #in frames (last 9 seconds)
df_jointatt$gaze_top<-fun_ident_hits(df=df_jointatt,
                                      position=head_position,
                                      time_window=head_timing)

head_before_position<-c(0.4,0.6,0.025,0.325) # in relative position - 0.2 * 0.3
head_before_timing<-c(451,1050) #in frames (1.66 seconds)
df_jointatt$gaze_head_before<-fun_ident_hits(df=df_jointatt,
                                      position=head_before_position,
                                      time_window=head_before_timing)

head_after_position<-c(0.4,0.6,0.025,0.325) # in relative position - 0.2 * 0.3
head_after_timing<-c(2700,3300) #in frames (2 seconds)
df_jointatt$gaze_head_after<-fun_ident_hits(df=df_jointatt,
                                      position=head_after_position,
                                      time_window=head_after_timing)

stimulus_left_position<-c(0.175,0.375,0.75,1) # in relative position - 0.2 * 0.25
stimulus_left_timing<-c(450,3300) #in frames
df_jointatt$gaze_stimulus_left<-fun_ident_hits(df=df_jointatt,
                                            position=stimulus_left_position,
                                            time_window=stimulus_left_timing)

stimulus_right_position<-c(0.625,0.825,0.75,1) # in relative position - 0.2 * 0.25
stimulus_right_timing<-c(450,3300) #in frames
df_jointatt$gaze_stimulus_right<-fun_ident_hits(df=df_jointatt,
                                            position=stimulus_right_position,
                                            time_window=stimulus_right_timing)


df_jointatt$gaze_stimulus<-df_jointatt$gaze_stimulus_right|df_jointatt$gaze_stimulus_left
df_jointatt$gaze_target<-(df_jointatt$gaze_stimulus_left & df_jointatt$position=='target_left') | (df_jointatt$gaze_stimulus_right & df_jointatt$position=='target_right')
df_jointatt$gaze_distractor<-(df_jointatt$gaze_stimulus_left & df_jointatt$position=='target_right') | (df_jointatt$gaze_stimulus_right & df_jointatt$position=='target_left')

df_jointatt$gazehead_when<-ifelse(df_jointatt$gaze_head,df_jointatt$ts_event,NA)
df_jointatt$gazetop_when<-ifelse(df_jointatt$gaze_top,df_jointatt$ts_event,NA)

df_jointatt$gazetarget_when<-ifelse(df_jointatt$gaze_target,df_jointatt$ts_event,NA)
df_jointatt$gazedistractor_when<-ifelse(df_jointatt$gaze_distractor,df_jointatt$ts_event,NA)


      # #measure diagnostics
      #   #target
      # ggplot(df_jointatt[df_jointatt$gaze_target,],aes(x=gazepos_x,y=gazepos_y))+
      #   geom_hex(bins=30)+scale_fill_gradientn(colours=rev(rainbow(3)))+
      #   ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+theme_bw()
      #
      #   #stimulus
      # ggplot(df_jointatt[df_jointatt$gaze_head,],aes(x=gazepos_x,y=gazepos_y))+
      #   geom_hex(bins=30)+scale_fill_gradientn(colours=rev(rainbow(3)))+
      #   ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+theme_bw()
      #
      #
      # ggplot(df_jointatt,aes(x=gazepos_x,y=gazepos_y))+
      #   geom_hex(bins=30)+scale_fill_gradientn(colours=rev(rainbow(3)))+
      #   ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+theme_bw()


## MAIN OUTCOME: DEFINE REACTIVE JOINT ATTENTION (RJA) ####

#split by trial per participant
split_by_trial<-split(df_jointatt,interaction(df_jointatt$id,df_jointatt$index_trial))

# Franchini 2016 - RJA definition
# was calculated as the number of correct gaze shifts between the face (or the hand in the pointing condition) and the AOI of
# the object the actor was looking at (one point per condition when present). Also, we calculated RJA errors as the number of
# gaze shifts between the face/hand and the AOI of the unreferenced object. A gaze shift was considered successful when the
# child’s gaze went directly from the face/hand to the object within 500 ms at least one time per video.

#RJA function - considers gazes on head/fixcross in last 500ms and gaze on target now - delivers mesaure of reactive joint attention
fun_find_rja<-function(x,samples_of_500ms=150,rja_start=900){

  # #testing
  # x<-split_by_trial[[1387]]
  # samples_of_500ms<-150
  # rja_start<-900

  #was gaze on head/fixcross in last 500ms?
  #last 500 ms or smaller range of sample if at the beginning of the vector (i.e. sample<150)
  samples_to_consider<-sapply(x$ts_event,function(y){
    first_sample<-ifelse((y-samples_of_500ms)<1,1,y-samples_of_500ms)
    last_sample<-y
    seq(first_sample,last_sample)
  })
  gazes_to_consider<-sapply(samples_to_consider,function(y){x$gaze_head[y]})
  #ALSO consider that these gazes were fixations
  gazes_with_fixation<-sapply(samples_to_consider,function(y){x$fixation[y]})

  #any fixation in last 500ms on head/fixcross
  fix_head_last500ms<-mapply(function(x,y){ifelse(any((x & y),na.rm=T),T,F)},x=gazes_to_consider,y=gazes_with_fixation)
  fix_head_last500ms<-as.logical(fix_head_last500ms) #set to logical - for empty trials
  #DEFINE RJA
  #- fixate target (gaze left if target left or gaze right if target right) & fixation_occured & fixated head/fixcross in last 500ms & rja started (when objects start to shake)
  #- return NA if any of parameters is na
  #rja <- ifelse(x$gaze_target & x$fixation & fix_head_last500ms & (x$ts_event > rja_start),T,F)
  rja <- ifelse(is.na(x$gaze_target)|is.na(x$fixation)|is.na(fix_head_last500ms),NA,
                ifelse(x$gaze_target & x$fixation & fix_head_last500ms & (x$ts_event > rja_start),T,F))


  # #DIAGNOSTICS
  # table(gaze_target)
  # table(x$fixation)
  # table(fix_head_last500ms)
  # table(rja)

  # also tja if last sample was rja (in case target is still fixated while fixation of head was longer than 500 ms ago)
  rja_last_sample<-c(F,rja[-(length(rja)-1)])
  rja <- ifelse(rja | (rja_last_sample & x$gaze_target),T,F)

  return(rja)

}


df_jointatt$rja<-unsplit(pbsapply(split_by_trial,fun_find_rja),
                         f=interaction(df_jointatt$id,df_jointatt$index_trial))

#additional variables
df_jointatt$rja_when<-ifelse(df_jointatt$rja,df_jointatt$ts_event,NA)
df_jointatt$rja_early<-ifelse(df_jointatt$rja_when<1200,T,F)
df_jointatt$rja_late<-ifelse(df_jointatt$rja_when>1200 & df_jointatt$rja_when<2100,T,F)


    # # # #measure diagnostics
    # hist(df_jointatt$rja_when)
    # with(df_jointatt,table(rja,condition,group))
    # with(df_jointatt,by(rja,interaction(group,condition),function(x){table(x)[2]/table(x)[1]}))
    # with(df_jointatt,by(rja,group,sum,na.rm=T))
    # with(df_jointatt,table(rja,group))
    #  with(df_jointatt,table(rja,timepoint)
    # #-->descriptive higher rja rate in TD

    # ggplot(df_jointatt,aes(x=rja_when,group=group,fill=group))+geom_histogram()+facet_wrap(~condition+group)+theme_bw()
    # ggplot(df_jointatt,aes(x=rja_when,group=group,fill=group))+geom_density()+facet_grid(vars(condition),vars(group))+theme_bw()
    # ggplot(df_jointatt,aes(x=rja_when,group=group,fill=group))+geom_density(alpha=0.5)+theme_bw()
    # ggplot(df_jointatt,aes(x=rja_when,group=group,fill=group))+geom_density()+theme_bw()
    #
    # summary(df_jointatt$rja)
    # with(df_jointatt,table(index_trial))
    # with(df_jointatt,table(droplevels(eventlog)))
    # with(df_jointatt,hist(gazepos_x))
    # with(df_jointatt,hist(gazepos_y))
    # with(df_jointatt,hist(rpd))
    # with(df_jointatt,table(rja))


#- DEMOGRAPHIC DATA #####

#demographic data is created with sample_overview_Jan2024.Rmd
#TODO: need new data export
load(paste0(home_path,project_path,"/data/demogr_data_270722"))
load(paste0(home_path,"/PowerFolders/data_AFFIP/demogr_total_baseline_0623.Rdata"))

    unique(df_jointatt$pic)
    demogr$id
    ###->
    unique(df_jointatt$pic)[!(unique(df_jointatt$pic) %in% substr(demogr$id,1,3))]
    ### one has ET data but no demographics (14.03.2024)
    hist(table(df_jointatt$id))

unique(demogr$id)
demogr$id<-substr(demogr$id,1,3)

demogr$group<-ifelse(substr(demogr$id,1,1)=='9','TD','ASD')

    ##-->currently FU is missing
    table(demogr$group,demogr$t_IQ)
    table(demogr$group,demogr$t)


    ## - exclude IDs  with assessments problems (n=4)####

    # Es gibt einige IDs, die ich ausgeschlossen hatte, zu denen es aber Eye-Tracking-Daten gibt:
    #   -036 und 085 hatten die ADI-R/ADOS-Cutoffs nicht erreicht
    # -925_K war sehr auffällig Richtung ADHS, Eye-Tracking hatte kaum geklappt
    # -942_K hat mittlerweile eine ASS-Diagnose, nimmt jetzt als 096 an A-FFIP-Studie teil; 941_K ist das zugehörige Geschwisterkind (laut Mutter mittlerweile auch auffällig)
    # -965_K war sehr auffällig in der Testung, FSK war auch auffällig
    #
    # Zusätzlich gab es beim Speichern Probleme mit Kind 076 und 108: Für diese sind aus irgendeinem Grund in den
    # Datensätzen jeweils Kinder, die an anderen Tagen erhoben wurden, mit abgespeichert.
    # Anbei Code zum Rausschmeißen der falschen Datenpunkte:
    #

    exclude_ids<-c('036','085','925','942','965')
    demogr<-demogr[!(demogr$id %in% exclude_ids),]

    ## - exclude IDs without eyetracking data (n=9) ####
    demogr<-demogr[(demogr$id %in% unique(df_jointatt$pic)),]

    ## - compare sample ####
    with(demogr,by(Geschlecht_Index,group,table))
    with(demogr,by(test_age,group,summary))
    with(demogr,by(IQ,group,summary))

    with(demogr,by(test_age,group,summary))
    sort(demogr$test_age[demogr$group=='ASD'])
    sort(demogr$test_age[demogr$group=='TD'])

    ## - remove IDs with test_age higher than max of ASD (helps in matching) ####
    #demogr<-demogr[demogr$test_age < max(demogr$test_age[demogr$group=='ASD']),]


    ###no matching required in final sample, when treatment versus control is unblinded

    ### --> match sample ####
    require(MatchIt)

    groupBoo<-with(demogr,ifelse(group=='ASD',1,0))

    #ALL - matching
    set.seed(100)
    all.match<-matchit(groupBoo~test_age,
                       data=demogr,
                       method='nearest',discard='both',
                       ratio=8, #match four controls to each ASD
                       replace=T,caliper=0.20)

    all.match

    summary(all.match)

    #remove unmatched cases - from eye tracking data
    all.match<-match.data(all.match)

    with(all.match,t.test(test_age~group))

    #df_jointatt<-df_jointatt[df_jointatt$pic %in% all.match$id,]

#--> DATA PLAUSIBILITY AFTER MATCHING ####

### - check data distribution  ####


#trial time
hist(df_jointatt$ts_event)

#initial data plausibility --> similar data scross trials
table(df_jointatt$index_trial)

#individual data sets
  #length(unique(df_jointatt$id)) #data of 171 measurement timepoints after matching
  #length(unique(df_jointatt$pic)) #data of 101 participants after matching)

length(unique(df_jointatt$id)) #data of 300 measurement timepoints before matching (march 2024)
length(unique(df_jointatt$pic)) #data of 143 participants before matching (march 2024)
by(df_jointatt$id,df_jointatt$timepoint,function(x){length(unique(x))})
##--> after matching (September 2022): K = 44, T2 = 51, T4 = 44, T6 = 32
##--> without matching (March 2023): K = 65, T2 = 62, T4 = 53, T6 = 47, K_FU = 20, FU2 = 15, FU3 = 2
##--> without matching (March 2024): K = 75, T2 = 62, T4 = 53, T6 = 47, K_FU = 26, FU2 = 25, FU3 = 12

#missing data
table(is.na(df_jointatt$rpd))[2]/sum(table(is.na(df_jointatt$rpd))) #--> 43.2 missing data in pupil data
table(is.na(df_jointatt$gazepos_x))[2]/sum(table(is.na(df_jointatt$gazepos_x))) #--> 51.8% missing data in gaze x
table(is.na(df_jointatt$gazepos_y))[2]/sum(table(is.na(df_jointatt$gazepos_y))) #--> 56.0% missing data in gaze y

table(df_jointatt$saccade) #roughly 10% saccade data is plausible
table(df_jointatt$fixation) #roughly 90% fixation data is plausible

#data distribution --> exclude implausible
par(mfrow = c(3, 2))
hist(df_jointatt$pd,30,main='pupil size (mm)')
#hist(df_jointatt$baseline_pd,30,main='baseline pupil size (mm)') #pupil size in first 500ms of a trial
hist(df_jointatt$rpd[df_jointatt$rpd<2 & df_jointatt$rpd>-2],30,main='pupillary response') #--> implausible outlier #pupil response (baseline pd - pupil size)
hist(df_jointatt$gaze_speed,30,main='gaze velocity (degrees/s)')
hist(df_jointatt$gaze_accel,30,main='gaze acceleration (degrees/s²)')
hist(df_jointatt$gazepos_y,30,main='gaze location on x-axis (0-1)')
hist(df_jointatt$gazepos_x,30,main='gaze location on y-axis (0-1)')
par(mfrow = c(1, 1))

#VALID TRIALS by participants
hist(as.numeric(with(df_jointatt[df_jointatt$fixation,],
                     by(index_trial,id,function(x){length(unique(x))}))),xlab='trials',main='trials with fixation data by participant')

### --> visualize gaze behavior (for INSAR 2023 abstract) ####

#gaze behavior
require(ggplot2)
require(hexbin)
# ggplot(df_jointatt,aes(x=gazepos_x,y=gazepos_y))+
#   geom_hex(bins=30)+
#   scale_fill_gradientn(colours=rev(rainbow(3)))+
#   ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)

ggplot(df_jointatt,aes(x=gazepos_x,y=gazepos_y))+
  stat_density_2d(aes(fill=after_stat(density)), n=50, geom = "raster", contour = FALSE)+
  scale_fill_gradientn(colours=rev(rainbow(3)))


#before rja cueing onset
ggplot(df_jointatt[df_jointatt$ts_event<1000,],aes(x=gazepos_x,y=gazepos_y))+
  #geom_hex(bins=30)+
  stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours=rev(rainbow(3)))+
  ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
  #facet_wrap(~timepoint)+
  theme_void()

#after RJA cueing onset
ggplot(df_jointatt[df_jointatt$ts_event>1300 & df_jointatt$ts_event<1600,],aes(x=gazepos_x,y=gazepos_y))+
  #geom_hex(bins=30)+
  stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours=rev(rainbow(3)))+
  ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
  #facet_wrap(~timepoint)+
  theme_void()

#late RJA
ggplot(df_jointatt[df_jointatt$ts_event>2000 & df_jointatt$ts_event<3000,],aes(x=gazepos_x,y=gazepos_y))+
  #geom_hex(bins=30)+
  stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours=rev(rainbow(3)))+
  ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
  facet_wrap(~timepoint)+
  theme_bw()

#late RJA
ggplot(df_jointatt[df_jointatt$ts_event>2700 & df_jointatt$ts_event<3300,],aes(x=gazepos_x,y=gazepos_y))+
  #geom_hex(bins=30)+
  stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours=rev(rainbow(3)))+
  ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
  #facet_wrap(~timepoint)+
  theme_void()



ggplot(df_jointatt[df_jointatt$ts_event>1000 & df_jointatt$ts_event<1100,],aes(x=gazepos_x,y=gazepos_y))+
  geom_bin_2d()+
  scale_fill_gradientn(colours=rev(rainbow(3)))+
  ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
  facet_wrap(~timepoint)+
  theme_bw()


      #TODO: animate gaze behavior
      require(gganimate)
      require(gifski) #gif renderer - animate()
      #transition variable in animation:
      frame_rate<-300
      df_jointatt$time_s<-with(df_jointatt,round(ts_event/frame_rate,1))



      #BETWEEN GROUPS
      animated_gaze<-ggplot(df_jointatt,aes(x=gazepos_x,y=gazepos_y))+
        stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
        #geom_bin_2d(bins=100)+
        #scale_fill_gradientn(colours=rev(rainbow(3)))+
        scale_fill_gradientn(colours = c("black","white","red","orange","yellow","green","blue"),
                               values = c(1.0,0.5,0.2,0.15,0.1,0.05,0))+
        ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
        facet_wrap(~timepoint)+
        theme_bw()+
        theme(legend.position="none")+
        #animation specific
        #transition_time(time_s)+facet_wrap(~timepoint)
        transition_time(time_s)+
        ggtitle('Gaze behavior between groups','time: {frame_time}')

      gif_1<-animate(animated_gaze, duration = 10, fps = 10, width = 800, height = 800, renderer = gifski_renderer())
      anim_save(file=paste0(home_path,project_path,'/output/gaze_animation_time_2ddensity.gif'),gif_1)

      hist(df_jointatt$time_s)
      #ACROSS GROUPS
      animated_gaze<-ggplot(df_jointatt[df_jointatt$time_s>9 & df_jointatt$time_s<11,],aes(x=gazepos_x,y=gazepos_y))+
        stat_density_2d(aes(fill = ..density..), n=50, geom = "raster", contour = FALSE) +
        #geom_bin_2d(bins=100)+
        #scale_fill_gradientn(colours=rev(rainbow(3)))+
        scale_fill_gradientn(colours = c("black","white","red","orange","yellow","green","blue"),
                             values = c(1.0,0.5,0.2,0.15,0.1,0.05,0))+
        ylim(1,0)+xlim(0,1)+coord_fixed(ratio = 9/16)+
        #facet_wrap(~timepoint)+
        theme_bw()+
        theme(legend.position="none")+
        #animation specific
        #transition_time(time_s)+facet_wrap(~timepoint)
        transition_time(time_s)+
        ggtitle('Gaze behavior between groups','time: {frame_time}')

      gif_1<-animate(animated_gaze, duration = 10, fps = 10, width = 800, height = 800, renderer = gifski_renderer())
      anim_save(file=paste0(home_path,project_path,'/output/gaze_animation_las2s.gif'),gif_1)


### --> visualize RJA ####

##RJA
ggplot(df_jointatt,aes(x=rja_when/300))+
  geom_density(fill='orange')+
  #geom_freqpoly()+
  labs(x='video duration (s)',y='RJA occurances (density)')+
  xlim(c(0,12))+
  geom_vline(xintercept=1)+geom_vline(xintercept=3.5,linetype='dashed')+geom_vline(xintercept=9.5,linetype='dashed')+
  theme_bw()

#change facotr ordering
#df_jointatt$timepoint <- factor(df_jointatt$timepoint, levels = c("T2", "T4", "T6", "FU2", "FU3","K","K_FU"))
##RJA between timepoints
df_jointatt$group<-with(df_jointatt,ifelse(grepl('K',timepoint),'NTC','ASD'))
ggplot(df_jointatt[df_jointatt$timepoint!='FU3',],aes(x=rja_when/300,group=timepoint,fill=group))+
  geom_density()+
  #geom_freqpoly()+
  labs(x='video duration (s)',y='RJA occurances (density)')+
  xlim(c(0,11))+
  geom_vline(xintercept=1)+geom_vline(xintercept=3.5,linetype='dashed')+geom_vline(xintercept=9.5,linetype='dashed')+
  facet_grid(vars(timepoint))+theme_bw()

##head gaze between timepoints
ggplot(df_jointatt,aes(x=gazehead_when/300,group=timepoint,fill=timepoint))+
  geom_density()+
  #geom_freqpoly()+
  labs(x='video duration (s)',y='head gaze (density)')+
  xlim(c(0,11))+
  geom_vline(xintercept=1)+geom_vline(xintercept=4,linetype='dashed')+geom_vline(xintercept=10,linetype='dashed')+
  facet_grid(vars(timepoint))+theme_bw()

##head gaze between timepoints
ggplot(df_jointatt,aes(x=gazetarget_when/300,group=timepoint,fill=timepoint))+
  geom_density()+
  #geom_freqpoly()+
  labs(x='video duration (s)',y='target gaze (density)')+
  xlim(c(0,11))+
  geom_vline(xintercept=1)+geom_vline(xintercept=4,linetype='dashed')+geom_vline(xintercept=10,linetype='dashed')+
  facet_grid(vars(timepoint))+theme_bw()


##RJA between conditions
ggplot(df_jointatt[!is.na(df_jointatt$gaze_top),],aes(x=ts_event,fill=gaze_top))+geom_bar(position = "fill")+
  labs(x='sample (1/300s)',y='proportion of all samples')+
  geom_vline(xintercept=450,lty=1,color='black')+
  geom_vline(xintercept=900,lty=2,color='blue')+
  geom_vline(xintercept=2700,lty=2,color='blue')+
  facet_grid(vars(position),vars(condition))+theme_bw()

ggplot(df_jointatt[df_jointatt$fixation & !is.na(df_jointatt$gaze_target),],aes(x=ts_event,fill=gaze_target))+geom_bar(position = "fill")+
  labs(x='sample (1/300s)',y='proportion of all samples')+
  geom_vline(xintercept=900,lty=2,color='blue')+
  geom_vline(xintercept=2700,lty=2,color='blue')+
  geom_hline(yintercept=0.25,lty=1,color='black')+
  facet_grid(vars(condition),vars(position))+
  theme_bw()

### --> visualize pupillary response progression (for INSAR 2023 abstract) ####
#pupil size progression - between conditions

#--> overall pupillary response
ggplot(df_jointatt[sample(1:nrow(df_jointatt),size=nrow(df_jointatt)/100),],aes(x=ts_event/300,y=rpd))+
  geom_smooth()+theme_bw()+labs(title='pupil size within trial',x='time (samples)',y='pupil size (mm)')+
  geom_vline(xintercept=1.1)+geom_vline(xintercept=3,linetype=2)+geom_vline(xintercept=9,linetype=2)
  ###--> pupillary response is not a six degree polynomial

#--> between timepoints
df_jointatt<-df_jointatt[df_jointatt$timepoint!='FU3',] #remove FU3
ggplot(df_jointatt[sample(1:nrow(df_jointatt),size=nrow(df_jointatt)/5),],aes(x=ts_event/300,y=rpd,group=timepoint,color=timepoint))+
  geom_smooth(aes(fill=timepoint))+theme_bw()+labs(title='pupil size within trial',x='video duration (s)',y='pupil size change (mm)')+
  geom_vline(xintercept=1)+geom_vline(xintercept=3,linetype='dashed')+geom_vline(xintercept=9,linetype='dashed')

#--> between condition
ggplot(df_jointatt[sample(1:nrow(df_jointatt),size=nrow(df_jointatt)/10),],aes(x=ts_event/300,y=rpd,group=condition,color=condition))+
  geom_smooth()+theme_bw()+labs(title='pupil size within trial',x='time (samples)',y='pupil size (mm)')+
  geom_vline(xintercept=1)+geom_vline(xintercept=3,linetype='dashed')+scale_color_brewer(palette='Dark2')

#--> between condition x timepoint
ggplot(df_jointatt[sample(1:nrow(df_jointatt),size=nrow(df_jointatt)/10),],aes(x=ts_event,y=rpd,group=timepoint,color=timepoint))+
  geom_smooth()+theme_bw()+labs(title='pupil size within trial',x='time (samples)',y='pupil size (mm)')+
  geom_vline(xintercept=450)+geom_vline(xintercept=900)+facet_wrap(~condition)
  #--> intense: only t2 has attenuated response / mild: K response is stronger than t2/t4/t6

#-> condition x group
ggplot(df_jointatt,aes(x=ts_event,y=rpd,group=interaction(condition,group),color=condition,linetype=group))+
  geom_smooth(method='lm', formula = y ~ x + poly(x,6))+theme_bw()+labs(title='pupil size within trial',x='time (samples)',y='pupil size (mm)')+
  geom_vline(xintercept=450)+geom_vline(xintercept=900)

### --> create TRIAL AGGREGATED DATA set (df_trial) #####

#split by trial per participant
split_by_trial<-split(df_jointatt,interaction(df_jointatt$id,df_jointatt$index_trial))

cueing_onset<-1050
#time to target
time_to_target<-pbsapply(split_by_trial,function(x,stimulus_onset=cueing_onset,frequency=300){

  relevant_data<-x$ts_event[x$fixation & x$gaze_target & x$ts_event>stimulus_onset]

  ifelse(all(is.na(relevant_data)),
         NA,
         (min(relevant_data,na.rm=T)-stimulus_onset)/frequency)

})

#time to distractor
time_to_distractor<-pbsapply(split_by_trial,function(x,stimulus_onset=cueing_onset,frequency=300){

  relevant_data<-x$ts_event[x$fixation & x$gaze_distractor & x$ts_event>stimulus_onset]

  ifelse(all(is.na(relevant_data)),
         NA,
         (min(relevant_data,na.rm=T)-stimulus_onset)/frequency)

})


        ##DIAGNOSTICS
        # hist(time_to_target,50,main='gaze: time to target',xlab='duration (s)')
        # hist(time_to_distractor,50,main='gaze: time to distractor',xlab='duration (s)')
        # summary(time_to_target)
        # summary(time_to_distractor)
        # table(time_to_target)
        #
        # table(is.na(time_to_target))

#pupillary response - rpd is already corrected for baseline
## --> cutoff defined by pupillary respons progression
rpd_early<-pbsapply(split_by_trial,function(x){mean(x$rpd[x$ts_event>300 & x$ts_event<600],na.rm=T)})
rpd_middle<-pbsapply(split_by_trial,function(x){mean(x$rpd[x$ts_event>1350 & x$ts_event<1650],na.rm=T)})
rpd_late<-pbsapply(split_by_trial,function(x){mean(x$rpd[x$ts_event>3000 & x$ts_event<3300],na.rm=T)})

#gaze durations - based on fixations
frequency<-300
time_per_sample<-1/frequency
gazeduration_fixcross<-pbsapply(split_by_trial,function(x){table(x$gaze_fixcross[x$fixation])[2]*time_per_sample})
gazeduration_head<-pbsapply(split_by_trial,function(x){table(x$gaze_head[x$fixation])[2]*time_per_sample})
gazeduration_headbefore<-pbsapply(split_by_trial,function(x){table(x$gaze_head_before[x$fixation])[2]*time_per_sample})
gazeduration_headafter<-pbsapply(split_by_trial,function(x){table(x$gaze_head_after[x$fixation])[2]*time_per_sample})
gazeduration_stimulusleft<-sapply(split_by_trial,function(x){table(x$gaze_stimulus_left[x$fixation])[2]*time_per_sample})
gazeduration_stimulusright<-sapply(split_by_trial,function(x){table(x$gaze_stimulus_right[x$fixation])[2]*time_per_sample})
gazeduration_stimulus<-ifelse(is.na(gazeduration_stimulusleft),gazeduration_stimulusright,gazeduration_stimulusleft)
gazeduration_target<-sapply(split_by_trial,function(x){table(x$gaze_target[x$fixation])[2]*time_per_sample})
gazeduration_distractor<-sapply(split_by_trial,function(x){table(x$gaze_distractor[x$fixation])[2]*time_per_sample})
target_missed<-ifelse(gazeduration_target<0.1,T,F)

gazehead_t_mean<-sapply(split_by_trial,function(x){mean(x$gazehead_when,na.rm=T)})
gazetop_t_mean<-sapply(split_by_trial,function(x){mean(x$gazetop_when,na.rm=T)})
gazetarget_t_mean<-sapply(split_by_trial,function(x){mean(x$gazetarget_when,na.rm=T)})
gazedistractor_t_mean<-sapply(split_by_trial,function(x){mean(x$gazedistractor_when,na.rm=T)})


#rja parameters
rja_any<-sapply(split_by_trial,function(x){any(x$rja,na.rm=T)})
rja_any_early<-sapply(split_by_trial,function(x){any(x$rja_early,na.rm=T)})
rja_any_late<-sapply(split_by_trial,function(x){any(x$rja_late,na.rm=T)})

rja_duration<-sapply(split_by_trial,function(x){table(x$rja[x$fixation])[2]*time_per_sample})
rja_duration_early<-sapply(split_by_trial,function(x){table(x$rja_early[x$fixation])[2]*time_per_sample})
rja_duration_late<-sapply(split_by_trial,function(x){table(x$rja_late[x$fixation])[2]*time_per_sample})

rja_occurances<-sapply(split_by_trial,function(x){table(diff(x$rja)==1)[2]}) #diff of logical return 1 for TRUE after FALSE - how often equals occurances of rja
rja_t_mean<-sapply(split_by_trial,function(x){mean(x$rja_when,na.rm=T)})
rja_t_sd<-sapply(split_by_trial,function(x){sd(x$rja_when,na.rm=T)})

#saccade duration, fixation duration
saccade_duration<-sapply(split_by_trial,function(x){table(x$saccade)[2]*time_per_sample})
fixation_duration<-sapply(split_by_trial,function(x){table(x$fixation)[2]*time_per_sample})

#aggregate factor data (processing with R3.6 does not copy level labels)
id<-as.factor(with(df_jointatt,by(id,interaction(id,index_trial),unique)))
pic<-as.factor(with(df_jointatt,by(pic,interaction(id,index_trial),head,n=1)))
group<-as.factor(with(df_jointatt,by(group,interaction(id,index_trial),head,n=1)))
timepoint<-as.factor(with(df_jointatt,by(timepoint,interaction(id,index_trial),head,n=1)))
index_trial<-as.factor(with(df_jointatt,by(index_trial,interaction(id,index_trial),head,n=1)))
condition<-as.factor(with(df_jointatt,by(condition,interaction(id,index_trial),head,n=1)))
stimulus<-as.factor(with(df_jointatt,by(stimulus,interaction(id,index_trial),head,n=1)))
position<-as.factor(with(df_jointatt,by(position,interaction(id,index_trial),head,n=1)))

#aggregated mean numeric data
baseline_pd<-as.numeric(with(df_jointatt,by(baseline_pd,interaction(id,index_trial),mean,na.rm=T)))

#data quality
no_data<-sapply(split_by_trial,function(x){table(is.na(x$gazepos_x))[2]/length(x$gazepos_x)})
#no_data<-sapply(split_by_trial,function(x){sum(is.na(x$gazepos_x))/length(x$gazepos_x)})

df_trial<-data.frame(id,pic,group,timepoint,index_trial,condition,stimulus,position,
                     baseline_pd,rpd_early,rpd_middle,rpd_late,
                     time_to_target,saccade_duration,fixation_duration,
                     gazeduration_fixcross,gazeduration_head,gazeduration_headafter,gazeduration_headbefore,
                     gazeduration_stimulus,gazeduration_target,gazeduration_distractor,target_missed,
                     gazehead_t_mean,gazetop_t_mean,gazetarget_t_mean,gazedistractor_t_mean,
                     rja_any,rja_any_early,rja_any_late,rja_duration,rja_duration_early,rja_duration_late,
                     rja_occurances,rja_t_mean,rja_t_sd,no_data)


###---> merge demographics and df_trial - corrects for matched participants ####

warning('did you match data, then look here + latest dmeographic file required')
# #for matched
# df_trial<-merge(df_trial,all.match[,!(names(all.match)=='group')],by.x='pic',by.y='id',all.x=T)
# df_jointatt<-merge(df_jointatt,all.match[,!(names(all.match)=='group')],by.x='pic',by.y='id')

#for unmatched --> latest demographic file required (CURRENTLY NOT AVAILABLE)
df_trial<-merge(df_trial,demogr,by.x='pic',by.y='id',all.x=T)
df_jointatt<-merge(df_jointatt,demogr,by.x='pic',by.y='id',all.x = T)

###--> save final data frame ####

save(df_jointatt, df_trial, demogr, file=paste0(home_path,project_path,"/data/all_data_preprocessed_FINAL_jointatt_140324.Rdata"))

  ### --------------- ####

  ## CREATE DF_TIMEPOINT --> calcualte age at timepoints and date of measurement####

  #data frame aggegated by timepoint
  unique_id<-unique(df_trial$id)
  unique_id<-unique_id[complete.cases(unique_id)] #remove NA
  pic<-substr(unique_id,1,3)
  timepoint<-substr(unique_id,5,8)
  timepoint<-ifelse(grepl('t2',timepoint) | grepl('T2',timepoint),'T2',
                    ifelse(grepl('t4',timepoint) | grepl('T4',timepoint),'T4',
                           ifelse(grepl('t6',timepoint) | grepl('T6',timepoint),'T6',
                                  ifelse(grepl('fu2',timepoint) | grepl('FU2',timepoint),'FU2',
                                         ifelse(grepl('fu3',timepoint) | grepl('FU3',timepoint),'FU3',
                                                ifelse(grepl('K_fu',timepoint) | grepl('K_FU',timepoint) | grepl('k_FU',timepoint) | grepl('k_fu',timepoint),'K_FU','K'))))))


  df_timepoint<-data.frame(unique_id,pic,timepoint)
  #df_timepoint<-merge(df_timepoint,all.match,by.x='pic',by.y='id')
  df_timepoint<-merge(df_timepoint,demogr,by.x='pic',by.y='id')

  df_timepoint$timepoint <- factor(df_timepoint$timepoint, levels = c("FU2", "FU3","K","K_FU","T2", "T4", "T6"))


  table(df_timepoint$timepoint)

  ###--> per participant data-frame

  #### - date of measurement #####

  date_of_measurement<-as.numeric(with(df_jointatt,by(timestamp,id,median)))/1000000 #microsecond format to seconds
  date_of_measurement<-as.Date(as.POSIXct(date_of_measurement, origin="1970-01-01"))
  id_of_measurement<-as.character(with(df_jointatt,by(id,id,head,n=1)))
  df_date_measure<-data.frame(id=id_of_measurement,date=date_of_measurement)

  hist(date_of_measurement,20)
  df_timepoint<-merge(df_timepoint,df_date_measure,by.x='unique_id',by.y = 'id')

  df_timepoint$timepoint <- factor(df_timepoint$timepoint, levels = c("T2", "T4", "T6", "FU2", "FU3","K","K_FU"))

  ggplot(df_timepoint[df_timepoint$timepoint!='FU3',],aes(x=date,fill=timepoint))+geom_histogram(bins=50)+theme_bw()

  ggplot(df_timepoint[df_timepoint$timepoint!='FU3',],aes(x=pic,y=date,color=timepoint))+geom_point()+theme_bw()+labs(x='participant')



  #### - caluclate age at date ####

  list_initial_date<-with(df_timepoint,by(date,pic,min,simplify=FALSE))
  initial_date<-as.Date(sapply(list_initial_date,paste))
  df_initial_date<-data.frame(pic=names(list_initial_date),initial_date)
  df_timepoint<-merge(df_timepoint,df_initial_date,by='pic')
  difference_in_months<-(df_timepoint$date-df_timepoint$initial_date)/30

  df_timepoint$age_months_at_date<-as.numeric(df_timepoint$age_months+difference_in_months)

  #### --> merge df_timepoint to df_trial ####

  df_timepoint_merge<-df_timepoint[,!(names(df_timepoint) %in% names(df_trial))]
  df_trial<-merge(df_trial,df_timepoint_merge,by.x='id',by.y='unique_id',all.x=T)

  ## DESCRIPTIVE TABLE ####

  fun_return_descriptive<-function(variable,group){
    mean_values<-by(variable,group,function(x){round(mean(x,na.rm=T),2)})
    sd_values<-by(variable,group,function(x){round(sd(x,na.rm=T),2)})
    paste0(mean_values,'/',sd_values)
  }

  #order factors
  #df_timepoint$timepoint <- factor(df_timepoint$timepoint, levels = c("T2", "T4", "T6", "FU2", "FU3","K","K_FU"))

  #standard descriptives
  #n_timepoint<-with(df_trial,by(id,timepoint,function(x){length(unique(x))}))
  n_timepoint<-with(df_timepoint,by(pic,timepoint,function(x){length(unique(x))}))
  gender<-with(df_timepoint,by(Geschlecht_Index,timepoint,function(x){paste0(table(x)[2],'/',table(x)[1])}))
  test_age<-with(df_timepoint,fun_return_descriptive(variable=test_age,group=timepoint))
  age<-with(df_timepoint,fun_return_descriptive(variable=age_months,group=timepoint))
  age_at_date<-with(df_timepoint,fun_return_descriptive(variable=age_months_at_date,group=timepoint))
  iq<-with(df_timepoint,fun_return_descriptive(variable=IQ,group=timepoint))
  srs_16<-with(df_timepoint,fun_return_descriptive(variable=srs_sum,group=timepoint))
  rbsr<-with(df_timepoint,fun_return_descriptive(variable=RBSR_ges,group=timepoint))
  cbcl_ges<-with(df_timepoint,fun_return_descriptive(variable=CBCL_T_GES,group=timepoint))

  #main variable descriptives
  rja_duration<-with(df_trial,fun_return_descriptive(variable=rja_duration,group=timepoint))
  fixation_duration<-with(df_trial,fun_return_descriptive(variable=fixation_duration,group=timepoint))
  rpd_middle<-with(df_trial,fun_return_descriptive(variable=rpd_middle,group=timepoint))
  no_data<-with(df_trial,fun_return_descriptive(variable=no_data,group=timepoint))


  descriptives_table<-rbind(n_timepoint,gender,age_at_date,age,test_age,iq,srs_16,rbsr,cbcl_ges,
                            rja_duration,fixation_duration,rpd_middle,no_data)

  #remove FU3
  descriptives_table<-descriptives_table[,-5]


  #descriptives_table<-descriptives_table[,c(5,6,7,1,3,4)]

  row_names<-c('n','gender (F/M)','age at date (months)','age (at T2; months)','test age (at T2; months)','IQ (at T2)','SRS-16','RBSR','CBCL total (T)','RJA duration (s)','fixation duration (s)','pupillary response (z)','missing data (%)')


  descriptives_table<-cbind(row_names,descriptives_table)
  descriptives_table<-ifelse(descriptives_table=='0','<0.001',descriptives_table)


  require(kableExtra)
  table_sample<-descriptives_table %>%
    kbl(caption = "Sample description",
        col.names = c('','ASD (T2)','ASD (T4)','ASD (T6)','ASD (FU2)','NTC (T2)','NTC (FU2)'),
        row.names = F) %>%
    kable_classic(full_width = F, html_font = "Cambria")

  table_sample



  ####--> create df_ids (one row per individual)####

  table(df_timepoint$pic,df_timepoint$timepoint)

  #which measurement timepoin ts are available
  which_timepoints<-ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,4]==1 & table(df_timepoint$pic) %in% c(4,5),'T2+T4+T6+FU2',
                        ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,3]==1 & table(df_timepoint$pic)==3,'T2+T4+T6',
                           ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,2]==1 & table(df_timepoint$pic)==2,'T2+T4',
                                  ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,1]==1 & table(df_timepoint$pic)==1,'T2',
                                         ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,7]==1 & table(df_timepoint$pic)==2,'K+FU',
                                                ifelse(table(df_timepoint$pic,df_timepoint$timepoint)[,6]==1 & table(df_timepoint$pic)==1,'K','Other'))))))

  pic_participants<-names(table(df_timepoint$pic))
  df_ids<-data.frame(pic=pic_participants,which_timepoints)

  table(df_ids$which_timepoints)

  df_trial<-merge(df_trial,df_ids,by='pic',all.x=T) #add which_timepoints to df_trial
  df_ids<-merge(df_ids,df_timepoint[!duplicated(df_timepoint$pic),],by='pic')

  #add eye tracking metrics
  rja_any_p<-aggregate(df_trial$rja_any,by=list(df_trial$pic), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})
  rja_early_p<-aggregate(df_trial$rja_any_early,by=list(df_trial$pic), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})[,2]
  rja_late_p<-aggregate(df_trial$rja_any_late,by=list(df_trial$pic), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})[,2]
  rja_duration<-aggregate(df_trial$rja_duration,by=list(df_trial$pic), FUN=mean, na.rm=T)[,2]
  rpd_middle<-aggregate(df_trial$rpd_middle,by=list(df_trial$pic), FUN=mean, na.rm=T)[,2]
  pd_500<-aggregate(df_trial$baseline_pd,by=list(df_trial$pic), FUN=mean, na.rm=T)[,2]
  gazeduration_head<-aggregate(df_trial$gazeduration_head,by=list(df_trial$pic), FUN=mean, na.rm=T)[,2]

  df_ids_et<-data.frame(rja_any_p,rja_early_p,rja_late_p,rja_duration,rpd_middle,gazeduration_head,pd_500)
  names(df_ids_et)[1:2]<-c('pic','rja_any_p')

  df_ids<-merge(df_ids,df_ids_et,by='pic')

  ###calculate ET change variables (T6-T2)

  df_trial$timepoint<-with(df_trial,ifelse(timepoint==1,'T2',
                                           ifelse(timepoint==2,'T4',
                                                  ifelse(timepoint==3,'T6',
                                                         ifelse(timepoint==4,'FU2',
                                                                ifelse(timepoint==5,'FU3',
                                                                       ifelse(timepoint==6,'K',
                                                                              ifelse(timepoint==7,'K_FU',NA))))))))


  rja_any_p_T6<-with(df_trial[df_trial$timepoint=='T6',],aggregate(rja_any,by=list(pic), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])}))
  rja_any_p_T2<-with(df_trial[df_trial$timepoint=='T2',],aggregate(rja_any,by=list(pic), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])}))
  rja_any_p_change<-merge(rja_any_p_T6,rja_any_p_T2,by='Group.1')
  names(rja_any_p_change)<-c('pic','rja_any_p_T6','rja_any_p_T2')
  rja_any_p_change$rja_change_T6T2<-rja_any_p_change$rja_any_p_T6-rja_any_p_change$rja_any_p_T2

  rpd_middle_T6<-with(df_trial[df_trial$timepoint=='T6',],aggregate(rpd_middle,by=list(pic), FUN=mean, na.rm=T))
  rpd_middle_T2<-with(df_trial[df_trial$timepoint=='T2',],aggregate(rpd_middle,by=list(pic), FUN=mean, na.rm=T))
  rpd_middle_change<-merge(rpd_middle_T6,rpd_middle_T2,by='Group.1')
  names(rpd_middle_change)<-c('pic','rpd_middle_T6','rpd_middle_T2')
  rpd_middle_change$rpd_change_T6T2<-rpd_middle_change$rpd_middle_T6-rpd_middle_change$rpd_middle_T2

  gazeduration_head_T6<-with(df_trial[df_trial$timepoint=='T6',],aggregate(gazeduration_head,by=list(pic), FUN=mean, na.rm=T))
  gazeduration_head_T2<-with(df_trial[df_trial$timepoint=='T2',],aggregate(gazeduration_head,by=list(pic), FUN=mean, na.rm=T))
  gazeduration_head_change<-merge(gazeduration_head_T6,gazeduration_head_T2,by='Group.1')
  names(gazeduration_head_change)<-c('pic','gazeduration_head_T6','gazeduration_head_T2')
  gazeduration_head_change$gazeduration_change_T6T2<-gazeduration_head_change$gazeduration_head_T6-gazeduration_head_change$gazeduration_head_T2

  df_ids<-df_ids[,-(39:47)]

  df_ids<-merge(df_ids,rja_any_p_change,by='pic',all.x=T)
  df_ids<-merge(df_ids,rpd_middle_change,by='pic',all.x=T)
  df_ids<-merge(df_ids,gazeduration_head_change,by='pic',all.x=T)


  hist(df_ids$rja_change_T6T2)
  hist(df_ids$rpd_change_T6T2)
  hist(df_ids$gazeduration_change_T6T2)

  #-->predict intervention group
  df_ids$pic[df_ids$rja_change_T6T2>0 & df_ids$rpd_change_T6T2>0 & df_ids$gazeduration_change_T6T2>0]

  ###--> add ET data to df_timepoint ####
  rja_any_p<-aggregate(df_trial$rja_any,by=list(df_trial$id), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})
  rja_early_p<-aggregate(df_trial$rja_any_early,by=list(df_trial$id), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})[,2]
  rja_late_p<-aggregate(df_trial$rja_any_late,by=list(df_trial$id), FUN=function(x){table(x)[2]/(table(x)[2]+table(x)[1])})[,2]
  rja_duration<-aggregate(df_trial$rja_duration,by=list(df_trial$id), FUN=mean, na.rm=T)[,2]
  rpd_middle<-aggregate(df_trial$rpd_middle,by=list(df_trial$id), FUN=mean, na.rm=T)[,2]
  pd_500<-aggregate(df_trial$baseline_pd,by=list(df_trial$id), FUN=mean, na.rm=T)[,2]
  gazeduration_head<-aggregate(df_trial$gazeduration_head,by=list(df_trial$id), FUN=mean, na.rm=T)[,2]

  df_timepoint_et<-data.frame(rja_any_p,rja_early_p,rja_late_p,rja_duration,rpd_middle,gazeduration_head,pd_500)
  names(df_timepoint_et)[1:2]<-c('unique_id','rja_any_p')

  df_timepoint<-merge(df_timepoint,df_timepoint_et,by='unique_id')

  ###--> analysis of dropouts (ET data) ####

  #standard descriptives
  n_timepoint<-as.numeric(table(df_ids$which_timepoints))
  gender<-with(df_ids,by(Geschlecht_Index,which_timepoints,function(x){paste0(table(x)[2],'/',table(x)[1])}))
  test_age<-with(df_ids,fun_return_descriptive(variable=test_age,group=which_timepoints))
  age<-with(df_ids,fun_return_descriptive(variable=age_months,group=which_timepoints))
  iq<-with(df_ids,fun_return_descriptive(variable=IQ,group=which_timepoints))

  #clinical variables
  ADOS_CSS_SA<-with(df_ids,fun_return_descriptive(variable=Severity_SA,group=which_timepoints))
  ADOS_CSS_RRB<-with(df_ids,fun_return_descriptive(variable=Severity_RRB,group=which_timepoints))
  ADOS_CSS_TOTAL<-with(df_ids,fun_return_descriptive(variable=Severity_Gesamt,group=which_timepoints))
  cbcl_int<-with(df_ids,fun_return_descriptive(variable=CBCL_T_INT,group=which_timepoints))
  cbcl_ext<-with(df_ids,fun_return_descriptive(variable=CBCL_T_EXT,group=which_timepoints))
  cbcl_ges<-with(df_ids,fun_return_descriptive(variable=CBCL_T_GES,group=which_timepoints))

  #main variable descriptives
  rja_duration<-with(df_trial,fun_return_descriptive(variable=rja_duration,group=which_timepoints))
  fixation_duration<-with(df_trial,fun_return_descriptive(variable=fixation_duration,group=which_timepoints))
  rpd_middle<-with(df_trial,fun_return_descriptive(variable=rpd_middle,group=which_timepoints))
  no_data<-with(df_trial,fun_return_descriptive(variable=no_data,group=which_timepoints))

  descriptives_table_ids<-rbind(n_timepoint,gender,test_age,age,iq,
                            ADOS_CSS_SA,ADOS_CSS_RRB,ADOS_CSS_TOTAL,cbcl_int,cbcl_ext,cbcl_ges,
                            rja_duration,fixation_duration,rpd_middle,no_data)

  row_names<-c('n','gender (F/M)','test age (months)','age (months)','IQ',
               'ADOS CSS SA','ADOS CSS RRB','ADOS CSS TOTAL','CBCL int (T)','CBCL ext (T)','CBCL total (T)',
               'RJA duration (s)','fixation duration (s)','pupillary response (z)','missing data (%)')


  descriptives_table_ids<-cbind(row_names,descriptives_table_ids)
  descriptives_table_ids<-ifelse(descriptives_table_ids=='0','<0.001',descriptives_table_ids)

  descriptives_table_ids

  require(kableExtra)
  table_sample<-descriptives_table_ids %>%
    kbl(caption = "Sample description",
        col.names = c('','NTC (T2)','NTC (FU)','ASD (Other)','ASD (T2)','ASD (T2+T4)','ASD (T2+T4+T6)','ASD (T2+T4+T6+FU)'),
        row.names = F) %>%
    kable_classic(full_width = F, html_font = "Cambria")

  table_sample

  require(ggplot2)
  ggplot(df_ids[df_ids$which_timepoints %in% c('T2','T2+T4','T2+T4+T6'),],aes(y=test_age,x=which_timepoints,fill=which_timepoints))+
    geom_boxplot()+scale_fill_brewer(palette="Dark2")+theme_bw()

  ggplot(df_ids[df_ids$which_timepoints %in% c('T2','T2+T4','T2+T4+T6'),],aes(y=Severity_SA,x=which_timepoints,fill=which_timepoints))+
    geom_boxplot()+scale_fill_brewer(palette="Dark2")+theme_bw()

  ggplot(df_ids[df_ids$which_timepoints %in% c('T2','T2+T4','T2+T4+T6'),],aes(y=Severity_RRB,x=which_timepoints,fill=which_timepoints))+
    geom_boxplot()+scale_fill_brewer(palette="Dark2")+theme_bw()

  ggplot(df_trial[df_trial$which_timepoints %in% c('T2','T2+T4','T2+T4+T6'),],aes(y=rja_duration,x=which_timepoints,fill=which_timepoints))+
    geom_boxplot()+scale_fill_brewer(palette="Dark2")+theme_bw()


  #### --> DATA ANALYSIS  ####

 hist(df_jointatt$ts_event)
 hist(df_jointatt$pd)
 hist(df_jointatt$rpd)
 hist(df_trial$no_data) #excluded trials with less than 50% of data
 hist(df_trial$fixation_duration)

 table(df_jointatt$condition,df_jointatt$stimulus)
 table(df_jointatt$condition,df_jointatt$position)
 table(df_jointatt$stimulus,df_jointatt$position)

  ##--> correlation tables ####
 cor_matrix<-cor(df_trial[,c('rja_duration','rja_occurances','time_to_target','baseline_pd',
                           'rpd_early','rpd_middle','rpd_late','no_data',
                           'BOSCC_Average_Total','test_age','age_months','CBCL_T_GES',
                           'CBCL_T_INT','CBCL_T_EXT','srs_sum','RBSR_ges')],use='complete.obs')[]

 cor_matrix_format<-cor_matrix
 cor_matrix_format<-round(cor_matrix_format,2)
 cor_matrix_format[upper.tri(cor_matrix)]<-''
 cor_matrix_format<-as.data.frame(cor_matrix_format)

 cor_matrix_format

  # --> visualize pupillary size + response stability (Emmy Noether proposal) ####

 stat_sum_df <- function(fun, geom="crossbar", ...) {
   stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
 }

 df_trial$timepoint <- factor(df_trial$timepoint, levels = c("T2", "T4", "T6", "FU2", "FU3","K","K_FU"))
 #ggplot(df_trial[df_trial$pic %in% unique(df_trial$pic[df_trial$timepoint=='T6']),],
 ggplot(df_trial[df_trial$which_timepoints=='T2+T4+T6',])+
  aes(x=timepoint,y=baseline_pd,group=pic,color=pic)+
   stat_sum_df("mean_cl_normal",geom='point')+
   stat_summary(fun="mean",geom="line")+theme_bw()

 ggplot(df_trial[df_trial$pic %in% unique(df_trial$pic[df_trial$timepoint=='T6']),],
        aes(x=baseline_pd))+geom_density(fill='grey')+theme_bw()


 #select only individuals that have T6 data
 df_trial_repeated<-df_trial[df_trial$pic %in% unique(df_trial$pic[df_trial$timepoint=='T6']),]

 #get baseline means for T2,T4,T6 --> assess correlations and ICC
 baseline_pd_t2<-with(df_trial_repeated[df_trial_repeated$timepoint=='T2',],aggregate(baseline_pd,by=list(pic),FUN=mean,na.rm=T))
 baseline_pd_t4<-with(df_trial_repeated[df_trial_repeated$timepoint=='T4',],aggregate(baseline_pd,by=list(pic),FUN=mean,na.rm=T))
 baseline_pd_t6<-with(df_trial_repeated[df_trial_repeated$timepoint=='T6',],aggregate(baseline_pd,by=list(pic),FUN=mean,na.rm=T))
 names(baseline_pd_t2)<-c('pic','t2')
 names(baseline_pd_t4)<-c('pic','t4')
 names(baseline_pd_t6)<-c('pic','t6')

 baseline_pd<-merge(baseline_pd_t2,baseline_pd_t4,by='pic')
 baseline_pd<-merge(baseline_pd,baseline_pd_t6,by='pic')

 cor.test(baseline_pd$t2,baseline_pd$t4)
 cor.test(baseline_pd$t4,baseline_pd$t6)
 cor.test(baseline_pd$t2,baseline_pd$t6)

 psych::ICC(baseline_pd[,-1])

 #get pupillary response means for T2,T4,T6 --> assess correlations and ICC
 rpd_t2<-with(df_trial_repeated[df_trial_repeated$timepoint=='T2',],aggregate(rpd_middle,by=list(pic),FUN=mean,na.rm=T))
 rpd_t4<-with(df_trial_repeated[df_trial_repeated$timepoint=='T4',],aggregate(rpd_middle,by=list(pic),FUN=mean,na.rm=T))
 rpd_t6<-with(df_trial_repeated[df_trial_repeated$timepoint=='T6',],aggregate(rpd_middle,by=list(pic),FUN=mean,na.rm=T))
 names(rpd_t2)<-c('pic','t2')
 names(rpd_t4)<-c('pic','t4')
 names(rpd_t6)<-c('pic','t6')

 rpd<-merge(rpd_t2,rpd_t4,by='pic')
 rpd<-merge(rpd,rpd_t6,by='pic')

 cor.test(rpd$t2,rpd$t4)
 cor.test(rpd$t4,rpd$t6)
 cor.test(rpd$t2,rpd$t6)

 psych::ICC(rpd[,-1])


  ### RJA by age ####

 ggplot(df_timepoint,aes(x=age_months_at_date,y=rja_any_p,color=group,fill=group))+geom_smooth(method='lm', formula = y ~ x + poly(x,4))+theme_bw()+
   labs(x='age (months)',y='RJA likelihood (%)')
 ggplot(df_timepoint,aes(x=age_months_at_date,y=rja_any_p,color=timepoint))+geom_point()


 table(df_trial$timepoint)
 table(df_timepoint$age_months,df_timepoint$timepoint)

 with(df_timepoint,by(age_months,timepoint,mean))

 #define age at measurement timepoint
 df_trial$age_at_timepoint<-with(df_trial,ifelse(timepoint=='T4',age_months+6,
                                        ifelse(timepoint=='T6',age_months+12,age_months)))

 #social attention - head gaze duration
 ggplot(df_trial,aes(x=age_at_timepoint,y=gazeduration_head,group=timepoint,color=timepoint))+
   geom_point()+
   geom_smooth(method='lm',formula = y ~ x + I(x^2),aes(fill=timepoint))+
   theme_bw()

 #rja duration
 ggplot(df_trial,aes(x=age_at_timepoint,y=rja_duration,group=timepoint,color=timepoint))+
   geom_point()+
   geom_smooth(method='lm',formula = y ~ x + I(x^2),aes(fill=timepoint))+
   theme_bw()

 ggplot(df_trial,aes(x=srs_sum,y=rja_duration,group=timepoint,color=timepoint,fill=timepoint))+
   geom_point()+
   geom_smooth(method='lm',formula = y ~ x + I(x^2))+
   theme_bw()

 ggplot(df_trial,aes(x=age_at_timepoint,y=rja_any))+
   geom_bar(stat='identity')+
   theme_bw()


 #LINEAR MIXED MODELS ####

 with(df_trial,ftable(rja_any,timepoint))
 hist(df_trial$rja_duration,30)
 hist(df_trial$gazeduration_headbefore)

 with(df_trial,round(table(timepoint,rja_any)/colSums(table(rja_any,timepoint)),2)) %>%
   kbl(caption = "RJA likelihood per trial") %>%
   kable_classic(full_width = F, html_font = "Cambria")

  ##--> Polynomial time models - unaggregated data ####

 #bin 100ms timespans in new variable
 frame_rate<-300
 df_jointatt$time_s<-with(df_jointatt,round(ts_event/frame_rate,1))

 #-->pupillary response plot looks like 6th degree polynomial
 lmm<-lmer(scale(rpd)~timepoint*
             (scale(time_s)+scale(I(time_s^2))+scale(I(time_s^3))+scale(I(time_s^4))+scale(I(time_s^5))+scale(I(time_s^6)))+ #temporal progression
             condition+stimulus+
              (1|index_trial)+(1|pic),data=df_jointatt[df_jointatt$timepoint!='FU3',])
 ###--> takes a few minutes
 anova(lmm)

 # ##task effects
 contrast(emmeans(lmm,~condition),'pairwise') #--> intense and neutral associated with higher rpd than mild and point
 contrast(emmeans(lmm,~stimulus),'pairwise') #--> intense and neutral associated with higher rpd than mild and point

 contrast(emmeans(lmm,~timepoint|time_s,
                  at=list(time_s =c(0,0.5,1,1.5,2,2.5,3,3.5,4.0,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10)),
                  params='deg'),'pairwise')


 #plot interaction
 model_plot_data<-as.data.frame(emmeans(lmm,~timepoint+time_s,
                                            at=list(time_s =c(0,0.5,1,1.5,2,2.5,3,3.5,4.0,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10)),
                                            params='deg'))
  #--> in the 4-6 second range T6 and K have higher pupillary response than T2


 ggplot(model_plot_data,aes(x=as.factor(time_s),y=emmean,group=timepoint,color=timepoint))+
   # geom_boxplot(aes(fill=timepoint,
   #                  middle=emmean,
   #                  lower=emmean-SE,
   #                  upper=emmean+SE,
   #                  ymin=emmean-2*SE,
   #                  ymax=emmean+2*SE),stat = "identity")+
   geom_errorbar(aes(ymin=emmean-SE,ymax=emmean+SE))+
   geom_line()+
   geom_point()+
   #geom_ribbon(aes(ymin = emmean-SE, ymax = emmean+SE, fill=timepoint), alpha = 0.1, linetype=2, outline.type = 'full')+
   theme_bw()+
   xlab('time (s)')+ylab('effect size (z)')
   #scale_fill_discrete(name = "timepoint", labels = c("K", "T2","T4","T6"))+ #change legend labels
   #theme(legend.position = c(0.2, 0.2),legend.box.background = element_rect(color="black",size=1.1)) #legend position and box


 ##---> data quality between group --> not used as covariate ####

 lmm<-lmer(no_data~timepoint+
             condition+
             (1|pic),data=df_trial)
 anova(lmm)
 emmeans(lmm,~timepoint)
 contrast(emmeans(lmm,~timepoint),'pairwise')
 ###--> no differences in data quality between groups, although significant effect

 lmm<-lmer(scale(rpd_middle)~scale(no_data)+
             condition+
             (1|pic),data=df_trial)
 anova(lmm)
 ###---> pupillary response is unrelated to data quality --> no used as covariate

 lmm<-glmer(rja_any~scale(no_data)+
              (1|pic),data=df_trial,family='binomial')
 summary(lmm)
 exp(fixef(lmm)['scale(no_data)'])
 ###--> more missing data is assoicated with lower likelihood of rja --> thus used in rja models, but not in rpd models


 ##--> RJA - MAIN FINDING - on per-trial level ####

  df_lmm<-df_trial[df_trial$timepoint!='FU3',]
  with(df_lmm,table(rja_any,timepoint))

  lmm<-glmer(rja_any~timepoint+
                #condition+stimulus+
                scale(no_data)+
                (1|index_trial)+(1|pic),data=df_lmm,family=binomial(link = "logit"))

   summary(lmm)
   anova(lmm)

   require(lattice)
   qqmath(ranef(lmm)) ##qqplot of trials and PIC


   round(anova(lmm),2) %>%
     kbl(caption = "RJA model") %>%
     kable_classic(full_width = F, html_font = "Cambria")


   emmeans(lmm,~timepoint)

   confint(emmeans(lmm,~timepoint))
   contrast(emmeans(lmm,~timepoint),'revpairwise')

   #translated to Odd ratio (and ) - T2 vs T6
   exp(-(confint(contrast(emmeans(lmm,~timepoint),'pairwise'))[14,2])) #T2 - T6
   exp(-(confint(contrast(emmeans(lmm,~timepoint),'revpairwise'))[4,2])) #FU2 - T2


   plot(emmeans(lmm,~timepoint))+theme_bw()
   ###--> control children more joint attention than all ASD groups
   ##--> ASD group improves over time (better joint attention after one year)

   #create data frame with predicted RJA - log odds
   predicted_loglik_RJA<-predict(lmm)
   predicted_RJA_mean<-as.numeric(with(df_lmm[!is.na(df_lmm$no_data),],by(predicted_loglik_RJA,id,mean)))
   timepoint<-as.factor(with(df_lmm[!is.na(df_lmm$no_data),],by(timepoint,id,head,n=1)))
   group<-as.factor(with(df_lmm[!is.na(df_lmm$no_data),],by(group,id,head,n=1)))

   df_plot_predicted_RJA<-data.frame(predicted_RJA_mean,timepoint,group)

            #data frame of estimated marginalized means
            emm_data<-data.frame(plot(emmeans(lmm,~timepoint))[[1]])
            emm_data$timepoint<-factor(emm_data$timepoint, levels = c("T2", "T4", "T6", "FU2", "FU3","K","K_FU"))
            emm_data$group<-with(emm_data,ifelse(grepl('K',timepoint),'NTC','ASD'))


            ggplot(data=emm_data,aes(x=timepoint))+

              #add conventional error bars
              geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),width=0.2)+


              #overplot predicted values
              geom_jitter(data=df_plot_predicted_RJA,
                          aes(x=timepoint,
                              y=predicted_RJA_mean,color=group),alpha=0.8,width =0.4)+


              geom_boxplot(aes(fill=group,
                               middle=the.emmean,
                               lower=the.emmean-1.5*SE,
                               upper=the.emmean+1.5*SE,
                               ymin=asymp.LCL,
                               ymax=asymp.UCL),stat = "identity",alpha=0.7)+

              #add significance comaprison in figure
              geom_signif(stat="identity",
                          data=data.frame(x=c(1,1,3.5,3,3,4), xend=c(3.5,1,3.5,4,3,4), y=c(2,2,2,1.75,1.75,1.75), yend=c(2,1.5,1.75,1.75,1.5,1.5), annotation="*"),
                          aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation))+


              labs(x='groups',y='RJA likelihood (z)')+
              #scale_fill_discrete(name = "timepoint", labels = c("T2","T4","T6","FU2","K","K_FU"))+
              scale_x_discrete(limits = c('T2','T4','T6','FU2','K','K_FU'))+ #order x axis
              theme_bw()
              #theme(legend.position="none")


            ggplot(df_timepoint,aes(x=timepoint,y=rja_any_p,color=pic))+geom_point()+geom_line()

            stat_sum_df <- function(fun, geom="crossbar", ...) {
              stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
            }


            #group with intervention?
            ggplot(df_timepoint[df_timepoint$pic %in% unique(df_timepoint$pic[df_timepoint$timepoint=='T6']) & df_timepoint$timepoint %in% c('T2','T6'),],
                   aes(x=timepoint,y=rja_any_p,group=pic,color=pic))+
              stat_sum_df("mean_cl_normal",geom='point')+
              stat_summary(fun="mean",geom="line")+theme_bw()


            ggplot(df_timepoint[df_timepoint$pic %in% unique(df_timepoint$pic[df_timepoint$timepoint=='T6']) & df_timepoint$timepoint %in% c('T2','T6'),],
                   aes(x=timepoint,y=rpd_middle,group=pic,color=pic))+
              stat_sum_df("mean_cl_normal",geom='point')+
              stat_summary(fun="mean",geom="line")+theme_bw()

            ggplot(df_timepoint[df_timepoint$pic %in% unique(df_timepoint$pic[df_timepoint$timepoint=='T6']) & df_timepoint$timepoint %in% c('T2','T6'),],
                   aes(x=timepoint,y=gazeduration_head,group=pic,color=pic))+
              stat_sum_df("mean_cl_normal",geom='point')+
              stat_summary(fun="mean",geom="line")+theme_bw()




    with(df_trial,by(rja_duration,timepoint,psych::describe))
    with(df_trial,by(time_to_target,timepoint,psych::describe))
    ###--> normalize in time-to-target
        with(df_trial,by(fixation_duration,timepoint,psych::describe))
        with(df_trial,by(gazeduration_head,timepoint,psych::describe))
        with(df_trial,by(gazeduration_headbefore,timepoint,psych::describe))
        with(df_trial,by(gazeduration_target,timepoint,psych::describe))
        ###--> seems like shorter gazes on target and longer gazes on head after therapy

        ## RJA early --> in first 4 seconds of video
        df_lmm<-df_trial[df_trial$timepoint!='FU3',]
        lmm<-glmer(rja_any_early~timepoint+
                     scale(no_data)+
                     (1|index_trial)+(1|pic),data=df_lmm,family=binomial(link = "logit"))
        anova(lmm)
        confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
        exp(-(confint(contrast(emmeans(lmm,~timepoint),'pairwise'))[14,2])) #T2 - T6


        #RJA late --> after 4 seconds until 7 seconds
        lmm<-glmer(rja_any_late~timepoint+
                     scale(no_data)+condition+
                     (1|index_trial)+(1|pic),data=df_lmm,family=binomial(link = "logit"))
        anova(lmm)
        confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
        ###--> T6 improves compared to T2 in LATE RJA


        lmm<-lmer(rja_duration~timepoint*condition+
                    scale(no_data)+(1|pic),data=df_trial)

        round(anova(lmm),3)%>%
          kbl(caption = "RJA duration") %>%
          kable_classic(full_width = F, html_font = "Cambria")


        plot(contrast(emmeans(lmm,~timepoint),'eff'))+theme_bw()
        ###--> longer rja duration in K versus T2,T4,T6

        lmm<-lmer(rja_occurances~timepoint*condition+
                    scale(no_data)+(1|pic),data=df_trial)

        round(anova(lmm),3)%>%
          kbl(caption = "RJA occurances") %>%
          kable_classic(full_width = F, html_font = "Cambria")

        plot(contrast(emmeans(lmm,~timepoint),'eff'))+theme_bw()
        ###--> less occurances of rja in T6 compared to T2/K


    ### RJA per trial ####

        ##--> does not converge for rja_any

    lmm<-lmer(scale(gazeduration_head)~as.numeric(index_trial)*timepoint+
                 scale(no_data)+
                 (1|pic),data=df_trial)
        anova(lmm)


        lmm<-lmer(scale(rja_duration)~as.numeric(index_trial)*timepoint+
                    scale(no_data)+
                    (1|pic),data=df_trial)
        anova(lmm)




        ### --> other measures (e.g.: social attention) ####
      df_lmm<-df_trial[df_trial$timepoint!='FU3',]
      lmm<-lmer(gazeduration_head~timepoint*condition+
                scale(no_data)+(1|pic),data=df_lmm)

        round(anova(lmm),3)%>%
          kbl(caption = "gaze duration head") %>%
          kable_classic(full_width = F, html_font = "Cambria")

    plot(contrast(emmeans(lmm,~timepoint),'eff'))+theme_bw()
    confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
    ###--> control children more joint attention than all ASD groups
    ##--> ASD group improves over time (more attention to head)
    plot(contrast(emmeans(lmm,~condition),'eff'))
    ###--> mild condition associated with most gazes on head
    with(df_trial,by(gazeduration_head,timepoint,psych::describe))

    #when does gaze on head occur --> later after therapy
    lmm<-lmer(scale(gazetop_t_mean)~timepoint*condition+
                scale(no_data)+(1|pic),data=df_trial)
    anova(lmm)
    confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
    plot(contrast(emmeans(lmm,~timepoint),'eff'))


    lmm<-lmer(gazeduration_headafter~timepoint*condition+
                scale(no_data)+(1|pic),data=df_lmm)
    anova(lmm)
    plot(contrast(emmeans(lmm,~timepoint),'eff'))
    plot(contrast(emmeans(lmm,~condition),'eff'))
    confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
    #--> only T2 has shortr head after gazeduration than K
    #--> thus normalizes in T4 and T6

    lmm<-lmer(gazeduration_headbefore~timepoint*condition+
                scale(no_data)+(1|pic),data=df_trial)
    anova(lmm)
    plot(contrast(emmeans(lmm,~timepoint),'eff'))
    confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
    #-->  K has longer gazeduration on head before cueing onset

          lmm<-glmer(target_missed~timepoint+
                       scale(no_data)+(1|pic),data=df_trial,family=binomial(link = "logit"))
          anova(lmm)
          table(df_trial$target_missed)
          ###--> not a good metric as ceiling effect


          lmm<-lmer(gazeduration_target~timepoint*condition+
                      scale(no_data)+(1|pic),data=df_trial)
          anova(lmm)
          #-->  no differences in gazeduration on target


          lmm<-lmer(scale(rja_t_sd)~timepoint*condition+
                       scale(no_data)+(1|pic),data=df_trial)
          anova(lmm)
          plot(contrast(emmeans(lmm,~timepoint),'pairwise'))


          lmm<-lmer(no_data~timepoint*condition+
                      (1|pic),data=df_trial)
          anova(lmm)
          confint(contrast(emmeans(lmm,~timepoint),'pairwise'))
          ###--> T6 has less data than T2, but does not differ to K

  ##--> pupillary response by timepoint and condition ####


    lmm<-lmer(scale(baseline_pd)~timepoint+
                condition+stimulus+
                (1|pic),data=df_trial)

    round(anova(lmm),3) %>%
      kbl(caption = "baseline pupil size") %>%
      kable_classic(full_width = F, html_font = "Cambria")


    lmm<-lmer(scale(rpd_early)~timepoint+condition+
                scale(no_data)+(1|pic),data=df_trial)
    anova(lmm)
    ggplot(df_trial,aes(x=condition,y=rpd_early,group=interaction(timepoint,condition),fill=timepoint))+geom_boxplot()

    df_lmm<-df_trial[df_trial$timepoint!='FU3',]
    lmm<-lmer(scale(rpd_middle)~timepoint+
                (1|pic),data=df_lmm)

    round(anova(lmm),3) %>%
      kbl(caption = "pupillary response after cueing") %>%
      kable_classic(full_width = F, html_font = "Cambria")

    plot(contrast(emmeans(lmm,~timepoint),'eff'))
    ###--> marginally significant: higher middle rpd in K compared to T2
    ggplot(df_trial,aes(x=condition,y=rpd_middle,group=timepoint,fill=timepoint))+geom_violin()


    lmm<-lmer(scale(rpd_late)~timepoint+condition+
                scale(no_data)+(1|pic),data=df_trial)
    anova(lmm)


    # require(effectsize)
    # d_to_r(0.5)

  ##--> association with clinical phenotypes ####

  #Across groups
  #RJA
  lmm<-glmer(rja_any~scale(srs_sum)+(1|pic),data=df_trial[df_trial$timepoint %in% c('K','T2'),],family=binomial(link = "logit"))
  summary(lmm)
  exp(fixef(lmm)[2])

  lmm<-glmer(rja_any~scale(RBSR_ges)+(1|pic),data=df_trial[df_trial$timepoint %in% c('K','T2'),],family=binomial(link = "logit"))
  summary(lmm)
  exp(fixef(lmm)[2])

  lmm<-glmer(rja_any~scale(CBCL_T_GES)+(1|pic),data=df_trial[df_trial$timepoint %in% c('K','T2'),],family=binomial(link = "logit"))
  summary(lmm)
  exp(fixef(lmm)[2])
  ##--> higher quantitiative symptoms associated with lower likelihodd of rja

  #visualization

  require(gridExtra)

  ## extract legend
  g_legend<-function(x){
    tmp <- ggplot_gtable(ggplot_build(x))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

  my_legend<-g_legend(ggplot(df_ids,aes(srs_sum,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw())


  grid.arrange(
  ggplot(df_ids,aes(srs_sum,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()+theme(legend.position = 'none'),
  ggplot(df_ids,aes(CBCL_T_GES,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()+theme(legend.position = 'none'),
  ggplot(df_ids,aes(srs_sum,gazeduration_head))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()+theme(legend.position = 'none'),
  ggplot(df_ids,aes(CBCL_T_GES,gazeduration_head))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()+theme(legend.position = 'none'),
  my_legend,
  bottom='clinical phenotype',left= 'social attention')


  #Within ASD
  lmm<-glmer(rja_any~scale(BOSCC_Average_Total)+(1|pic),data=df_trial,family=binomial(link = "logit"))
  summary(lmm)
  exp(fixef(lmm)[2])

  ggplot(df_ids[df_ids$group=='ASD',],
         aes(BOSCC_Average_Total,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()

  ggplot(df_ids[df_ids$group=='ASD',],
         aes(Severity_Gesamt,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()

  ggplot(df_ids[df_ids$group=='ASD',],
         aes(IQ,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()

  ggplot(df_ids[df_ids$group=='ASD',],
         aes(age_months,rja_any_p))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()



  #GAZE HEAD
  lmm<-lmer(scale(gazeduration_head)~scale(srs_sum)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(srs_sum)')
  ggplot(df_trial,aes(x=srs_sum,y=gazeduration_head))+geom_smooth(method='lm')+geom_jitter(aes(color=timepoint))+theme_bw()
  #--> higher srs sum associated with less gaze duration on head

  lmm<-lmer(scale(gazeduration_head)~scale(RBSR_ges)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(srs_sum)')
  #--> bit not with RBSR - makes sense

  lmm<-lmer(scale(gazeduration_head)~scale(CBCL_T_GES)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(CBCL_T_GES)')
  ggplot(df_trial,aes(x=CBCL_T_GES,y=gazeduration_head))+geom_smooth(method='lm')+geom_jitter(aes(color=timepoint))+theme_bw()
  #--> also higher psychopathology is associated with shorter gaze duration on head

  #RPD
  lmm<-lmer(scale(rpd_middle)~scale(srs_sum)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(srs_sum)')
  ggplot(df_trial,aes(x=srs_sum,y=rpd_middle))+
    geom_smooth(method='lm')+
    geom_jitter(aes(color=timepoint))+theme_bw()
  ##--> higher pupillary response lower SRS

  ggplot(df_ids,
         aes(srs_sum,rpd_middle))+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()

  lmm<-lmer(scale(rpd_middle)~scale(RBSR_ges)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(RBSR_ges)')

  lmm<-lmer(scale(rpd_middle)~scale(CBCL_T_GES)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(CBCL_T_GES)')
  ggplot(df_trial,aes(x=CBCL_T_GES,y=rpd_middle))+
    geom_smooth(method='lm')+
    geom_jitter(aes(color=timepoint))+theme_bw()
  ##--> higher pupillary response lower CBCL


  #BASELINE PD
  lmm<-lmer(scale(baseline_pd)~scale(srs_sum)+(1|pic),data=df_trial)
  anova(lmm)

  lmm<-lmer(scale(baseline_pd)~scale(RBSR_ges)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(RBSR_ges)')

  lmm<-lmer(scale(baseline_pd)~scale(CBCL_T_GES)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(CBCL_T_GES)')
  ggplot(df_trial,aes(x=CBCL_T_GES,y=gazeduration_head))+geom_smooth(method='lm')+geom_jitter(aes(color=timepoint))+theme_bw()


  ##--> assoication of RJA with pupillary response ####
    with(df_trial,hist(baseline_pd)) ### mean pd in first 500ms of trial
    with(df_trial,hist(rpd_early,30)) ###change in mm compared to baseline
    with(df_trial,hist(rpd_middle,30)) ###change in mm compared to baseline
    with(df_trial,hist(rpd_late,30)) ###change in mm compared to baseline


  #prediction of rja by pupillary response
  lmm<-glmer(rja_any~scale(baseline_pd)+condition+age_months+
               (1|pic),data=df_trial,family=binomial(link = "logit"))
  summary(lmm)

  round(anova(lmm),3) %>%
    kbl(caption = "RJA likelihood") %>%
    kable_classic(full_width = F, html_font = "Cambria")

  fixef(lmm)['scale(baseline_pd)']
  confint(lmm,parm='scale(baseline_pd)')
  #-->a higher baseline pd is associated with less RJA
  ggplot(df_ids,aes(pd_500,rja_any_p))+geom_point()+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+theme_bw()



  lmm<-glmer(rja_any~scale(rpd_middle)+age_months_at_date+
               (1|pic),data=df_trial,family=binomial(link = "logit"))
  summary(lmm)

  round(anova(lmm),3) %>%
    kbl(caption = "RJA likelihood") %>%
    kable_classic(full_width = F, html_font = "Cambria")


  fixef(lmm)['scale(rpd_middle)']
  confint(lmm,parm='scale(rpd_middle)')
  #--> higher rps associated with higher likelihood of RJA (logOdds)
  ggplot(df_ids,aes(rpd_middle,rja_any_p))+geom_point()+geom_point(aes(color=group))+geom_smooth(color='grey',method='lm')+
    theme_bw()+
    coord_cartesian(xlim=c(-0.2,0.5))+
    labs(x='pupillary response (z)',y='RJA likelihood')


  lmm<-glmer(rja_any~scale(rpd_early)+timepoint+condition+
               (1|pic),data=df_trial,family=binomial(link = "logit"))
  anova(lmm)

  lmm<-glmer(rja_any~scale(rpd_late)+timepoint+condition+
               (1|pic),data=df_trial,family=binomial(link = "logit"))
  summary(lmm)
  ###--> early and late not assoicated with RJA likelihood

        lmm<-glmer(rja_any_early~scale(baseline_pd)+timepoint+
                     (1|pic)+(1|index_trial),data=df_trial,family=binomial(link = "logit"))
        summary(lmm)

        lmm<-glmer(rja_any_late~scale(baseline_pd)+timepoint+
                     (1|pic)+(1|index_trial),data=df_trial,family=binomial(link = "logit"))
        summary(lmm)
        ####--> higher baseline pd is negatively assoicated with premature rja and cued rja

        lmm<-glmer(rja_any_early~scale(rpd_middle)+timepoint+
                     (1|pic),data=df_trial,family=binomial(link = "logit"))
        summary(lmm)

        lmm<-glmer(rja_any_late~scale(rpd_middle)+timepoint+
                     (1|pic),data=df_trial,family=binomial(link = "logit"))
        summary(lmm)
        fixef(lmm)['scale(rpd_middle)']
        ####--> higher rpd  pd is positively assoicated with cued rja and not premature rja


  lmm<-glmer(rja_any~scale(rpd_middle)+timepoint+scale(no_data)+(1|pic),data=df_trial[!is.na(df_trial$rpd_middle),],family=binomial(link = "logit"))
  lmm2<-glmer(rja_any~timepoint+condition+scale(no_data)+(1|pic),data=df_trial[!is.na(df_trial$rpd_middle),],family=binomial(link = "logit"))
  anova(lmm,lmm2)
  ###--> lower baseline PD and higher middle rpd are associated with higher likelihood of RJA
  ####-> (especially a deliberate late rja)

  lmm<-lmer(scale(rja_duration)~scale(baseline_pd)*timepoint+condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)

  lmm<-lmer(scale(rja_duration_late)~scale(baseline_pd)*timepoint+condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)

  lmm<-lmer(scale(rja_duration_early)~scale(baseline_pd)*timepoint+condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)



  lmm<-lmer(scale(gazeduration_head)~scale(baseline_pd)*timepoint*condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)
  fixef(lmm)['scale(baseline_pd)']

  lmm<-lmer(scale(gazeduration_head)~scale(rpd_middle)*timepoint*condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)
  ggplot(df_trial,aes(x=rpd_middle,y=gazeduration_head))+geom_jitter(alpha=0.5,aes(color=timepoint))+geom_smooth(method='lm')+coord_fixed(ratio = 1/4)+theme_bw()


  ###--> association of head fixation with rja when and target duration ####
  lmm<-lmer(scale(rja_t_mean)~scale(gazeduration_headbefore)*timepoint*condition+
              scale(no_data)+(1|pic),data=df_trial)
  anova(lmm)
  confint(lmm,parm='scale(gazeduration_headbefore)')
  ggplot(df_trial,aes(x=gazeduration_headbefore,y=rja_t_mean,color=timepoint))+geom_jitter()+geom_smooth(method='lm')
  ###--> gaze on head duration predicts faster time on target

###-------- PRELIMINARY FINDINGS ####

  #why is there so many cases with NA in no_data? --> cases that have no gaze data

  #-->use rja estimation principle to code different Gaze shifts

  #FINDINGS:

    #main RJA FINDING:
    # T6 ist associated with higher likelihood of RJA and longer gaze duration on head compared to T2/T4
      #--> this is driven by more RJA that occur after cueing onset

    #PUPILLARY RESPONSE DIFFERENCES
    # margianlly higher pupillary response in K compared to T2
      #--> this is more pronounced in polynomial model with gradual effect K>T6>T4>T2

    # quantitative associations across groups:
    # RJA, gazeduration on head, and middle rpd are associated with quantitative symptoms (higher indices with lower scores)

    #data quality as covariate in RJA models, but not in RPD models - data quality is only associated with RJA likelihood

    #pupillary response association
    #higher baseline PD snd lower rpd middle are associated with less likely RJA
      #---> but not early rpd
        ####--> higher baseline pd is negatively assoicated with premature rja and cued rja
        ####--> higher rpd middle is positively assoicated with cued rja and not premature rja



    # faster RJA if longer head fixation before joint attention initiation

    # differential pupillary response in timepoints - in T2 a higher pupilalry response associated with faster RJA
    ###--> lower baseline PD and higher middle rpd are associated with higher likelihood of RJA
    ####-> (especially a deliberate late rja)


