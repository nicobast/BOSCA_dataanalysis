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

  ## B) delete data around blinks ####
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

### limit data to first 11 seconds ####

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


### MAIN OUTCOME: DEFINE REACTIVE JOINT ATTENTION (RJA) ####

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


###--> create TRIAL AGGREGATED DATA set (df_trial) #####

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


###--> save final data frame ####

save(df_jointatt, df_trial, demogr, file=paste0(home_path,project_path,"/data/all_data_preprocessed_FINAL_jointatt_140324.Rdata"))

