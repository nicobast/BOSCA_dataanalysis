# HEADER ####

## Script purpose: PREPROCESS ODDBALL DATA - BOSCA PROJECT
##
##    #only applies to voddball (prior ETbat07)
##
## Author: Nico Bast
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
## Email: nico.bast@kgu.de

sessionInfo()

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

# Load Data ####
start_time <- Sys.time()

datapath<-paste0(home_path,data_path) #preprocessed on LINUX machine

datapath<-"G:/BACKUP_Polzer_ETBattery/data_AFFIP/data"

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

#data.files<-data.files[c(grep('t2_gazedata.csv',data.files),grep('K_gazedata.csv',data.files))] #exclude other test points
#data.time<-data.time[c(grep('t2_timestamps.csv',data.time),grep('K_timestamps.csv',data.time))] #exclude other test points
#data.events<-data.events[c(grep('t2_event.csv',data.events),grep('K_event.csv',data.events))] #exclude other test points
# -> all three file lists (gazedata, time, events) should have the same length

#remove double entry of 094_t6
data.events<-data.events[-134]

# #read all gaze data into list class
# df.list.data<-lapply(data.files,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
# df.list.time<-lapply(data.time,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
# df.list.events<-lapply(data.events,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)

#read csv is very slow
test_list<-list(0)
for(i in 1:length(data.files)){
  test_list[[i]]<-fread(data.files[i])
  print(paste0('read: ',i))
}
###16 gb fails after 118 entries --> requires at least 32GB (June 2022)
### for loop is slower, but does not kill memory

df.list.data<-test_list
#df.list.data<-lapply(data.files,read.csv)
df.list.time<-lapply(data.time,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)
df.list.events<-lapply(data.events,read.csv,header = F, sep=",", dec=".", stringsAsFactors = T)

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
id.names<-substr(data.files,nchar(datapath)+2,nchar(datapath)+7) #08.07.19: now independent of relative to datapath
names(df)<-id.names

## --> SAVE MERGED DATA ####

#remove lists: clear ram
rm(df.list.data, df.list.events, df.list.time)

save(df,file="F:/temp_data_AFFIP/all_data_merged_040722.Rdata")

end_time <- Sys.time()
end_time - start_time

# --> ##load("F:\temp_data_AFFIP\all_data_merged_040722.Rdata")

#---------------------------------------------------------------------------------------#
## ----------------------- #####

## CREATE TIMESTAMP VARIABLE ####

#create long format variable: ts (in seconds format)
ts<-lapply(df,function(x){x<-x$timestamp})
ts<-lapply(ts,function(x){abs(head(x,n=1)-x)/1000000})

#add ts and change name of the variable
df<-mapply(cbind,df,ts,SIMPLIFY = FALSE) #simplify needs to be in capital letters
df<-lapply(df,fun_rename,variable_position=29,new_name='ts')

## --> SELECT TASK DATA: ODDBALL ####

  # #diagnostics - find right eventlogs
  # table(df[[100]]['eventlog'])
  # table(df[[220]]['eventlog'])
  # list_eventlog<-lapply(df,function(x){x['eventlog']})
  # df_eventlog<-dplyr::bind_rows(list_eventlog)
  # table(droplevels(df_eventlog$eventlog[grepl('voddball_exptrials_start',df_eventlog$eventlog)]))
  # table(droplevels(df_eventlog$eventlog))

  onset_voddball<-sapply(df,function(x){x$ts[which(x$eventlog=='voddball_exptrials_start')[1]]})
  offset_voddball<-sapply(df,function(x){x$ts[which(x$eventlog=='voddball_exptrials_end')[1]]})

  list_voddball<-mapply(function(x,y,z){x[x$ts>=y & x$ts<z,]},x=df,y=onset_voddball,z=offset_voddball, SIMPLIFY=F)

      # list_eventlog<-lapply(list_voddball,function(x){x['eventlog']})
      # df_eventlog<-dplyr::bind_rows(list_eventlog)
      # table(droplevels(df_eventlog$eventlog))

## --> CREATE VARIABLES: phase, trial, ts_event  ####

start_time <- Sys.time()

### - PHASE
### create a phase variable (target, fixcross, intro, outro)
fun_phase<-function(x){ifelse(grepl('attmovie',x[,'eventlog']),'attmovie',
                              ifelse(grepl('tar_',x[,'eventlog']),'target',
                                     ifelse(grepl('nov_',x[,'eventlog']),'novel',
                                            ifelse(grepl('dis_',x[,'eventlog']),'distractor',
                                                   ifelse(grepl('colscreen_col:1_dur:0.2',x[,'eventlog']),'fixcross',
                                                      ifelse(grepl('fixationcross_col:0_dur:0.2',x[,'eventlog']),'fixcross',
                                                            ifelse(grepl('colscreen_col:1_dur:2',x[,'eventlog']),'baseline','NA')))))))}


phase_data<-sapply(list_voddball,fun_phase)

#phase_data<-sapply(phase_data,as.factor)
list_voddball<-mapply(data.frame,list_voddball,phase_data, SIMPLIFY=F) #add
list_voddball<-lapply(list_voddball,fun_rename,variable_position=30,new_name='trial_phase') #rename


### - TRIAL + TS_EVENT

#identify trials --> function that based on change to identify trials and add trial sequence variable
fun_define_trials<-function(block_data){

  #create index variable - indicates when event changes #
  event.change<-lapply(block_data,function(x){which(diff(as.numeric(as.factor(x$trial_phase)))!=0)}) #index event change per participant
  length.datasets<-sapply(block_data,nrow) #length of each data set (per participant)
  event.change<-mapply(function(x,y){x<-c(0,x,y)},event.change,length.datasets) #add start and end event to index
  event.change<-lapply(event.change,diff) #create a difference value

  #create index for these new trials
  index.trial<-lapply(event.change,function(x){rep(seq_along(x),times=x)})

  #sequence over each new trial
  ts.event<-lapply(index.trial,function(x){as.numeric(do.call(c,by(x,x,seq_along, simplify=F)))}) #tts over each trial

  #add index and sequence data to list
  block_data<-mapply(function(x,y,z){data.frame(x,y,z)},block_data,index.trial,ts.event,SIMPLIFY = F)
  block_data<-lapply(block_data,fun_rename,variable_position=31,new_name='index_trial')
  block_data<-lapply(block_data,fun_rename,variable_position=32,new_name='ts_event')
  return(block_data)
}

#apply functions that add trial and trial sequence data
list_voddball<-fun_define_trials(list_voddball)

#---------------------------------------------------------------------------------------------------#
# DATA PREPROCESSING ####

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
list_voddball<-lapply(list_voddball,fun_screen_att)

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
list_voddball<-lapply(list_voddball,fun_screen_dist)

## -- drop unecessary data + particiapnts without data --> subsequent preprocessing requires far less RAM####
fun_required_necessary_data<-function(x){

  #drop raw eye tracking data
  x<-x[,!(grepl('gazepos2D',names(x)) | grepl('gazepos3D',names(x)) | grepl('eyepos',names(x)))]
  return(x)

}
list_voddball<-lapply(list_voddball,fun_required_necessary_data)

## -- remove participants without trial change - no data
dropped_participants<-which(sapply(list_voddball,function(x){max(x$index_trial,na.rm=T)})==1)
list_voddball<-list_voddball[-dropped_participants]

paste('date:',paste(Sys.Date()),',data of participants: n=',length(list_voddball))

# ## ----------> GAZE PREPROCESSING (Nyström, 2010 - CURRENTLY NOT IMPLEMENTED) ---------------------- ####
#
# start_time <- Sys.time()
#
# #variables for gaze preprocessing
# screen_width<-510 #mm on a 23 inch screen wiht 16:9 aspect ratio (Tobii TX 300 screen)
# screen_height<-290 #mm on a 23 inch screen wiht 16:9 aspect ratio (Tobii TX 300 screen)
# degrees_by_radian<-57.296 #fixed conversion facor
# velocity_cutoff<-1000 #visual degress per second
# acceleration_cutoff<-100000 #visual degress per second
# initial_velocity_cutoff<-200
# median_samples_trial<-1200
#
# #--> Savitksy Golay filter of length 15 ~ 50ms --> see for coefficients: http://www.statistics4u.info/fundstat_eng/cc_savgol_coeff.html
# #filter_sg15<-c(-78,-13,42,87,122,147,162,167,162,147,122,87,42,-13,-78)/1105
# filter_sg21<-c(-171,-76,9,84,149,204,249,284,309,324,329,324,309,284,249,204,149,84,9,-76,-171)/3059
#
# #required functions for gaze preprocessing:
#
# #A. - blink identification function - within 75ms to 250ms interval
# #--> INFO: identifies blinks and NAs 8 samples before and after it (~25ms)
# fun_blink_cor <- function(signal,lower_threshold=23,upper_threshold=75,samples_before=8,samples_after=8) {
#   #change NA to 999 for rle()-function
#   findna <- ifelse(is.na(signal),999,signal)
#   #find blinks:
#   #output of rle(): how many times values (NA) are repeated
#   repets <- rle(findna)
#   #stretch to length of PD vector for indexing
#   repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
#   #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
#   #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
#   #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
#   repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
#   #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
#   repets[findna!=999 & repets==1] <- 0
#   #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
#   changes <- c(diff(repets),0)
#   #define start (interval before blink/missing data)
#   changes.start<-which(changes==1) #where NA-sequence starts
#   #gives out row numbers of NA (blink) and previous 8 frames
#   start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
#   repets[start.seq]<-1
#   #define end (interval after blink/missing data)
#   changes.end<-which(changes==-1)+1 #where NA.sequence ends
#   #gives out row numbers of NA (blink) and subsequent 8 frames
#   end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
#   repets[end.seq]<-1
#   #replace PD data in blink interval (start to end) with NA
#   signal[repets==1]<-NA
#   return(signal)
# }
#
# #B. - return velocity
# fun_return_speed<-function(a,b,time){
#   time<-unlist(time)
#   gaze_diff_x<-diff(a)
#   gaze_diff_y<-diff(b)
#   gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
#   gaze_speed<-gaze_diff/diff(time)
#   gaze_speed<-c(NA,gaze_speed)
#   return(gaze_speed)
# }
#
# #C. - return acceleration
# fun_return_accel<-function(x,time){
#   time<-unlist(time)
#   diff_speed<-diff(x)
#   gaze_accel<-diff_speed/diff(time)
#   gaze_accel<-c(NA,gaze_accel)
#   return(gaze_accel)
# }
#
# #D. - function data driven velocity threshold to identify saccades --> see Nyström et al., 2010
# fun_velocity_threshold<-function(x){
#
#   cutoff<-initial_velocity_cutoff
#   diff_cutoff<-cutoff
#   k<-x
#
#   while(diff_cutoff>1){
#
#     k<-k[k<cutoff]
#     mean_speed<-median(k,na.rm=T)
#     sd_speed<-mean(k,na.rm=T)
#
#     cutoff_new<-mean_speed+(3*sd_speed)
#     if(is.na(cutoff_new)){
#       cutoff<-NA
#       break}
#
#     diff_cutoff<-cutoff-cutoff_new
#     if(cutoff_new>cutoff){cutoff_new<-cutoff}
#     cutoff<-cutoff_new
#
#   }
#
#   # #testing
#   # return(cutoff)
#
#   saccade_peak<-ifelse(x>cutoff,T,F)
#   return(saccade_peak)
#   #also return cutoffs in a list
#
# }
#
# #idnetififes most often values in a sequence
# fun_most_values <- function(x) {
#   ux <- unique(x)
#   tab <- tabulate(match(x, ux))
#   return(ux[tab == max(tab)])
#
# }
#
# #function gaze preprocessing (requires functions above)
# #--> inspired by Nyström et al. 2010 - denoising with Savitsky-Golay filter and adaptive velocity threshold filter
# fun_gaze_preprocess<-function(x){
#
#   #testing
#   x<-list_voddball[[10]]
#
#   #preprocess based on per trial level
#   x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
#   split_by_trial<-split(x,as.factor(x$index_trial))
#
#   #split_by_trial<-split(x,as.factor(x$index_trial))
#   timestamp<-sapply(split_by_trial,function(x){x<-x['ts']})
#   screen_distance<-sapply(split_by_trial,function(x){x<-x['screen_dist']})
#
#   #1. blink correction --> as noise
#   gazepos_x<-sapply(split_by_trial,function(x){fun_blink_cor(x$gazepos.x)})
#   gazepos_y<-sapply(split_by_trial,function(x){fun_blink_cor(x$gazepos.y)})
#
#   #3. drop trials with less than 50% of data
#   gazepos_x<-sapply(gazepos_x,function(x){if(sum(is.na(x))>0.5*median_samples_trial){x<-as.numeric(rep(NA,length(x)))}else{return(x)}})
#   gazepos_y<-sapply(gazepos_y,function(x){if(sum(is.na(x))>0.5*median_samples_trial){x<-as.numeric(rep(NA,length(x)))}else{return(x)}})
#
#   #2. convert relative gaze to degrees visual angle --> degrees from point of origin
#   gazepos_x_deg<-mapply(function(x,y){x<-x*atan(screen_width/y)*degrees_by_radian},x=gazepos_x,y=screen_distance,SIMPLIFY = F)
#   gazepos_y_deg<-mapply(function(x,y){x<-x*atan(screen_height/y)*degrees_by_radian},x=gazepos_y,y=screen_distance,SIMPLIFY = F)
#
#   #put together
#   unsplitting_factor<-as.factor(x$index_trial)
#   unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
#   gazepos_x<-unsplit(gazepos_x,f=unsplitting_factor)
#   gazepos_y<-unsplit(gazepos_y,f=unsplitting_factor)
#   gazepos_x_deg<-unsplit(gazepos_x_deg,f=unsplitting_factor)
#   gazepos_y_deg<-unsplit(gazepos_y_deg,f=unsplitting_factor)
#
#   x<-data.frame(unsplit_by_trials,gazepos_x,gazepos_y,gazepos_x_deg,gazepos_y_deg)
#   x$index_trial[x$index_trial==0]<-NA
#   return(x)
#
# }
#
# list_voddball<-lapply(list_voddball,fun_gaze_preprocess)
#
# #saccade identification based on velocity --> difficult with noisy data
# fun_saccade_ident<-function(x){
#
#   #x<-fun_gaze_preprocess(data_block1[[100]])
#
#   #preprocess based on per trial level
#   x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
#   split_by_trial<-split(x,as.factor(x$index_trial))
#   timestamp<-sapply(split_by_trial,function(x){x<-x['ts']})
#   gazepos_x_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_x_deg']})
#   gazepos_y_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_y_deg']})
#
#
#   # A. VELOCITY: FILTERING, DENOISING, PEAK IDENTIFICATION
#
#   #3. return velocity and acceleration variable --> in degrees per second
#   gaze_speed<-mapply(fun_return_speed,a=gazepos_x_deg,b=gazepos_y_deg,time=timestamp,SIMPLIFY = F)
#   gaze_accel<-mapply(fun_return_accel,x=gaze_speed,time=timestamp,SIMPLIFY = F)
#
#   #4. remove biologically implausible acceleration and velocity values
#   gaze_speed<-sapply(gaze_speed,function(x){ifelse(abs(x)>velocity_cutoff,NA,x)})
#   gaze_accel<-sapply(gaze_accel,function(x){ifelse(abs(x)>acceleration_cutoff,NA,x)})
#
#   #5. Denoising by Savitzky-Golay filter (length 21)
#   gaze_speed<-lapply(gaze_speed,function(x){
#     ifelse(length(x)<length(filter_sg21),
#            k<-x,
#            k<-as.numeric(stats::filter(x,filter_sg21)))
#     return(k)})
#
#   #--> SOME FORM OF DENOISING HAS TO BE APPLIED
#   #literature of signal denoise
#   #use a Savitzky Golay filter with length 15 ~ 50ms --> has been proven to be superior, see below
#
#
#   #filter/denoise benchmarking #
#   # sg7<-c(-2,3,6,7,6,3,-2)/21 #Savitsky Golay filter of length 7
#   # sg5<-c(-3,12,17,12,-3)/35 #Savitsky Golay filter of length 5
#   # sg21<-c(-171,-76,9,84,149,204,249,284,309,324,329,324,309,284,249,204,149,84,9,-76,-171)/3059
#   # sg15<-c(-78,-13,42,87,122,147,162,167,162,147,122,87,42,-13,-78)/1105
#   # ##--> see for coefficients: http://www.statistics4u.info/fundstat_eng/cc_savgol_coeff.html
#   #
#   # x<-gaze_speed[[15]]
#   # par(mfrow=c(3,2))
#   # plot(x,main='data')
#   # plot(as.numeric(filter(x,sg5)),main='S-G5 filter')
#   # plot(as.numeric(filter(x,sg15)),main='S-G15 filter')
#   # plot(as.numeric(filter(x,sg21)),main='S-G21 filter')
#   # #plot((na.approx(x,maxgap = interpolate_cutoff)),main='linear interpol')
#   # plot(pracma::savgol(na.approx(x),fl=21),main='S-G long')
#   # plot(spline(x,n=length(x)),main='spline')
#   # par(mfrow=c(1,1))
#   # ###--> compare different denoise / filtering methods
#   # ##--> Savitzky-Golay filter seems to be the winner - use a length of 21
#
#   #6.data-driven velocity threshold and identify peaks (saccades) based on it
#   velocity_peak<-lapply(gaze_speed,fun_velocity_threshold)
#
#   # # #  # #testing
#   # plot(gaze_speed[[10]],col=as.factor(unlist(velocity_peak[[10]])))
#   # plot(unlist(gaze_speed),col=as.factor(unlist(velocity_peak)))
#
#   #B.: DISPERSION
#
#   # fun_dispersion_filter<-function(a,b){
#   #
#   #   gaze_diff_x<-diff(a)
#   #   gaze_diff_y<-diff(b)
#   #   gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
#   #   gaze_diff<-c(NA,gaze_diff)
#   # }
#   #
#   # gaze_diff<-mapply(fun_dispersion_filter,a=gazepos_x,b=gazepos_y,SIMPLIFY=F)
#
#   #C.: DURATION
#
#   #modus in a window of ten samples
#   velocity_continues<-frollapply(velocity_peak,n=10,fun_most_values,align='center')
#   velocity_continues<-lapply(velocity_continues,as.logical)
#
#   #D.: DEFINE SACCADE
#
#   saccade<-mapply(function(x,y){ifelse(x & y,T,F)},x=velocity_peak,y=velocity_continues,SIMPLIFY=F)
#
#   #return trial-wise list to data.frame
#   unsplitting_factor<-as.factor(x$index_trial)
#   unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
#   gaze_speed<-unsplit(gaze_speed,f=unsplitting_factor)
#   gaze_accel<-unsplit(gaze_accel,f=unsplitting_factor)
#   velocity_peak<-unsplit(velocity_peak,f=unsplitting_factor)
#   velocity_continues<-unsplit(velocity_continues,f=unsplitting_factor)
#   saccade<-unsplit(saccade,f=unsplitting_factor)
#
#
#   x<-data.frame(unsplit_by_trials,gaze_speed,gaze_accel,velocity_peak,velocity_continues,saccade)
#   x$index_trial[x$index_trial==0]<-NA
#   return(x)
# }
#
# list_voddball<-lapply(list_voddball,fun_saccade_ident)
#
# #fixation identification --> difficult with noisy data
# fun_fixation_ident<-function(x,degree_fixation_cutoff=1,duration_fixation_cutoff=30){
#
#   #testing
#   # x<-data_block1[[100]]
#   # x<-fun_gaze_preprocess(x)
#   # x<-fun_saccade_ident(x)
#
#   x$index_trial<-ifelse(is.na(x$index_trial),0,x$index_trial) #consider NA as trial=0
#   split_by_trial<-split(x,as.factor(x$index_trial))
#   gazepos_x_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_x_deg']})
#   gazepos_y_deg<-sapply(split_by_trial,function(x){x<-x['gazepos_y_deg']})
#
#
#   #gaze difference from sampel to sample
#   fun_gaze_diff<-function(a,b){
#
#     gaze_diff_x<-diff(a)
#     gaze_diff_y<-diff(b)
#     gaze_diff<-sqrt((gaze_diff_x^2)+(gaze_diff_y^2))
#     gaze_diff<-c(NA,gaze_diff)
#
#   }
#
#   gaze_diff<-mapply(fun_gaze_diff,a=gazepos_x_deg,b=gazepos_y_deg,SIMPLIFY=FALSE)
#
#   #identify significant movement in subsequent or preceding samples
#
#   fun_movement_ident<-function(x){
#
#     no_movement_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(mean(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='left')
#     no_movement_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(mean(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='right')
#     #no_movement_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(all(x<degree_fixation_cutoff),T,F)},align='left')
#     #no_movement_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(all(x<degree_fixation_cutoff),T,F)},align='right')
#     no_movement<-ifelse(no_movement_next_samples | no_movement_last_samples,T,F)
#
#   }
#
#   no_movement<-sapply(gaze_diff,fun_movement_ident)
#
#   #identify significant drifts in subsequent or preceding samples
#
#   fun_gaze_diff_abs<-function(x){
#
#     gaze_diff_abs<-diff(x)
#     gaze_diff_abs<-c(NA,gaze_diff_abs)
#
#   }
#
#   gaze_diff_abs_x<-sapply(gazepos_x_deg,fun_gaze_diff_abs)
#   gaze_diff_abs_y<-sapply(gazepos_y_deg,fun_gaze_diff_abs)
#
#   fun_drift_ident<-function(x){
#
#     no_drift_next_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(sum(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='left')
#     no_drift_last_samples<-frollapply(x,n=duration_fixation_cutoff,function(y){ifelse(sum(y,na.rm=T)<degree_fixation_cutoff,T,F)},align='right')
#     no_drift<-ifelse(no_drift_next_samples | no_drift_last_samples,T,F)
#
#   }
#
#   no_drift_x<-sapply(gaze_diff_abs_x,fun_drift_ident)
#   no_drift_y<-sapply(gaze_diff_abs_y,fun_drift_ident)
#
#
#   #--> define a fixation sample if no movement or drift in the last or next 100ms ~ 30 samples
#   fixation<-mapply(function(x,y,z){ifelse(x & y & z,T,F)},x=no_movement,y=no_drift_x,z=no_drift_y)
#
#   #add to data
#   unsplitting_factor<-as.factor(x$index_trial)
#   unsplit_by_trials<-unsplit(split_by_trial,f=unsplitting_factor)
#   fixation<-unsplit(fixation,f=unsplitting_factor)
#   x<-data.frame(unsplit_by_trials,fixation)
#   x$index_trial[x$index_trial==0]<-NA
#   return(x)
#
# }
#
# list_voddball<-lapply(list_voddball,fun_fixation_ident)
#
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
list_voddball<-lapply(list_voddball, func_pd_preprocess)

### - select relevant data ####

fun_select_data<-function(x){
  x<-x[,names(x) %in% c('eventlog','timestamp','ts','trial_phase','index_trial',
                        'ts_event','gazepos.x','gazepos.y',
                        'pd','center_deviation','screen_dist','gazepos_x','gazepos_y','gaze_speed','gaze_accel',
                        'velocity_peak','velocity_continues','saccade','fixation')]
  return(x)
}

list_voddball<-lapply(list_voddball,fun_select_data)

### - create ID variable - required before concatenating blocks ####
fun_retrieve_id<-function(x,id){
  k<-nrow(x)
  id<-rep(id,k)
  x[,'id']<-id
  return(x)
}

list_voddball<-mapply(fun_retrieve_id,x=list_voddball,id=names(list_voddball),SIMPLIFY=F)

### - calculate baseline pupil size - corrected pupil size ####
fun_baseline<-function(x){

  split_by_trial<-split(x,as.factor(x$index_trial))
  baseline_data<-sapply(split_by_trial,function(x){mean(x$pd[x$ts_event<=150],na.rm=T)}) #select between event
  rpd<-unlist(mapply(function(x,y){x$pd-y},x=split_by_trial,y=baseline_data)) #correct for baseline by subtraction
  baseline_pd<-rep(baseline_data,sapply(split_by_trial,nrow)) #calculate baseline pd

  x[,'rpd']<-rpd
  x[,'baseline_pd']<-baseline_pd
  return(x)

}

list_voddball<-lapply(list_voddball,fun_baseline)

###------ save preprocessed data (single sample level) #####
df_oddball<-dplyr::bind_rows(list_voddball)

save(df_oddball,file=paste0(home_path,project_path,"/data/all_data_preprocessed_voddball_050922.Rdata"))


### -- CREATE PIC; GROUP, and TIMEPOINT VARIABLE ####

#load(paste0(home_path,project_path,"/data/all_data_preprocessed_voddball_050922.Rdata"))

#create group variable
df_oddball$group<-ifelse(grepl('_K',df_oddball$id),'TD','ASD')

#create timepoint variable
df_oddball$timepoint<-ifelse(grepl('_t2',df_oddball$id) | grepl('_T2',df_oddball$id),'T2',
                     ifelse(grepl('_t4',df_oddball$id) | grepl('_T4',df_oddball$id),'T4',
                            ifelse(grepl('_t6',df_oddball$id) | grepl('_T6',df_oddball$id),'T6','K')))

#create individual id variable
df_oddball$pic<-substr(df_oddball$id,1,3)

#### --> DATA ANALYSIS  ####

hist(df_oddball$ts_event)
hist(df_oddball$pd)
hist(df_oddball$rpd)
table(df_oddball$trial_phase)

with(df_oddball,by(pd,trial_phase,summary))

require(ggplot2)
df_oddball_core<-df_oddball[df_oddball$trial_phase %in% c('distractor','novel','target'),]
ggplot(df_oddball_core[df_oddball_core$ts_event<750,],aes(x=ts_event,y=rpd,group=trial_phase,color=trial_phase))+geom_smooth()+theme_bw()
####--> does not work - stimulus onset produces PLR

