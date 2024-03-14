# HEADER ####

## Script purpose: PREPROCESS BOSCA DATA - loads data and merges events with eye-tracking data
##
##
## Author: Nico Bast
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Nico Bast, `r paste(format(Sys.Date(), "%Y"))`
## Email: nico.bast@kgu.de
##
##  INFO: extracted stimulus onset info from splitting videos intro frames with custom python script
##  e.g.: CMD line = C:\Users\Nico\PowerFolders\project_video_salience\code\python_code_salience_extraction\split_video_to_images.py "neutral_rabbit_left" "C:\Users\Nico\PowerFolders\ETbat\Jointattention\stimuli\blockB"
##  # tested with four videos: length ~ 550 frames = 11 s, face onset ~ 64-72 frames ~ 1.35 sec, cueing onset ~ 165-180 frames ~  3.46 sec, stimulus_shaking ~ 150 frames = 3 sec
##
##  SYSTEM REQUIREMENTS: at least 64GB of RAM to read data of 314 sessions (13.03.2024)
##
##  line 1-170 is independent of the task and describes generel preprocessing of BOSCA data


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


#--> now found on NAS

end_time <- Sys.time()
end_time - start_time

# --> ##load("F:\temp_data_AFFIP\all_data_merged_040722.Rdata")
