# setup ####
require(ggplot2)
require(readxl)
require(lme4)
require(lmerTest)
require(emmeans)
require(MuMIn)

#TO DO
 #-include demographics
 #-estimate pupillary response

#CONCLUSION:
  #--> cut down to 32 trials
  #--> make the blocks abortable

#detect system
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

# load preprocessed data (see code/data_preprocessing_0421.R) ####
load(file=paste0(home_path,'/PowerFolders/data_AFFIP/data_preprocessed_300421.Rdata'))
#load(file=paste0(home_path,'/PowerFolders/data_AFFIP/data_preprocessed_290622.Rdata'))
#note: pd is already preprocessed

df_demo<-read_xlsx(paste0(home_path,'/PowerFolders/data_AFFIP/demographic_data_1.xlsx'))

table(df$id)

#show diferent evens
sort(table(df$eventlog))
hist(df$pd)


# peprocess demographic data ####

names(df_demo)

table(df_demo$Geschlecht_Index)
sex<-as.factor(with(df_demo,ifelse(Geschlecht_Index=='weiblich','female','male')))


#calculate age
  table(is.na(df_demo$Datum_BayleyIII))
  table(is.na(df_demo$T1_WPPSI_3_0_3_11_Datum))
  table(is.na(df_demo$T1_WPPSI_4_0_7_2_Datum))

  #get date from different variables
  assessment_date<-with(df_demo,ifelse(!is.na(Datum_BayleyIII),Datum_BayleyIII,
                                       ifelse(!is.na(T1_WPPSI_3_0_3_11_Datum),T1_WPPSI_3_0_3_11_Datum,
                                                     ifelse(!is.na(T1_WPPSI_4_0_7_2_Datum),T1_WPPSI_4_0_7_2_Datum,NA))))

  #transfer attribute to retain date class:
  attributes(assessment_date)<-attributes(df_demo$Datum_BayleyIII)

  age<-with(df_demo,assessment_date-Geburt_Index)
  age<-ifelse(!is.na(age),as.numeric(age/365),df_demo$ADI_alter)

#cognitive ability

  table(is.na(df_demo$Kog_EQ))
  table(is.na(df_demo$T1_WPPSI_3_0_3_11_Handlungsteil_IQ))
  table(is.na(df_demo$T1_WPPSI_4_0_7_2_Handlungsteil_IQ))

  iq<-ifelse(!is.na(df_demo$Kog_EQ),df_demo$Kog_EQ,
             ifelse(!is.na(df_demo$T1_WPPSI_3_0_3_11_Handlungsteil_IQ),df_demo$T1_WPPSI_3_0_3_11_Handlungsteil_IQ,
                           ifelse(!is.na(df_demo$T1_WPPSI_4_0_7_2_Handlungsteil_IQ),df_demo$T1_WPPSI_4_0_7_2_Handlungsteil_IQ,NA)))

#ADOS_CSS

  ados_css<-with(df_demo,ifelse(!is.na(ADOS_2_M1_comp_Score),ADOS_2_M1_comp_Score,
                                ifelse(!is.na(ADOS_2_M2_comp_Score),ADOS_2_M2_comp_Score,NA)))

#ADI

  adi_si<-df_demo$ADI_alg_SOZ_INT
  adi_com<-df_demo$ADI_alg_KOM
  adi_rb<-df_demo$ADI_alg_RITUALE

#srs total

  srs_total<-with(df_demo,ifelse(!is.na(T1_SRS_LM_Gesamtwert_N_k_RW),T1_SRS_LM_Gesamtwert_N_k_RW,
                                 ifelse(!is.na(T1_SRS_LV_Gesamtwert_N_k_RW),T1_SRS_LV_Gesamtwert_N_k_RW,
                                        ifelse(!is.na(T1_SRS_LV_Gesamtwert_N_k_RW),T1_SRS_LV_Gesamtwert_N_k_RW,
                                               ifelse(!is.na(T2_SRS_LM_Gesamtwert_N_k_RW),T2_SRS_LM_Gesamtwert_N_k_RW,
                                                      ifelse(!is.na(T2_SRS_LV_Gesamtwert_N_k_RW),T2_SRS_LV_Gesamtwert_N_k_RW,
                                                             ifelse(!is.na(T2_SRS_LV_Gesamtwert_N_k_RW),T2_SRS_LV_Gesamtwert_N_k_RW,NA)))))))

#rbsr_r total

  rbsr_total<-with(df_demo,ifelse(!is.na(T1_RBS_LM_RBSR_Lam_Ges),T1_RBS_LM_RBSR_Lam_Ges,
                                 ifelse(!is.na(T1_RBS_LV_RBSR_Lam_Ges),T1_RBS_LV_RBSR_Lam_Ges,
                                        ifelse(!is.na(T1_RBS_LV_RBSR_Lam_Ges),T1_RBS_LV_RBSR_Lam_Ges,
                                               ifelse(!is.na(T2_RBS_LM_RBSR_Lam_Ges),T2_RBS_LM_RBSR_Lam_Ges,
                                                      ifelse(!is.na(T2_RBS_LV_RBSR_Lam_Ges),T2_RBS_LV_RBSR_Lam_Ges,
                                                             ifelse(!is.na(T2_RBS_LV_RBSR_Lam_Ges),T2_RBS_LV_RBSR_Lam_Ges,NA)))))))

  #create data frame (sample description)
  df_describe<-data.frame(df_demo$ID_Bado,df_demo$ID_Studie,df_demo$Gruppe,sex,age,iq,ados_css,adi_si,adi_com,adi_rb,srs_total,rbsr_total)
  names(df_describe)[1:3]<-c('id_bado','id_core','group')

  #streamline id variable (remove "_K" from id for controls)
  levels(df_describe$id_core)[which(nchar(levels(df_describe$id_core))>3)]<-substr(levels(df_describe$id_core)[which(nchar(levels(df_describe$id_core))>3)],1,3)
  levels(df_describe$group)<-c('ASD','TD')

  #table description
  with(df_describe,by(age,group,psych::describe))
  with(df_describe,t.test(age~group))

  with(df_describe,by(sex,group,table))

  with(df_describe,by(iq,group,psych::describe))
  with(df_describe,t.test(iq~group))

  psych::describe(df_describe$ados_css)

  with(df_describe,by(srs_total,group,psych::describe))
  with(df_describe,t.test(srs_total~group))

  with(df_describe,by(rbsr_total,group,psych::describe))
  with(df_describe,t.test(rbsr_total~group))

#------- TASK: VISUAL SEARCH ####
vissearch<-df[grep('target:',df$eventlog),]
hist(table(vissearch$id))
rm(df)

# onscreen time
vissearch_split<-split(vissearch,vissearch$id)
nobs<-sapply(vissearch_split,nrow) #number of observations
possible_nobs<- 469*16*4 # 16 trials per block * 4 blocks (469=number timestamps per stimulus "slide")
nobs_percent<-nobs/possible_nobs*100
round(nobs_percent,2)
hist(nobs_percent)
psych::describe((round(nobs_percent,2)))
rm(vissearch_split)

#percent of onscreen time by group
with(vissearch,by(nobs_percent,grepl('_t',names(nobs)),psych::describe))

## 1. CALCULATE VARIABLES ####

##  Hits on target

#target positions in visual search in (0,1) space
target.x<-(960+c(-375,-375,-125,-125,+125,+125,+375,+375))/1920
target.y<-(540+c(-125,+125,-375,+375,-375,+375,-125,+125))/1080

rect.size<-165 #taken from matlab
rect.size.x<-rect.size/1920
rect.size.y<-rect.size/1080
target.rect.xmin<-target.x-rect.size.x/2
target.rect.xmax<-target.x+rect.size.x/2
target.rect.ymin<-target.y-rect.size.y/2
target.rect.ymax<-target.y+rect.size.y/2

#vissearch - target.position factor variable
#target.pos<-data.frame(target.x,target.y)
target.x.fac<-as.factor(substr(vissearch$eventlog,9,13))
levels(target.x.fac)<-c('1235','485','735','985')
target.y.fac<-as.factor(substr(vissearch$eventlog,14,18))
levels(target.y.fac)<-c('315','565','315','565','65','65','815','815') #more levels are necessary as string with '1235' changes length of string
target.fac<-interaction(target.x.fac,target.y.fac)
target.fac<-droplevels(target.fac)

#HITS in vissearch
target.hits.var<-rep(0,nrow(vissearch))
target.list<-c('_485_315_','_485_565_','_735_65_','_735_815_','_985_65_','_985_815_','_1235_315_','_1235_565_')
for(i in 1:length(target.list)){
  #rows with hit of target 1
  target.hits<-which((vissearch$gazepos.x>target.rect.xmin[i])&(vissearch$gazepos.x<target.rect.xmax[i])&(vissearch$gazepos.y>target.rect.ymin[i])&(vissearch$gazepos.y<target.rect.ymax[i]))
  #rows with aim of target 1
  target.trials<-grep(target.list[i],vissearch$eventlog)
  #rows with hit of target 1 with aim of target 1
  target.hits<-target.hits[target.hits %in% target.trials]
  target.hits.var[target.hits]<-1}

target_hit<-target.hits.var

vissearch <- data.frame(vissearch, target_hit)

    # numer of hits
    table(vissearch$target_hit)
    #hists by participant
    hist(with(vissearch,by(target_hit,id,function(x){table(x)[2]/length(x)})))

    par(mfrow=c(1,2))
    hist(with(vissearch[vissearch$group=='ASD',],by(target_hit,id,function(x){table(x)[2]/length(x)})),10,main='ASD')
    hist(with(vissearch[vissearch$group=='TD',],by(target_hit,id,function(x){table(x)[2]/length(x)})),10,main='TD')

    with(vissearch,by(ts.event,droplevels(interaction(id,index.event)),max,na.rm=T))


# hits on stimuli locations/areas (independent from target/distractor category)
locations <- rep(0, 8)
for(i in 1:8) {
  locations[i] <- length(which(vissearch$gazepos.x>target.rect.xmin[i] & vissearch$gazepos.x<target.rect.xmax[i] & vissearch$gazepos.y>target.rect.ymin[i] & vissearch$gazepos.y<target.rect.ymax[i]))
}
location.hits <- sum(locations)
possible.location.hits <- nrow(vissearch)
round(location.hits/possible.location.hits*100,2)  # percentage of gazeposition in stimulus locations

##gaze center deviation
center_deviation<-with(vissearch,sqrt((gazepos.x-0.5)^2+(gazepos.y-0.5)^2))
hist(center_deviation)

## 2. CREATE data frame per trial ####

#aggregate by trial
int_trial<-with(vissearch,interaction(id,eventlog,index.event)) #eventlog needed as eventlogs occure two times per id

factors_trial<-with(vissearch,aggregate(data.frame(group,id,index.event),
                                       by=list(int_trial),head,n=1))

metrics_mean_trial<-with(vissearch,aggregate(data.frame(pd,ts,ts.event,target_hit,gazepos.x,gazepos.y,center_deviation),
                                            by=list(int_trial),mean,na.rm=T))

table(df_vissearch_trial$target)

#center fixation at start

cent_dev_start<-as.numeric(by(center_deviation[vissearch$ts.event<50],
                              droplevels(int_trial)[vissearch$ts.event<50],
                              mean,na.rm=T))

#hit latency (seconds to hit)
  #TODO: alternative calculation: respective to a moving avergae, if above threshold
hit_latency<-as.numeric(by(vissearch$ts.event[vissearch$target_hit==1],
                          droplevels(int_trial)[vissearch$target_hit==1],
                          function (x) {min(x, na.rm=T)*0.0033}))

  #plot(density(hit_latency,na.rm=T))

###relative pupil size

#ggplot(vissearch,aes(x=ts.event,y=scale(pd),group=group,color=group))+geom_smooth()

pd_peak<-with(vissearch,as.numeric(by(pd[ts.event>125 & ts.event<=175],
                                      droplevels(int_trial)[ts.event>125 & ts.event<=175],
                                      mean,na.rm=T)))

pd_baseline<-with(vissearch,as.numeric(by(pd[ts.event<=50],
                                          droplevels(int_trial)[ts.event<=50],
                                          mean,na.rm=T)))

pd_response<-pd_peak-pd_baseline

#data quality variable
recorded_samples_per_trial<-table(as.numeric(droplevels(int_trial)))
max_recorded_samples<-469
valid_data<-as.numeric(recorded_samples_per_trial/max_recorded_samples)

#concatenate to data.frame
df_vissearch_trial<-data.frame(factors_trial,metrics_mean_trial[,-1],valid_data,hit_latency,pd_response,cent_dev_start,pd_baseline)
df_vissearch_trial$valid_data_cat<-ifelse(valid_data>0.5,T,F)

#create additional variables
df_vissearch_trial<-df_vissearch_trial[order(df_vissearch_trial$id,df_vissearch_trial$index.event),] #sort order for calculating trial number
#trial_number<-as.numeric(unlist(with(df_vissearch_trial,by(index.event,id,seq_along))))

    #create trial number
          #problem: intermediate trials might be missing, thus calculate change value and use cumsum to create trial number
           trial_change<-as.numeric(with(df_vissearch_trial,round(diff(index.event)/11))) #differences in index.event indicate trial change
           trial_change<-c(1,trial_change)
           trial_change[trial_change<0]<-1 # <0 == change in participants (id) in sequence
           trial_change<-ifelse(trial_change>5,1,trial_change) #control for changes in blocks
    trial_number<-as.numeric(unlist(with(df_vissearch_trial,by(trial_change,id,cumsum))))

    #extract condition variables
    target<-as.factor(substr(df_vissearch_trial$Group.1,14,14))
    target_location<-as.factor(substr(df_vissearch_trial$Group.1,16,23))
      #correct levels:
      levels(target_location)<-ifelse(substr(levels(target_location),8,8)!='5',substr(levels(target_location),1,7),levels(target_location))
      levels(target_location)<-ifelse(substr(levels(target_location),nchar(levels(target_location)),nchar(levels(target_location)))!='5',substr(levels(target_location),1,6),levels(target_location))


## 3. merge data.frame with df_describe ####
 df_vissearch_trial<-data.frame(df_vissearch_trial,trial_number,target,target_location)

 par(mfrow=c(1,1))
 hist(df_vissearch_trial$target_hit)

 #merge data frame with demographics

 df_vissearch_trial<-df_vissearch_trial[,-grep('group',names(df_vissearch_trial))] #remove group as also in df_describe
 df_vissearch_trial$id_core<-substr(df_vissearch_trial$id,1,3) #allign id variable
 df_vissearch_trial<-merge(df_vissearch_trial,df_describe,by='id_core')

###--> PRELIMINARY RESULTS (unaggregated): ####

 ###Gliga 2015 was not controlled for cognitive ability or data qualitiy -->
 ###Gliga 2015: "proportion of trials on target"

 ##association of pd_baseline and pd_response

 ggplot(df_vissearch_trial,aes(x=scale(pd_baseline),y=scale(pd_response),color=group))+geom_smooth()+geom_point(alpha=0.1)+xlim(c(-2,2))+ylim(c(-2,2))+theme_bw()

 ggplot(df_vissearch_trial[df_vissearch_trial$trial_number<65,],aes(x=trial_number,y=scale(pd_response)))+geom_smooth()

 ggplot(vissearch,aes(x=ts.event,y=scale(pd),group=group,color=group))+geom_smooth()


 ## define valid trials
 df_vissearch_trial_valid<-df_vissearch_trial[df_vissearch_trial$cent_dev_start<0.2 &
                                                df_vissearch_trial$valid_data>0.5,]

 ##estimate random effect structure

 re_model1<-lmer(target_hit~1+
                   (1|id_core),
                 data = df_vissearch_trial_valid,REML=F)

 re_model2<-lmer(target_hit~1+
                   (1|target_location)+(1|id_core),
                 data = df_vissearch_trial_valid,REML=F)

 re_model3<-lmer(target_hit~1+
                   (1|target)+(1|target_location)+(1|id_core),
                 data = df_vissearch_trial_valid,REML=F)

 re_model4<-lmer(target_hit~1+
                   (1|trial_number)+(1|target)+(1|target_location)+(1|trial_number)+(1|id_core),
                 data = df_vissearch_trial_valid,REML=F)


 anova(re_model1,re_model2,re_model3,re_model4)

 #---> trial number has no systematic effect on hit likelihood

 summary(lmer(target_hit~trial_number+
                (1|target)+(1|target_location)+(1|id_core),
              data = df_vissearch_trial_valid))

 ggplot(df_vissearch_trial[df_vissearch_trial$trial_number<65,],aes(x=trial_number,y=target_hit))+geom_smooth()

 #---> trial number has effect on hit latency
 summary(lmer(scale(hit_latency)~trial_number+
                (1|target)+(1|target_location)+(1|id_core),
              data = df_vissearch_trial_valid))

 ggplot(df_vissearch_trial[df_vissearch_trial$trial_number<65,],aes(x=trial_number,y=hit_latency))+geom_smooth()


 trials_per_id<-with(df_vissearch_trial,as.numeric(by(trial_number,id,function(x){length(unique(x))})))

 hist(trials_per_id)

 #do task conditions affect results
 df_vissearch_trial$hit_cat<-with(df_vissearch_trial,ifelse(target_hit<0.065,F,T)) #translates to 100ms
 summary(glmer(hit_cat~target+target_location+
                 (1|trial_number)+(1|id),family=binomial,
               data=df_vissearch_trial[df_vissearch_trial$valid_data_cat==T & df_vissearch_trial$trial_number<=32,]))

 #differences in hits between groups
    #RBSR associated with larger BPS and attenuated pupillar response across groups
    #In TD, a higher pupillary response is associated with faster hits on target (lower latency)

       lmm_model<-lmer(scale(pd_baseline)~scale(rbsr_total)+group+scale(valid_data)+
                         scale(age)+sex+scale(iq)+
                         (1|target)+(1|target_location)+(1|id),
                       data=df_vissearch_trial_valid)

       anova(lmm_model)
       cbind(fixef(lmm_model)['scale(rbsr_total)'],
             confint(lmm_model,parm='scale(rbsr_total)'))


       lmm_model<-lmer(scale(pd_response)~scale(rbsr_total)+group+scale(valid_data)+
                         scale(age)+sex+scale(iq)+
                         (1|target)+(1|target_location)+(1|id),
                       data=df_vissearch_trial_valid)

       anova(lmm_model)
       cbind(fixef(lmm_model)['scale(rbsr_total)'],
             confint(lmm_model,parm='scale(rbsr_total)'))

       lmm_model<-lmer(scale(hit_latency)~scale(rbsr_total)+group*scale(pd_response)+scale(valid_data)+
                         scale(age)+sex+scale(iq)+
                       (1|target)+(1|target_location)+(1|id),
                     data=df_vissearch_trial_valid)

              anova(lmm_model)
              require(emmeans)
              fixef(lmm_model)
              emtrends(lmm_model,~group,var='pd_response')

        lmm_model2<-lmer(scale(target_hit)~scale(rbsr_total)+group*scale(pd_response)+scale(valid_data)+
                           scale(age)+sex+scale(iq)+
                       (1|target)+(1|target_location)+(1|id),
                      data=df_vissearch_trial_valid)

        anova(lmm_model2)
        fixef(lmm_model2)

        df_vissearch_trial_valid$hit_cat<-with(df_vissearch_trial_valid,ifelse(df_vissearch_trial_valid$target_hit<0.065,F,T)) #translates to 100ms
        glmm_model<-glmer(hit_cat~group*scale(pd_response)+scale(rbsr_total)+scale(valid_data)+
                            scale(age)+sex+scale(iq)+
                        (1|target)+(1|target_location)+(1|id),family=binomial,
               data=df_vissearch_trial_valid)

        summary(glmm_model)
        emmeans(glmm_model,~group)

        table(df_demo$Gruppe)

        #--> findings not very different when cut by half
        summary(lmer(scale(hit_latency)~group+scale(valid_data)+
                       (1|target)+(1|target_location)+(1|id),
                     data=df_vissearch_trial[df_vissearch_trial$trial_number<=32,]))

        summary(lmer(scale(target_hit)~group+scale(valid_data)+
                       (1|target)+(1|target_location)+(1|id),
                     data=df_vissearch_trial[df_vissearch_trial$trial_number<=32,]))

        df_vissearch_trial$hit_cat<-with(df_vissearch_trial,ifelse(target_hit<0.065,F,T)) #translates to 100ms
        summary(glmer(hit_cat~group+scale(valid_data)+
                        (1|trial_number)+(1|id),family=binomial,
                      data=df_vissearch_trial[df_vissearch_trial$trial_number<=32,]))


## 4. aggregate to per_id data ####

        factors_id<-with(df_vissearch_trial,aggregate(data.frame(id_core,id_bado,group,sex),
                                                by=list(id_core),head,n=1))

        metrics_mean_id<-with(df_vissearch_trial,aggregate(data.frame(pd,valid_data,cent_dev_start,gazepos.x,gazepos.y,hit_latency,pd_response,pd_baseline,
                                                                      valid_data_cat,hit_cat,
                                                                      age,iq,
                                                                      ados_css,adi_com,adi_si,adi_rb,srs_total,rbsr_total),
                                                     by=list(id_core),mean,na.rm=T))

        df_vissearch_id<-data.frame(factors_id,metrics_mean_id)


###--> preliminary analysis (aggregated data) ####

        with(df_vissearch_id,by(pd_baseline,group,psych::describe))
        with(df_vissearch_id,by(pd_response,group,psych::describe))
        with(df_vissearch_id,by(hit_latency,group,psych::describe))
        with(df_vissearch_id,by(hit_cat,group,psych::describe))

        summary(lm(pd_response~group+age+iq+sex,df_vissearch_id))
        summary(lm(pd_baseline~group+age+iq+sex,df_vissearch_id))
        summary(lm(hit_latency~group+age+iq+sex,df_vissearch_id))
        summary(lm(hit_cat~group+age+iq+sex,df_vissearch_id))

        summary(lm(pd_response~group+age+iq+sex,df_vissearch_id))
        summary(lm(pd_baseline~group+age+iq+sex,df_vissearch_id))
        summary(lm(hit_latency~group+age+iq+sex,df_vissearch_id))
        summary(lm(hit_cat~group+age+iq+sex,df_vissearch_id))

        #associations of repetitive behavior across groups with pupillary metrics
        summary(lm(scale(pd_response)~scale(srs_total)+scale(rbsr_total)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(pd_baseline)~scale(srs_total)+scale(rbsr_total)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(hit_latency)~scale(srs_total)+scale(rbsr_total)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(hit_cat)~scale(srs_total)+scale(rbsr_total)+age+iq+sex,df_vissearch_id))

        summary(lm(scale(hit_latency)~scale(pd_response)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(hit_cat)~scale(pd_response)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(hit_latency)~scale(pd_baseline)+age+iq+sex,df_vissearch_id))
        summary(lm(scale(hit_cat)~scale(pd_baseline)+age+iq+sex,df_vissearch_id))

        #in ASD group
        df_vissearch_id_ASD<-df_vissearch_id[df_vissearch_id$group=='ASD',]

        summary(lm(scale(hit_latency)~scale(pd_response)+age+iq+sex,df_vissearch_id_ASD))
        summary(lm(scale(hit_cat)~scale(pd_response)+age+iq+sex,df_vissearch_id_ASD))
        summary(lm(scale(hit_latency)~scale(pd_baseline)+age+iq+sex,df_vissearch_id_ASD))
        summary(lm(scale(hit_cat)~scale(pd_baseline)+age+iq+sex,df_vissearch_id_ASD))

        summary(lm(scale(pd_response)~scale(ados_css)+scale(adi_si)+scale(adi_com)+scale(adi_rb)+age+sex,df_vissearch_id_ASD))
        summary(lm(scale(pd_baseline)~scale(ados_css)+scale(adi_si)+scale(adi_com)+scale(adi_rb)+age+sex,df_vissearch_id_ASD))
        summary(lm(scale(hit_latency)~scale(ados_css)+scale(adi_si)+scale(adi_com)+scale(adi_rb)+age+sex,df_vissearch_id_ASD))
        summary(lm(scale(hit_cat)~scale(ados_css)+scale(adi_si)+scale(adi_com)+scale(adi_rb)+age+sex,df_vissearch_id_ASD))

###--> INFORMATIVE PLOTS ####

# heat map of gaze positions across participants
d<-ggplot(vissearch,aes(x=gazepos.x,y=gazepos.y))
d+geom_hex(aes(fill=..density..))+ xlim(0,1)+ylim(1,0) + scale_fill_gradientn(colours=rev(rainbow(3)))+
  labs(x="Scaled x-coordinate", y="Scaled y-coordinate", title="visual search gaze position") +
  coord_fixed(ratio=1080/1920)+
  theme_bw()

# plot gazeposition with visualization of target hits
target.hits.var$target.hits <- as.factor(target.hits.var$target.hits)
levels(target.hits.var$target.hits) <- c("no hit", "hit")
d<-ggplot(vissearch,aes(x=gazepos.x,y=gazepos.y,color=target.hits.var$target.hits))

# #f?r Probanden getrennt
# d + xlim(0,1)+ylim(1,0) + geom_point(alpha=0.5) + annotate('rect',xmin = target.rect.xmin,xmax = target.rect.xmax, ymin = target.rect.ymin,ymax = target.rect.ymax, alpha = 0.8) +
#   labs(x="gazeposition on x-axis", y="gazeposition on y-axis", title="gazeposition and target hits per participant", color="target hits") + facet_wrap(~id)
# # ?ber Probanden hinweg -> keine vorher definierten "K?sten" f?r targets, sondern crosses
# d + xlim(0,1)+ylim(1,0)+ geom_point(alpha=0.5) + theme_bw() + annotate('point',x = target.x, y = target.y, size = 5, shape = 3, color = 'black') +
#   labs(x="Scaled x-coordinate", y="Scaled y-coordinate", title="gazeposition and target hits across participants", color="color") #+ coord_cartesian(xlim=c(0,1),ylim=c(0,1))
# # across participants mit "K?sten":
# d + xlim(0,1)+ylim(1,0)+ geom_point(alpha=0.5) + annotate('rect',xmin=target.rect.xmin, xmax = target.rect.xmax, ymin = target.rect.ymin,ymax = target.rect.ymax, alpha=0.5) +
#   labs(x="Scaled x-coordinate", y="Scaled y-coordinate", title="gazeposition and target hits across participants", color="target hits") #+ coord_cartesian(xlim=c(0,1),ylim=c(0,1))

# histogram: hits vs. no hits for participants
ggplot(target.hits.var, aes(x=target.hits)) + geom_histogram(binwidth=1) +
  labs(x = "hits vs. no hits", y = "number of hits", title ="target hits for each participant") +
  facet_wrap(~id)

# histogram: percentage hits per participant
target.hits.var_split <- split(target.hits.var, target.hits.var$id) # split into list
nobs<-sapply(target.hits.var_split,nrow)
number.of.hits.perid <- lapply(target.hits.var_split, function(x) {length(which(x$target.hits == "hit"))}) # number of hits per participants
percentage.hits <- lapply(1:length(number.of.hits.perid), function(x) {number.of.hits.perid[[x]]/nobs[x]}) # percentage of hits per participant
percentage.hits.df <- do.call(rbind, percentage.hits)
percentage.hits.df <- data.frame(levels(target.hits.var$id), percentage.hits.df)
names(percentage.hits.df) <- c("id", "percentage.hits")
ggplot(percentage.hits.df, aes(x=id, y=percentage.hits)) + geom_bar(stat="identity") +
  labs(x="Participant", y="Percentage hits", title="target hits per participant")+
  theme_bw()

#percent of hits  by group
with(percentage.hits.df,by(percentage.hits,grepl('_t',levels(id)),psych::describe))

    ###--> groups do not differ by percentage of target hits on the group level
    df_vissearch_agg<-data.frame(percentage.hits.df,nobs_percent)
    df_vissearch_agg$group<-as.factor(ifelse(grepl('_t',df_vissearch_agg$id)==T,'ASD','TD'))
    summary(lm(scale(percentage.hits)~group*scale(nobs_percent),data=df_vissearch_agg))


#--> 3. OUTCOME MEASURES

# number of target hits fÃ¼r einzelne target locations
by(target.hits.var$target.hits,target.fac,table)

# relative HITS by target.position
# erst mean pro proband pro condition errechnen und dann pro condition, damit nicht mehr werte fÃ¼r einen probanden eingehen
target.hits.var <- data.frame(target.hits.var, target.fac)
levels(target.hits.var$target.hits) <- c(0,1)
target.hits.var$target.hits <- as.numeric(target.hits.var$target.hits)
mean.per.id.cond <- by(target.hits.var$target.hits, interaction(target.hits.var$target.fac, target.hits.var$id), mean)
hits.per.id.cond <- data.frame(rep(levels(target.hits.var$id), each=length(levels(target.fac))),
                               rep(levels(target.hits.var$target.fac), length(levels(target.hits.var$id))),
                               as.numeric(mean.per.id.cond),
                               names(mean.per.id.cond))
names(hits.per.id.cond) <- c("id", "cond","mean", "cond.id")
condmean <- by(hits.per.id.cond$mean, hits.per.id.cond$cond, mean, na.rm=T)
condsd <- by(hits.per.id.cond$mean, hits.per.id.cond$cond, sd, na.rm=T)

# plotte mean and sd in "kÃ¤sten"
ggplot(vissearch,aes(x=gazepos.x, y=gazepos.y)) + xlim(0,1) + ylim(1,0) + theme_bw() +
  annotate('rect',xmin=target.rect.xmin, xmax = target.rect.xmax,
           ymin = target.rect.ymin,ymax = target.rect.ymax, alpha=0.4) +
  annotate(geom="text", x=target.x, y=target.y, size=3,
           label=c(paste("M =", round(condmean[[2]],2), "\n SD =", round(condsd[[2]],2)),
                   paste("M =", round(condmean[[4]],2), "\n SD =", round(condsd[[4]],2)),
                   paste("M =", round(condmean[[5]],2), "\n SD =", round(condsd[[5]],2)),
                   paste("M =", round(condmean[[7]],2), "\n SD =", round(condsd[[7]],2)),
                   paste("M =", round(condmean[[6]],2), "\n SD =", round(condsd[[6]],2)),
                   paste("M =", round(condmean[[8]],2), "\n SD =", round(condsd[[8]],2)),
                   paste("M =", round(condmean[[1]],2), "\n SD =", round(condsd[[1]],2)),
                   paste("M =", round(condmean[[3]],2), "\n SD =", round(condsd[[3]],2)))) +
  labs(x="Scaled x-axis of the screen", y="Scaled y-axis of the screen", title="Mean and standard deviation of target hits")

#FRAGE: target.x und target.y stimmen nicht mit faktorenlevels/pixel Ã¼berein (weichen jeweils um 100 ab)

##------ TASK: EMOTION EXPRESSION DATA ####

## 1. Create variables ####

# cut to relevant data
emoexp<-df[grep('emoexp_',df$eventlog),]
#emoexp<-df[which(df$eventlog %in% c('CF1','CF14','CF3','GF1','GF14','GF3','OF1','OF14','OF3','SF1','SF14','SF3')),]
# pd<-(emoexp$pupildil_L+emoexp$pupildil_R)/2

#drop unused levels
emoexp$eventlog<-droplevels(emoexp$eventlog)
table(emoexp$eventlog)

## calculate rpd (based on mean PD per participant)
df_rpd_id<-with(emoexp[emoexp$ts.event<50,],by(pd,id,mean,na.rm=T)) #mean of 200ms of every participant
df_rpd_id<-data.frame(as.numeric(df_rpd_id),names(df_rpd_id))
names(df_rpd_id)<-c('tpd','id')
emoexp<-merge(emoexp,df_rpd_id,by='id')
emoexp$rpd<-emoexp$pd/emoexp$tpd

## define categories
emo.cat<-data.frame(emoexp$id, rep(NA,nrow(emoexp)))
names(emo.cat) <- c("id", "emo.cat")
emo.cat$emo.cat[c(grep('CF1',emoexp$eventlog),grep('OF1',emoexp$eventlog),grep('SF1',emoexp$eventlog),grep('GF1',emoexp$eventlog))]<-1
emo.cat$emo.cat[c(grep('CF3',emoexp$eventlog),grep('OF3',emoexp$eventlog),grep('SF3',emoexp$eventlog),grep('GF3',emoexp$eventlog))]<-2
emo.cat$emo.cat[c(grep('CF14',emoexp$eventlog),grep('OF14',emoexp$eventlog),grep('SF14',emoexp$eventlog),grep('GF14',emoexp$eventlog))]<-3
emo.cat$emo.cat<-as.factor(emo.cat$emo.cat)
levels(emo.cat$emo.cat)<-c('happy','afraid','neutral')

emo_cat<-emo.cat$emo.cat
emoexp<-data.frame(emoexp,emo_cat)

## PCA
df_pca<-emoexp[,c('id','eventlog','ts.event','pd','emo_cat')]
df_pca<-df_pca[df_pca$ts.event<=1800,] #cut after 6 seconds (5*300Hz = 1800)
df_pca$ts.event_rounded<-round(df_pca$ts.event/10) #reduce data
PCAdata<-reshape2::dcast(df_pca,formula = id+emo_cat~ts.event_rounded, fun.aggregate = median, na.rm=T, value.var ='pd')

  #impute missing data
  PCAdata<-missMDA::imputePCA(PCAdata[,3:length(names(PCAdata))],ncp=2) #impute by expected PCA components (missMDA package)
  PCAmodel<-psych::principal(data.frame(PCAdata['completeObs']),nfactors = 2,rotate='varimax')

  #PCAmodel

  plot(PCAmodel$loadings[,1])
  plot(PCAmodel$loadings[,2])

  #merge to data frame (emoexp)
  df_loadings<-data.frame(1:181,PCAmodel$loadings[,1],PCAmodel$loadings[,2])
  names(df_loadings)<-c('ts.event_rounded','RC1_loading','RC2_loading')

  emoexp$ts.event_rounded<-round(emoexp$ts.event/10)
  emoexp<-merge(emoexp,df_loadings,by='ts.event_rounded',all.x=T) #lost data? --> at the end, after 6 seconds

  #define factor loadeding weighted rpd variable
  emoexp$rpd_RC1<-emoexp$rpd*emoexp$RC1_loading
  emoexp$rpd_RC2<-emoexp$rpd*emoexp$RC2_loading

  emoexp$rpd_RC1_z<-scale(emoexp$rpd*emoexp$RC1_loading)
  emoexp$rpd_RC2_z<-scale(emoexp$rpd*emoexp$RC2_loading)

##gaze center deviation
center_dev<-with(emoexp,sqrt((gazepos.x-0.5)^2+(gazepos.y-0.5)^2))
hist(center_dev)

###aggregate to trials

int_trial<-with(emoexp,interaction(id,eventlog)) #eventlog needed as eventlogs occure two times per id

factors_trial<-with(emoexp,aggregate(data.frame(group,id,eventlog,index.event),
                                        by=list(int_trial),head,n=1))

metrics_mean_trial<-with(emoexp,aggregate(data.frame(pd,ts,ts.event,gazepos.x,gazepos.y,center_dev,rpd_RC1,rpd_RC2),
                                             by=list(int_trial),mean,na.rm=T))

#data quality variable
recorded_samples_per_trial<-table(as.numeric(int_trial)) #sample per event
data_quality<-data.frame(recorded_samples_per_trial,factors_trial$eventlog)
names(data_quality)<-c('row','freq','eventlog')
max_samples<-with(emoexp,by(ts.event,eventlog,max,na.rm=T)) #max samples per video
data_quality2<-data.frame(as.numeric(max_samples),attributes(max_samples)$dimnames)
data_quality<-merge(data_quality,data_quality2,by='eventlog')
data_quality<-data_quality[order(data_quality$eventlog),] #sort in the same order as factors_trial
data_quality<-data_quality$freq/data_quality$as.numeric.max_samples.
#hist(data_quality)

#emotion category
emo_cat<-with(factors_trial,ifelse(grepl('14',eventlog),'neutral',
                          ifelse(grepl('3',eventlog),'afraid','happy')))

#stimulus category --> displayed individual

stim_cat<-with(factors_trial,ifelse(grepl('CF',eventlog),'man',
                                    ifelse(grepl('GF',eventlog),'girl1',
                                           ifelse(grepl('SF',eventlog),'boy','girl2'))))

#pupillary response --> based on visual inspection, see above
pd_peak<-with(emoexp,as.numeric(by(pd[ts.event>600 & ts.event<=750],
                                      droplevels(int_trial)[ts.event>600 & ts.event<=750],
                                      mean,na.rm=T)))


pd_baseline<-with(emoexp,as.numeric(by(pd[ts.event<=100],
                                          droplevels(int_trial)[ts.event<=100],
                                          mean,na.rm=T)))

pd_response<-pd_peak-pd_baseline

#relative pupillary response --> rpd is corrected by mean per participant

rpd_peak<-with(emoexp,as.numeric(by(rpd[ts.event>600 & ts.event<=750],
                                   droplevels(int_trial)[ts.event>600 & ts.event<=750],
                                   mean,na.rm=T)))


rpd_baseline<-with(emoexp,as.numeric(by(rpd[ts.event<=100],
                                       droplevels(int_trial)[ts.event<=100],
                                       mean,na.rm=T)))

rpd_response<-rpd_peak-rpd_baseline

##--> create data-frame
df_emoexp_trial<-data.frame(factors_trial,emo_cat,stim_cat,data_quality,pd_response,pd_baseline,rpd_response,metrics_mean_trial)

#create trial number
df_emoexp_trial<-df_emoexp_trial[order(df_emoexp_trial$id,df_emoexp_trial$index.event),] #sort order for calculating trial number
#problem: intermediate trials might be missing, thus calculate change value and use cumsum to create trial number
trial_change<-as.numeric(with(df_emoexp_trial,round(diff(index.event)/2))) #differences in index.event indicate trial change (divided by two as 2 additional events occur after each eventlog before next trial)
trial_change<-c(1,trial_change) #add 1 at beginning as diff is shorter than total length
trial_change[trial_change<0]<-1 # <0 == change in participants (id) in sequence
trial_change<-ifelse(trial_change>1,1,trial_change) #control for changes in blocks
trial_number<-as.numeric(unlist(with(df_emoexp_trial,by(trial_change,id,cumsum))))

df_emoexp_trial<-data.frame(df_emoexp_trial,trial_number)

## 2. merge with demographics ####

df_emoexp_trial<-df_emoexp_trial[,-grep('group',names(df_emoexp_trial))]
df_emoexp_trial$id_core<-substr(df_emoexp_trial$id,1,3)

table(df_emoexp_trial$id_core)
table(df_describe$id_core) #TODO: new data export, not all participants are included

df_emoexp_trial<-merge(df_emoexp_trial,df_describe,by='id_core')

##--> visualization ####

##pupillary response:

##a.) pupillary response by group and emotion
  ggplot(emoexp[emoexp$ts.event<1800,],aes(x=ts.event/300,y=rpd,group=interaction(emo_cat,group),color=emo_cat,linetype=group)) + geom_smooth() +
    labs (x = "time (s)", y = "realtive pupil dilation", title = "pupil dilation over emotion conditions")+theme_bw()

  #by group
  ggplot(emoexp[emoexp$ts.event<2000,],aes(x=ts.event/300,y=rpd,group=group,color=group)) + geom_smooth() +
    labs (x = "time (s)", y = "realtive pupil dilation", title = "pupil dilation over emotion conditions")+theme_bw()

  ggplot(emoexp[emoexp$ts.event<2000,],aes(x=ts.event/300,y=rpd)) + geom_smooth() +
    labs (x = "time (s)", y = "realtive pupil dilation", title = "pupil dilation over emotion conditions")+theme_bw()

#b.) heat map of gaze positions across participants
  d<-ggplot(emoexp,aes(x=gazepos.x,y=gazepos.y))
  d+geom_hex(aes(fill=..density..))+ xlim(0,1)+ylim(1,0) + scale_fill_gradientn(colours=rev(rainbow(3)))+
    labs(x="Scaled x-coordinate", y="Scaled y-coordinate", title="visual search gaze position") +
    coord_fixed(ratio=1080/1920)+
    theme_bw()

#c.) figure PCA loadings
  loads<-c(PCAmodel$loadings[,1],PCAmodel$loadings[,2])
  ncol_matrix<-ncol(data.frame(PCAdata['completeObs']))
  labels<-c(rep('RC1',ncol_matrix),rep('RC2',ncol_matrix))
  id<-c(rep(1:ncol_matrix,2))

  g1<-ggplot(data=data.frame(loads,labels,id),aes(x=id*10,y=loads,group=labels,color=labels,fill=labels))+
    geom_smooth(se=F) +theme_bw() +ggtitle('PCA: Factor Loadings') + theme(legend.position = "none",axis.title.x=element_blank()) +
    annotate('text', x=3000, y=0.8, color='#F8766D' ,label='RC1')+
    annotate('text', x=3000, y=0.5, color='#00BFC4' ,label='RC2')

  g2<-ggplot(emoexp,aes(x=ts.event,y=rpd_RC1))+geom_smooth(se=F,method='lm', formula = y ~ x + poly(x,3),color='#F8766D')+theme(axis.title.y =   element_blank(),axis.title.x = element_blank(),legend.position = "none")+ggtitle('Early Pupillary Response (PR1)')

  g3<-ggplot(emoexp,aes(x=ts.event,y=rpd_RC2))+geom_smooth(se=F,method='lm', formula = y ~ x + poly(x,3),color='#00BFC4')+theme(axis.title.y = element_blank(),axis.title.x =  element_blank(),legend.position = "none")+ggtitle('Late Pupillary Response (PR2)')

  grid.arrange(g1,g2,g3,ncol=3, bottom ='scene duration (ms)')

###--> preliminary analysis####

###-- a.) random effect structure ####

re_model1<-lmer(pd_response~1+(1|id),df_emoexp_trial,REML=F)
re_model2<-lmer(pd_response~1+(1|trial_number)+(1|id),df_emoexp_trial,REML=F)
re_model3<-lmer(pd_response~1+(1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial,REML=F)
anova(re_model1,re_model2,re_model3)

###-- b.) effect of task --> control for missing data (data quality) and center deviation ####

lmm_model<-lmer(data_quality~group+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

contrast(emmeans(lmm_model,~group),'pairwise') #higher data quality in TD

lmm_model<-lmer(center_dev~group+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

contrast(emmeans(lmm_model,~group),'pairwise') #higher center deviation in ASD

lmm_model<-lmer(pd_baseline~group+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model) #-->no group difference in PD baseline

###-- c.) effect of group ####

lmm_model<-lmer(rpd_response~emo_cat*group+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model)
contrast(emmeans(lmm_model,~group),'pairwise') #higher pupillary response in TD


lmm_model<-lmer(pd_response~emo_cat*group+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model)
fixef(lmm_model)
contrast(emmeans(lmm_model,~group),'pairwise') #higher pupillary response in TD

contrast(emmeans(lmm_model,~group|emo_cat),'pairwise',adjust='none') #between: response difference for neutral stimuli
plot(emmeans(lmm_model,~group|emo_cat))

contrast(emmeans(lmm_model,~emo_cat|group),'pairwise',adjust='none') #within: ASD response difference between afraid - neutral
##--> does this indicate that ASD is reactive to the emotional content but not the faces?

lmm_model<-lmer(pd_response~emo_cat*group+
                  data_quality+center_dev+
                  age+iq+sex+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model) #-->inclcusion of covariates does not alter effects
r.squaredGLMM(lmm_model)
plot(emmeans(lmm_model,~group+emo_cat))

emmeans(lmm_model,~group+emo_cat) #groups still differ in neutral

###

lmm_model<-lmer(pd_response~srs_total+rbsr_total+
                data_quality+center_dev+
                  age+iq+sex+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

summary(lmm_model)

lmm_model<-lmer(scale(pd_response)~scale(ados_css)*emo_cat+
                  data_quality+center_dev+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model)
emmeans(lmm_model,~emo_cat)
emtrends(lmm_model,~emo_cat,var='ados_css') #higherr ados associated with larger response on neutral stimuli

## response for pupillary components

lmm_model<-lmer(scale(rpd_RC1)~group*emo_cat+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model)
contrast(emmeans(lmm_model,~emo_cat),'pairwise',adjust='none')
#-->larger late weighted pupillary response for afraid compared to neutral stimuli

lmm_model<-lmer(scale(rpd_RC2)~group*emo_cat+
                  (1|stim_cat)+(1|trial_number)+(1|id),df_emoexp_trial)

anova(lmm_model)
emmeans(lmm_model,~emo_cat)
contrast(emmeans(lmm_model,~emo_cat),'pairwise',adjust='none')
#--> larger early weighted pupillary response for neutral compared to happy and afraid

# #--> OLD - emo exp ####
#
# # onscreen time ?ber kategorien hinweg
# samples_per_event<-with(emoexp,by(ts.event,droplevels(eventlog),max,na.rm=T))
#
# emoexp_split<-split(emoexp,emoexp$id)
# nobs<-sapply(emoexp_split,nrow) #number of observations
# #possible_nobs<- sum(samples_per_event[grep('emoexp_', names(samples_per_event))])
# possible_nobs<- sum(samples_per_event)
# nobs_percent<-nobs/possible_nobs*100
# round(nobs_percent,2)
#
# psych::describe((round(nobs_percent,2)))
#
# #ggplot(emoexp,aes(x=id,y=pd)) + geom_boxplot() #boxplot der rpd ?ber alle conditions hinweg (f?r ids einzeln)
#
# #### onscreen time f?r kategorien getrennt
# emoexp2 <- data.frame(emoexp, emo.cat)
# happy <- emoexp2[emoexp2$emo.cat == "happy",]
# happy_split <- split(happy, happy$id)             #list
# afraid <- emoexp2[emoexp2$emo.cat == "afraid",]
# afraid_split <- split(afraid, afraid$id)             #list
# neutral <- emoexp2[emoexp2$emo.cat == "neutral",]
# neutral_split <- split(neutral, neutral$id)             #list
#
# # percentage onscreen time: happy
# nobs<-sapply(happy_split,nrow) #number of observations
# possible_nobs<-sum(samples_per_event[c(grep("^emoexp_CF1$", names(samples_per_event)),
#                                        grep("^emoexp_OF1$", names(samples_per_event)),
#                                        grep("^emoexp_SF1$", names(samples_per_event)),
#                                        grep("^emoexp_GF1$", names(samples_per_event)))])
# nobs_percent_happy<-nobs/possible_nobs*100
# round(nobs_percent_happy,2)
# round(mean(nobs_percent_happy),2)
# psych::describe(nobs_percent_happy)
#
# # percentage onscreen time: afraid
# nobs<-sapply(afraid_split,nrow) #number of observations
# possible_nobs<-sum(samples_per_event[c(grep("^emoexp_CF3$", names(samples_per_event)),
#                                        grep("^emoexp_OF3$", names(samples_per_event)),
#                                        grep("^emoexp_SF3$", names(samples_per_event)),
#                                        grep("^emoexp_GF3$", names(samples_per_event)))])
# nobs_percent_afraid<-nobs/possible_nobs*100
# round(nobs_percent_afraid,2)
# round(mean(nobs_percent_afraid),2)
# psych::describe(nobs_percent_afraid)
#
#
# # percentage onscreen time: neutral
# nobs<-sapply(afraid_split,nrow) #number of observations
# possible_nobs<-sum(samples_per_event[c(grep("^emoexp_CF14$", names(samples_per_event)),
#                                        grep("^emoexp_OF14$", names(samples_per_event)),
#                                        grep("^emoexp_SF14$", names(samples_per_event)),
#                                        grep("^emoexp_GF14$", names(samples_per_event)))])
# nobs_percent_neutral<-nobs/possible_nobs*100
# round(nobs_percent_neutral,2)
# round(mean(nobs_percent_neutral),2)
# psych::describe(nobs_percent_neutral)
#
#
# #test if there are significant differences:
# #create data frame for inference statistics
# nobs_percent_cat <- data.frame(rep(levels(emoexp$id), 3),
#                                c(nobs_percent_afraid, nobs_percent_happy, nobs_percent_neutral),
#                                rep(c("afraid", "happy", "neutral"), each=length(levels(emoexp$id))))
# names(nobs_percent_cat) <- c("id", "nobs_percent", "category")
# require(lmer)
# require(lmerTest)
# #require(contrast)
# #require(lsmeans)
# mixed.model<-lmer(nobs_percent_cat$nobs_percent~nobs_percent_cat$category+(1|nobs_percent_cat$id),REML=T)
# #summary(mixed.model)
# anova(mixed.model)
# contrast(pmmeans(mixed.model,~category), method="pairwise")
#
# by(nobs_percent_cat$nobs_percent, nobs_percent_cat$category, psych::describe)
# # higher for neutral than happy and afraid
#
#
# #--> 2. INFORMATIVE PLOTS
#
# ggplot(data = emoexp2, aes(x = gazepos.x, y = gazepos.y))+ xlim(0,1)+ylim(1,0)+
#   geom_hex(aes(fill=..density..)) + scale_fill_gradientn(colours=rev(rainbow(3))) + scale_colour_gradientn(colours=rev(rainbow(3))) +
#   labs(x = 'Scaled x-coordinate', y='Scaled y-coordinate', title="valid fixation data") +
#   theme_bw()
#
#
#
# # # gazeposition nach Probanden plotten (EIN Plot)
# # ggplot(data = emoexp2, aes(x = gazepos.x, y = gazepos.y))+ xlim(0,1)+ylim(1,0)+
# #   geom_hex() + scale_fill_gradientn(colours=rev(rainbow(3))) +
# #   labs(x = 'scaled x-axis', y='scaled y-axis', title='valid fixation data', color= 'participants')
#
# # # heat map f?r probanden getrennt
# # g<-ggplot(emoexp2,aes(x=gazepos.x,y=gazepos.y))
# # g+geom_hex(bins=30,aes(fill=..density..))+xlim(0,1)+ylim(1,0)+
# #   facet_wrap(~id) + theme_bw()+ scale_fill_gradientn(colours=rev(rainbow(3))) +
# #   labs(x = 'screen x-axis', y='screen y-axis', title='heat map of valid fixation data')
#
# # Boxplot f?r rpd ?ber Kategorien (f?r Probanden getrennt)
# ggplot(data.frame(emoexp2,emo.cat, rpd),aes(x=emo.cat,y=rpd)) + geom_boxplot() + ylim(0.7, 1.5) +
#   facet_wrap(~id)  +
#   labs (x = "emotion condition", y = "pupil dilation", title = "pupil dilation over emotion conditions")
#
#
# # Plotte MW +-sd f?r Probanden getrennt in einen Plot
# ggplot(emoexp2,aes(x=emo.cat,y=rpd,group=id,color=id))+
#   stat_summary(fun.y = mean,
#                fun.ymin = function(x) mean(x) - sd(x),
#                fun.ymax = function(x) mean(x) + sd(x),
#                geom = "pointrange") +
#   stat_summary(fun.y = mean,
#                geom = "line") +
#   labs(x="emotion condition", y="pupil dilation", title="pupil dilation over emotion condition", color ="participant")
#
# # rpd over time
# shortest_video_time <- min(samples_per_event[grep("emoexp", names(samples_per_event))]) * 0.0033
# shortest_video <- which(samples_per_event[grep("emoexp", names(samples_per_event))] == shortest_video_time)
#
# #macht das Ã¼ber probanden hinweg Ã¼berhaupt sinn? da die mit mehr observations dann mehr eingehen?
# ggplot(emoexp2,aes(x=ts.event*0.0033,y=rpd,color=emo.cat))+ xlim(0,shortest_video_time) +
#   geom_smooth() +
#   theme_bw() +
#   labs(x='Duration (in seconds)',y='Relative pupil dilation', colour="Condition")
#
# #rpd over time for each participant
# ggplot(emoexp2,aes(x=ts.event*0.0033,y=rpd,color=emo.cat))+ xlim(0,shortest_video_time) +
#   geom_smooth() + facet_wrap(~id) + theme_bw() +
#   labs(x='duration (in seconds)',y='pupil dilation (in mm)')
#
# #--> 3. TASK OUTCOME MEASURE
#
# #f?r Probanden getrennt: decriptives der rpd nach condition getrennt
# with(emoexp,by(pd,emo.cat$emo.cat, psych::describe))
#
# #FRAGE: warum geht erste zeile nicht mehr
# rpd.byn<-as.numeric(by(emoexp2$rpd,droplevels(interaction(emoexp2$emo.cat,emoexp2$id,emoexp2$index.event)),mean, na.rm=T))
# rpd.byn<-as.numeric(by(rpd,droplevels(interaction(emoexp2$emo.cat,emoexp2$id,emoexp2$index.event)),function(x) {mean(x,na.rm=T)}))
# emocat.byn<-as.numeric(by(emoexp2$emo.cat,droplevels(interaction(emoexp2$emo.cat,emoexp2$id,emoexp2$index.event)),head,n=1))
# emocat.byn <- as.factor(emocat.byn)
# levels(emocat.byn) <- c("happy", "afraid", "neutral")
# id.byn<-as.numeric(by(emoexp2$id,droplevels(interaction(emoexp2$emo.cat,emoexp2$id,emoexp2$index.event)),head,n=1))
# levels(emocat.byn) <- c("happy", "afraid", "neutral")
#
# ggplot(data.frame(rpd.byn,emocat.byn),aes(x=emocat.byn,y=rpd.byn,group=emocat.byn, color=emocat.byn)) +  ylim(0.75, 1.25) +
#   geom_boxplot() + theme_bw() + labs(y="Relative pupil dilation", x="Condition", title="Relative pupil dilation per condition", color="Condition")
#
# id.byn <- as.factor(id.byn)
# by(rpd.byn, emocat.byn, psych::describe) # descriptives Ã¼ber Ã¼ber participants hinweg nach condition getrennt
#
# # #descriptives der pupil dilation ?ber alle Probanden hinweg -> wahrscheinlich wenig Sinn wg. untersch. Lichtverh?ltnisse
# # emo.cat2<-rep(NA,nrow(emoexp2))
# # emo.cat2[c(grep('CF1',emoexp$eventlog),grep('OF1',emoexp$eventlog),grep('SF1',emoexp$eventlog),grep('GF1',emoexp$eventlog))]<-1
# # emo.cat2[c(grep('CF3',emoexp$eventlog),grep('OF3',emoexp$eventlog),grep('SF3',emoexp$eventlog),grep('GF3',emoexp$eventlog))]<-2
# # emo.cat2[c(grep('CF14',emoexp$eventlog),grep('OF14',emoexp$eventlog),grep('SF14',emoexp$eventlog),grep('GF14',emoexp$eventlog))]<-3
# # emo.cat2<-as.factor(emo.cat2)
# # levels(emo.cat2)<-c('happy','afraid','neutral')
# # by(rpd,emo.cat2,psych::describe) # ?ber alle Probanden hinweg
# # ggplot(data.frame(rpd,emo.cat2), aes(x=emo.cat2, y=rpd)) + geom_boxplot()
#
# #--> 4. INFERENCE STATISTICS
#
# #unterschiede rpd ~ category
# require(lme4)
# require(lmerTest)
# id.byn <- as.factor(id.byn)
# mixed.model<-lmer(scale(rpd.byn)~emocat.byn+(1|id.byn),REML=T)
# anova(mixed.model)
# contrast(pmmeans(mixed.model,~emocat.byn), method="pairwise")
#
# #regression
# #1=happy, 2=afraid, 3=neutral
# colvar.inf <- data.frame(id.byn, rpd.byn,emocat.byn)
# colvar.inf$id.byn <- as.factor(colvar.inf$id.byn)
#
# diff.happyneutral <- as.numeric(by(colvar.inf, colvar.inf$id, function(x) {
#   mean(x$rpd.byn[x$emocat.byn==1], na.rm=T) - mean(x$rpd.byn[x$emocat.byn == 3], na.rm=T)}))
# diff.afraidneutral <- as.numeric(by(colvar.inf, colvar.inf$id, function(x) {
#   mean(x$rpd.byn[x$emocat.byn==2], na.rm=T) - mean(x$rpd.byn[x$emocat.byn == 3], na.rm=T)}))
#
# #difference happy vs. neutral ~ controlvar
# controlvar1 <- controlvar
# fit <- lm(diff.happyneutral ~ controlvar1$age_months +
#             controlvar1$ADOS_Vergleichswert +
#             controlvar1$cognitive_fun, data=data.frame(controlvar1, diff.happyneutral))
# summary(fit)
# #asbolute difference values
# fit <- lm(abs(diff.happyneutral) ~ controlvar1$age_months +
#             controlvar1$ADOS_Vergleichswert +
#             controlvar1$cognitive_fun, data=data.frame(controlvar1, diff.happyneutral))
# summary(fit)
#
# #difference afraid vs. neutral
# fit <- lm(diff.afraidneutral ~ controlvar1$age_months +
#             controlvar1$ADOS_Vergleichswert +
#             controlvar1$cognitive_fun, data=data.frame(controlvar1, diff.afraidneutral))
# summary(fit)
# #absolute difference values
# fit <- lm(abs(diff.afraidneutral) ~ controlvar1$age_months +
#             controlvar1$ADOS_Vergleichswert +
#             controlvar1$cognitive_fun, data=data.frame(controlvar1, diff.afraidneutral))
# summary(fit)
#
#
#
#

#--------- TASK: colscreen+fixcross
colscreen<-df[c(grep('fixationcross',df$eventlog),grep('colscreen',df$eventlog)),]
colscreen<-colscreen[-grep('dur:5',colscreen$eventlog),]
colscreen<-colscreen[-grep('dur:7.5',colscreen$eventlog),]

#colscreen<-df[c(grep('colscreen_col:1_dur:2',df$eventlog)),]

table(droplevels(colscreen$eventlog))

#reconstruct event --> create trial number
colscreen<-colscreen[order(colscreen$id,colscreen$timestamp),]
  #problem: intermediate trials might be missing, thus calculate change value and use cumsum to create trial number
  trial_change<-as.numeric(with(colscreen,diff(index.event))) #differences in index.event indicate trial change
  trial_change<-c(0,trial_change)
  trial_change<-ifelse(trial_change>1,1,0) #control for changes in blocks
  trial_number<-as.numeric(unlist(with(colscreen,by(trial_change,id,cumsum))))
  trial_number<-trial_number+1

  #TODO: trial specific timestamp (ts)
  #ts.trial<-lapply(trial_number,function(x){as.numeric(do.call(c,by(x,x,seq_along, simplify=F)))})
  test<-by(trial_number,interaction(colscreen$id,trial_number),seq_along)
  #already exists: --> see ts.event
  hist(colscreen$ts.event[colscreen$ts.event<200])

## calculate rpd (based on mean PD per participant)
df_rpd_id<-with(colscreen[colscreen$ts.event<50,],by(pd,id,mean,na.rm=T)) #mean of 200ms of every participant
df_rpd_id<-data.frame(as.numeric(df_rpd_id),names(df_rpd_id))
names(df_rpd_id)<-c('tpd','id')
tmp<-merge(colscreen[,c('id','pd')],df_rpd_id,by='id') #memory overflow

colscreen$rpd<-tmp$pd/tmp$tpd
colscreen$tpd<-tmp$tpd

#TODO: plot on machine with more memory
ggplot(colscreen,aes(x=ts.event,y=rpd)) + geom_smooth() +
  labs (x = "time (s)", y = "realtive pupil dilation", title = "pupil dilation over evvntlog")+theme_bw()

## ------------------ TASK: NATORIENT ######
natorient<-df[grep('natorient',df$eventlog),]
natorient<-natorient[natorient$eventlog %in% c('natorient_childrenparty','natorient_childrenplay','natorient_childrenplaying','natorient_childrenballbath'),]
natorient$eventlog<-droplevels(natorient$eventlog)
table(natorient$eventlog)

#percent data of natorient compared to total data:
nrow(natorient)/nrow(df)

#individual data per video
with(natorient,by(id,eventlog,function(x){length(unique(x))}))

