# -------- FAST TASK ANALYSIS ------- #
# from MIND data - data is provided upon reasonable request: nico.bast@kgu.de
  
  #info four-choice reaction time task with slow-baseline and fast-incentive condition

# 0. SETUP - required packages ####
 #data processing
require(plyr) #rbind.fill
require(zoo) #na.approx

  #visualization
require(ggplot2) #visualization
require(grid) #arrange ggplot objects
require(gridExtra) #arrange ggplot objects
 
 #analysis
require(nlme) #linear mixed models
require(lmerTest) #LMM p-values
require(emmeans) #post hoc testing in LMM

# A. READ DATA ####
## - GET data paths ####
setwd('C:/Users/Nico/PowerFolders/Supervision_SaraMIND/all data/sboxhoorn/')
data.files<-list.files(path=getwd(), full.names=TRUE,recursive=T) #lsit all data within folders/subfolders
data_ft<-data.files[grep('FastTask/log',data.files)] #retreieve fast task paths

## - GET participant IDs ####
## --> read the first 5 lines from all the WM textfiles to extract patient ids braingaze 
function_read_id<-function(x){
  k<-readLines(con=x, '5')
  pt.id<-k[grep('PatientId:',k)]
  pt.id<-substring(pt.id, 12,nchar(pt.id))
  pt.id<-sub('\"',"", pt.id) #saras code erased the last character of the patient id - corrected
  se.id<-k[grep('SessionId:',k)]
  se.id<-substring(se.id, 12,nchar(se.id))
  se.id<-sub('\"',"", se.id)
  return(pt.id)}
names<-as.character(lapply(data_ft,function_read_id))

length(names) #data of n=105

## - READ data to list ####    
function_read_data<-function(x){
  #k<-read.table(x,skip=23,fill=T,sep=';',header=T,dec=',')
  k<-try(read.table(x,skip=23,fill=T,sep=';',header=T,dec=','))
  if(!inherits(k,'try-error')) k #workaround to embedded nul issues
  #--> altered this function as in R > 3.4 read.table will give and error and stop for embedded nuls in file, maybe a file encoding issue? 
  return(k)}

start.time <- Sys.time()
list_data<-lapply(data_ft,function_read_data) #read all relevantdata files 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # 1.2 min for 105 datasets
names(list_data)<-names #save patient IDs to loaded data

# B. CREATE ESSENTIAL VARIABLES --> FORMAT LIST OF DATA TO MAIN DATA FRAME ####
## - EXCLUDE participants without patient IDs (mainly test data) ####
list_data<-list_data[-which(names(list_data)=='character(0)'|names(list_data)=='00000000-0000-0000-0000-000000000000')]
#29.08.19:: --> n=99 patient data

## - EXCLUDE participants with more than one session ####
duplicated_ids<-names(list_data)[duplicated(names(list_data))]
duplicated_ids
#case 1: "4a7ed1c6-237d-4e56-9bb0-3872f18408d0" #--> restarted session, exclude the one with few entries
size_list_data<-as.numeric(sapply(list_data,nrow))
excluded_session1<-which(names(list_data)=='4a7ed1c6-237d-4e56-9bb0-3872f18408d0' & size_list_data<3000)
#--> removes the session with few data

#case 2: "3272b2e4-db3a-4b3d-a71b-a60b929ae585" #--> probably did the test twice, exclude the one that was done later
function_read_date<-function(x){
  date_info<-readLines(con=x, '1')
  date_info<-date_info[grep('StartTime:',date_info)]
  date<-substring(date_info, 12,nchar(date_info)-7) # '-7' erases the time information
  #time<-substring(date_info, nchar(date_info)-5,nchar(date_info)-1)
  return(date)}
start_date<-as.character(sapply(data_ft,function_read_date)) #when was the data assessed?
excluded_session2<-which(names(list_data)=='3272b2e4-db3a-4b3d-a71b-a60b929ae585' & start_date=='2016.09.01')

list_data<-list_data[-c(excluded_session1,excluded_session2)]

## - CREATE an id and timestmap variable for flat data format ####
id.length<-sapply(list_data,nrow)
ts.var<-as.numeric(unlist(sapply(id.length,seq_len)))
id.label<-rep(names(id.length),id.length)

#individual trials are typegame x currentTrial
#trial.length<-unlist(sapply(list_data,function(x){table(x$currentTrial)}))
trial.length<-unlist(sapply(list_data,function(x){table(interaction(x$typeGame,x$currentTrial))}))

## - FORMAT list to data.frame (df_data) and add id + timestamp variable ####
#df_data<-do.call(rbind,list_data) #takes longer
df_data<-rbind.fill(list_data) #faster -- see: https://www.r-bloggers.com/the-rbinding-race-for-vs-do-call-vs-rbind-fill/
df_data<-cbind(id.label,ts.var,df_data)

## - CREATE timestamp within trial (ts.trial) variable and add to df_data ####
trial.id<-with(df_data,interaction(id.label,typeGame,currentTrial))
ts.trial<-unlist(sapply(trial.length,seq_len))

#remove first row of data of participants that contain no data
# otherwise cbind with ts.trial is not possible
df_data<-df_data[-which(is.na(trial.id)),]
df_data<-cbind(df_data,ts.trial)  

## - QC: RAW DATA QUALITY CONTROL ####
## - 1) turn -1 values into NA ###
df_data[df_data==-1]<-NA

## - 2) remove trial with length >360 (~12 seconds) ###
hist(trial.length[trial.length<600],50) #typical trial length for baseline (~100) and fasttask (~240)
cutoff<-360 #trials of 12 seconds length
trial.valid<-ifelse(trial.length>cutoff,F,T)
trial.valid<-rep(trial.valid,times=trial.length)
df_data<-df_data[trial.valid,]
  
## - 3) remove unecessary variables form workspace ###
rm(id.label,ts.var,ts.trial,trial.id,id.length,trial.length)

## --> FIRST DATA VISUALIZATION - UNPROCESSED DATA ####
with(df_data,hist(leftEyePupilSize,100))
with(df_data,hist(rightEyePupilSize,100))
### --> probabaly set size < 3mm to NA

# does the data make sense?
table(df_data$strState)
with(df_data[df_data$typeGame=='fastask',],hist(ts.trial,100))
with(df_data[df_data$typeGame=='baseline',],hist(ts.trial,100)) #--> why does this follow the same distribution?

with(df_data[df_data$strState=='foreperiod' & df_data$ts.trial<100,],hist(ts.trial,100))
with(df_data[df_data$strState=='isi' & df_data$ts.trial<100,],hist(ts.trial,100))
with(df_data[df_data$strState=='stimuli' & df_data$ts.trial<100,],hist(ts.trial,100))
#--> within trial ordering seems to be foreperiod --> stimuli --> isi

#when does the stimulus occur?
with(df_data[df_data$strState=='stimuli' & df_data$ts.trial<100 & df_data$typeGame=='fastask',],hist(ts.trial,100))
with(df_data[df_data$strState=='stimuli' & df_data$ts.trial<100 & df_data$typeGame=='baseline',],hist(ts.trial,100))
with(df_data[df_data$strState=='stimuli' & df_data$ts.trial<300 & df_data$typeGame=='baseline',],hist(ts.trial,100))
#--> in the fastask condition around the 30 sample mark (== 1s) as expected
#--> cant see stimulus signal in baseline ???

#difference between fastask and baseline - Task-evoked PD
g<-ggplot(df_data[df_data$ts.trial<60,],aes(x=ts.trial,y=rightEyePupilSize,group=typeGame,color=typeGame))
g+geom_smooth(se=F)+stat_summary(fun.data = mean_se, geom = "errorbar")
#--> substanitally higher PD during fastask and TEPR signal in tastask

    #baseline PD progression (baseline has to wait 8 seconds before reaction)
    g<-ggplot(df_data[df_data$typeGame=='baseline' & df_data$ts.trial<270,],aes(x=ts.trial,y=rightEyePupilSize))
    g+geom_smooth(se=F)+stat_summary(fun.data = mean_se, geom = "errorbar")
    #--> why is there a dip at 110sample == 3.5 seconds?
    ##--> participants look away from the center
    
    #fastask PD progression (baseline has to wait 8 seconds before reaction)
    g<-ggplot(df_data[df_data$typeGame=='fastask' & df_data$ts.trial<60,],aes(x=ts.trial,y=rightEyePupilSize))
    g+geom_smooth(se=F)+stat_summary(fun.data = mean_se, geom = "errorbar")

#Mean pupil dilation over trials
g<-ggplot(df_data[df_data$ts.trial<270,],aes(x=currentTrial,y=rightEyePupilSize,group=typeGame,color=typeGame))
g+geom_smooth(se=F)+stat_summary(fun.data = mean_se, geom = "errorbar")


# C. MERGE WITH DEMOGRAPHICS DATA - group, sex, age, etc. ####
#read demographics+group-identifier (dm) and merge to data 
#ET data that does not correspond to demographics data is excluded by merge
df_demo<-read.csv('C:/Users/Nico/PowerFolders/Supervision_SaraMIND/ParticipantsMiND.CSV',sep=';',header = T)
df_demo<-df_demo[,c(1:7,9,10)] #clear data - questionnaire data will be added later as it will use MICEd data
df_data<-merge(df_data,df_demo,by.x='id.label',by.y='PatientID..Braingaze.')

length(table(droplevels(df_data$id.label))) #after demographics data merging: n=97

# D. MERGE WITH QUESTIONNAIRE DATA (MICE imputed) ####

#load data from multiple imputated data set (MICE)
load(file='C:/Users/Nico/PowerFolders/Supervision_SaraMIND/descriptive_data/MiND_MIdata.Rdata')
df_quest<-all.new.data[[41]] #last data frame - i.e. latest iteration in MICE
#get ID labels from final data frame
load(file='C:/Users/Nico/PowerFolders/Supervision_SaraMIND/descriptive_data/MiNDfinalsample.RData')
id_label<-data$id.label
rm(data)
id_label<-droplevels(id_label)
df_quest<-cbind(id_label,df_quest)
df_quest<-df_quest[,-3] #exclude sex variable as also in df_data

df_data<-merge(df_data,df_quest,by.x='id.label',by.y='id_label')

  #excluded participant (n=91 --> n=90) difference from attention paper to fast paper
  #id_label[!(id_label %in% df_describe$id.label)]
  
length(table(droplevels(df_data$id.label))) #after questionnaire data merging: n=90

# E. MERGE WITH LOGFILES  #### 
## - read report files ####
data.files<-list.files(path=getwd(), full.names=TRUE,recursive=T)
path_logfile<-data.files[grep('FastTask/report',data.files)]

function_read_logfile<-function(x){
  k<-try(read.table(x,skip=24+5,fill=T,sep=';',header=F,dec=',')) #also skip the first 5 trials, thus +5
  if(!inherits(k,'try-error')) k} #gives an error put not problematic

list_logfile<-lapply(path_logfile,function_read_logfile)

#read ids from logfile data
function_read_id<-function(x){
  k<-readLines(con=x, '5')
  pt.id<-k[grep('PatientId:',k)]
  pt.id<-substring(pt.id, 12,nchar(pt.id))
  pt.id<-sub('\"',"", pt.id) #saras code erased the last character of the patient id - corrected
  return(pt.id)}
names_logfiles<-as.character(lapply(path_logfile,function_read_id))
names(list_logfile)<-names_logfiles #save patient IDs to loaded data

## - CLEAN Logfiles (erase test data and duplicates) ####
#remove unlabelled data
list_logfile<-list_logfile[-which(names(list_logfile)=='00000000-0000-0000-0000-000000000000')]
list_logfile<-list_logfile[-which(names(list_logfile)=='character(0)')]

#remove duplicates - participants with more than one logfile
    #duplicated_ids<-names(list_logfile)[duplicated(names(list_logfile))]
    #duplicated_ids #-->these participants have more than one logfile
    #sapply(list_logfile,nrow)[which(names(list_logfile) %in% duplicated_ids)]
# --> '3272b2e4-db3a-4b3d-a71b-a60b929ae585' did the battery twice
# --> '4a7ed1c6-237d-4e56-9bb0-3872f18408d0' did parts of it twice
list_logfile<-list_logfile[-c(excluded_session1,excluded_session2)]


## - form data.frame and create id variable ####
id.length<-unlist(sapply(list_logfile,nrow))
id.label<-rep(names(id.length),id.length)
df_logfile<-rbind.fill(list_logfile)
df_logfile<-cbind(id.label,df_logfile)
rm(id.length,id.label)

## - label logfile ####

    ## DESCRIPTION OF LOGFILE VARIABLES according to DESCRIPTION FILE (see link above)
    # - game type ["fastask"]: always is shown as " fastask " 
    # - trial[1,80]: Explains the current trial for baseline 
    # - trial block[1,80]: The same as trial. 
    # - target index[0,1,2,3]: Explains which circle turns in yellow color. 
    # - choseIndex[0,1,2,3]: Explains which circle patient chosed as a response. 
    # - posCircle0[vec2]: Explains the position for the circle 1 in pixels. 
    # - posCircle1[vec2]: Explains the position for the circle 2 in pixels. 
    # - posCircle2[vec2]: Explains the position for the circle 3 in pixels. 
    # - posCircle3[vec2]: Explains the position for the circle 4 in pixels. 
    # - averageTime[number]: Explains the average of the time of correct reponses (got in baseline part).  
    # - reactionTime[number]:  Explains response time in ms. 
    # - correctReponse[0,1]: Explains whether patient did a correct response or not.  

names(df_logfile)<-c('id.label','game_type','trial','trial_block','target_index','chosen_index','pos_circle0','pos_circle1','pos_circle2','pos_circle3','pos_circlex','average_time','reaction_time','correct_response')
#--> names taken from task description: file:///C:/Users/Nico/PowerFolders/Supervision_SaraMIND/Game%20log%20descriptions%20MIND_Jordi.pdf
#df_logfile[df_logfile==-1]<-NA #--> not necessary

## - CHECK logfile file data ####
attach(df_logfile)
table(target_index,chosen_index) #--> overall performance
hist(reaction_time[reaction_time<1000])
with(df_logfile[df_logfile$reaction_time<1000,],by(reaction_time,game_type,summary)) #--> differences in reaction time by game type
table(correct_response,game_type) ##--> descriptively higher error rate in fastask

hist(average_time) #--> not usable often just 0
table(pos_circlex) #--> circle variables conain no information

table(trial) ##--> trial 5 occurs 20 more than expected
  #trial_per_id<-table(id.label,trial)
  #trial_per_id[,1:5]
  
table(game_type)
length(table(id.label))

detach(df_logfile)

## - REMOVE unnecessary logfile data/variables ####
df_logfile<-df_logfile[,c('id.label','game_type','trial','target_index','chosen_index','average_time','reaction_time','correct_response')]

# - REMOVE 5 first trials of baseline in logfiles - they have been done twice
#--> we should expect 152 entries per participant (baseline=72 + fastask=80), but have 157 per participant
# with(df_logfile[df_logfile$trial<=20,],hist(trial))
# with(df_logfile[df_logfile$trial<=20 & df_logfile$game_type=='fastask',],hist(trial))
# with(df_logfile[df_logfile$trial<=20 & df_logfile$game_type=='baseline',],hist(trial))
#--> during baseline, the first four trials occur and especially trial 1 occurs more often
### --> NOW ALREADY EXCLUDED DURING READING OF LOGFILES

## - CREATE row identifier by interaction of id.label and currentTrial in log and report data ####

  #clean variabnles needed for ID MERGING VARAIBLES below
  df_data$typeGame<-droplevels(df_data$typeGame) #drop unused level
  df_data$id.label<-droplevels(df_data$id.label)
  
  id_data<-with(df_data,interaction(id.label,typeGame,currentTrial)) #trial identifier ID x type x trial
  id_logfile<-with(df_logfile,interaction(id.label,game_type,trial)) #trial identifier ID x type x trial

  df_data<-cbind(df_data,id_data)
  df_logfile<-cbind(df_logfile,id_logfile)

## - CALCULATE HIT/OMISSION/COMISSION errors ####
  
  with(df_logfile,table(chosen_index,target_index))
  
  hit<-with(df_logfile,ifelse(target_index==chosen_index,T,F))
  omission<-with(df_logfile,ifelse(chosen_index==-1,T,F))
  comission<-with(df_logfile,ifelse(chosen_index!=-1 & chosen_index!=target_index,T,F))
  
  table(hit)[2]/length(hit) #94% hits
  table(omission)[2]/length(omission) #0.2% omission errors
  table(comission)[2]/length(comission) #5.7% comission errors
  
  df_logfile<-data.frame(df_logfile,hit,omission,comission)
    
## - MERGE to df ####
  
  df<-merge(df_data,df_logfile,by.x='id_data',by.y='id_logfile')
    
# F. PUPIL DILATION PREPROCESSING ####
  
#reduce data as split takes very long  
df_pd<-df[,c('ts.trial','timestamp','leftEyePupilSize','rightEyePupilSize','leftValidity','rightValidity')]  
  
df_split<-split(df_pd,df$id_data) #takes around 90seconds to split
df_split<-lapply(df_split,function(x){x[order(x$ts.trial),]}) #order by ts.trial

  #timestamp differences
  mean_timestamp_diff<-round(sapply(df_split,function(x){mean(diff(x$timestamp),na.rm=T)}))/1000
  hist(mean_timestamp_diff[which(mean_timestamp_diff>0 & mean_timestamp_diff<70)],100,xlab='mean timestamp difference (ms)')
 
require(zoo)   
func_PD_basicpreprocessing<-function(x){
#define variables
Left_Diameter<-x$leftEyePupilSize
Right_Diameter<-x$rightEyePupilSize
Left_Validity<-x$leftValidity
Right_Validity<-x$rightValidity
RemoteTime<-x$timestamp
#plausible values
pl <- ifelse((Left_Diameter<2|Left_Diameter>8), NA, Left_Diameter)  
pl <- ifelse(Left_Validity<=2,pl,NA)
pr <- ifelse((Right_Diameter<2|Right_Diameter>8), NA, Right_Diameter)  
pr <- ifelse(Right_Validity<=2,pr,NA)
#linear interpolation
#Left
pl.smooth<-na.approx(pl,na.rm=F,rule=2,maxgap = 6) #impute missing values with interpolation
#Right
pr.smooth<-na.approx(pr,na.rm=F,rule=2,maxgap = 6) #impute missing values with interpolation
#take offset between left and right into account
pd.offset<-pl-pr  
pd.offset<-na.approx(pd.offset,rule=2)
#mean pupil dilation across both eyes
pl <- ifelse(is.na(pl.smooth)==FALSE, pl.smooth, pr.smooth+pd.offset)  
pr <- ifelse(is.na(pr.smooth)==FALSE, pr.smooth, pl.smooth-pd.offset)  
pd <- (pl+pr)/2  
return(pd)}
pd_list<-lapply(df_split,func_PD_basicpreprocessing)

  #see different mean trial length for baseline and fastask
  trial.length<-sapply(pd_list,length)
  hist(trial.length[trial.length<500],100)

  pd<-do.call(c,pd_list)
  df<-df[order(df$id_data,df$ts.trial),] #order for merging
  df<-data.frame(df,pd)
  

  #baseline pupil dilation - correct for mean in the first 200ms
    #pd_tonic<-sapply(pd_list,function(x){mean(x[1:6],na.rm=T)}) #first 200 ms <-- tonic corrected - far more missings
    #pd_tonic<-rep(pd_tonic,times=trial.length)
  pd_trial_baseline<-sapply(pd_list,function(x){mean(x[1:6],na.rm=T)}) #first 200 ms <-- tonic corrected - far more missings
  pd_trial_baseline<-rep(pd_baseline,times=trial.length)
  
  #phasic pupil dilation (in relation to trial)
    #rpd<-pd-pd_tonic
  rpd<-pd-pd_trial_baseline
    rpd[rpd>2]<-NA
    rpd[rpd<(-2)]<-NA
    hist(rpd,100) 
  df<-data.frame(df,rpd)
  
  #tonic pupil dilation (in relation to global value)
  pd_global_baseline<-with(df,round(mean(pd,na.rm=T),3))
  tpd<-with(df,pd-pd_global_baseline)
    tpd[tpd<(-2)]<-NA
    tpd[tpd>2]<-NA
    hist(tpd,100)
  df<-data.frame(df,tpd)  
  
# G. GAZE BEHAVIOR CALCULATION####
  hist(df$rightEyeRawX)
  hist(df$rightEyeRawY)
  hist(df$leftEyeRawX)
  hist(df$leftEyeRawY)
  
  attach(df)
  rx <- ifelse((rightEyeRawX<0|rightEyeRawX>1), NA, rightEyeRawX)  
  ry <- ifelse((rightEyeRawY<0|rightEyeRawY>1), NA, rightEyeRawY)  
  lx <- ifelse((leftEyeRawX<0|leftEyeRawX>1), NA, leftEyeRawX)  
  ly <- ifelse((leftEyeRawY<0|leftEyeRawY>1), NA, leftEyeRawY)  
  detach(df)
  rx <- ifelse(is.na(rx)==F, rx, lx)  
  lx <- ifelse(is.na(lx)==F, lx, rx)  
  ry <- ifelse(is.na(ry)==F, ry, ly)  
  ly <- ifelse(is.na(ly)==F, ly, ry)  
  gazepos.x <- (rx+lx)/2  
  gazepos.y <- (ry+ly)/2  
  
  
  g<-ggplot(data.frame(gazepos.x,gazepos.y)[df$strState=='foreperiod' & df$typeGame=='fastask',],aes(x=gazepos.x,y=gazepos.y))
  g+geom_hex(bins=30,aes(fill=..density..))
  
  g<-ggplot(data.frame(gazepos.x,gazepos.y)[df$strState=='stimuli' & df$typeGame=='fastask',],aes(x=gazepos.x,y=gazepos.y))
  g+geom_hex(bins=30,aes(fill=..density..))
  
  g<-ggplot(data.frame(gazepos.x,gazepos.y)[df$strState=='isi' & df$typeGame=='fastask',],aes(x=gazepos.x,y=gazepos.y))
  g+geom_hex(bins=30,aes(fill=..density..))
  
  #- calculate center deviation that can be applied to control for gaze behavior
  center_dev<-sqrt(abs(0.5-gazepos.x)^2+abs(0.5-gazepos.y)^2)
  df<-data.frame(df,center_dev)
  
    hist(center_dev[df$strState=='foreperiod' & df$typeGame=='fastask'])
    hist(center_dev[df$strState=='stimuli' & df$typeGame=='fastask']) #participant look at the stimuli
    hist(center_dev[df$strState=='isi' & df$typeGame=='fastask'])
    
  #gaze valid --> check this variable!
  gaze_valid<-ifelse((center_dev<0.1 & df$strState=='foreperiod')|(center_dev<0.1 & df$strState=='isi')|(df$strState=='stimuli'),T,F)  
  gaze_valid<-ifelse(is.na(gaze_valid),F,gaze_valid)
  df<-data.frame(df,gaze_valid,gazepos.x,gazepos.y)

# H. DATA CLEANING / RELABELLING ####
  
  # RELABEL group variable (group variable is derived from MICE data)
  df$groupASD<-as.factor(ifelse(df$groupASD==1,'ASD','no_ASD'))
  df$groupADHD<-as.factor(ifelse(df$groupADHD==1,'ADHD','no_ADHD'))
  
  #CHANGE SEX VARIABLE
  df$sex<-as.factor(ifelse(df$sex=='f','f','m'))
  
  #droplevels of id.label
  id.label<-df$id.label.x
  df<-subset(df, select=-c(id.label.x,id.label.y))
  df<-data.frame(id.label,df)
  
  #droplevels game state
  df$strState<-droplevels(df$strState)
  
  #remove long trials
  hist(sapply(pd_list,length),100)
  valid_trial_length<-which(sapply(pd_list,length)<400)
  df<-df[df$id_data %in% names(valid_trial_length),]
  
## --> save relevant object (df) ####
  
  #save(df,file='C:/Users/Nico/PowerFolders/Paper_MINDdata_fasttask/df_200919.Rdata')
  #--> save with R 3.4 processed - no embedded null warnings for file reads
  #save(df,file='C:/Users/Nico/PowerFolders/Paper_MINDdata_fasttask/df_230919.Rdata')
  #--> save with R 3.6.1 processed - 'embedded null' warnings and workaround for file reads
  # seems to drop certain participants
  #save(df,file='C:/Users/Nico/PowerFolders/Paper_MINDdata_fasttask/df_240919.Rdata')
  #--> save with R 3.4 processed - no embedded null warnings for file reads
  save(list=c('df','df_describe','df_agg'),file='C:/Users/Nico/PowerFolders/Paper_MINDdata_fasttask/df_021019.Rdata')
  #--> save with R 3.4 processed - no embedded null warnings for file reads
  