#setup ####

#NOTE: MAKE sure that lmerTEst is not loaded as it is not compatible with simr

#for R<3.4 - install archieved version of plotrix, which is compatible with R<3.5
# install_version("plotrix", version = "3.6", repos = "http://cran.us.r-project.org")
# --> plotrix is required by simr
#install.packages('simr')

require(lme4)
require(simr)

citation('simr')

#LOAD data ####

setwd('C:/Users/Nico/PowerFolders/Paper_MINDdata_fasttask')
#load('df_021019.Rdata')
load('df_261020.Rdata')


#CREATE per-trial data set - see also associated markdown file ####

#selected interaction:
int_trial<-with(df_valid,interaction(id.label,currentTrial,typeGame))

factors_trial<-with(df_valid,aggregate(data.frame(groupASD,groupADHD,id.label,currentTrial,typeGame,age,sex,pIQ),
                                       by=list(int_trial),head,n=1))

metrics_mean_trial<-with(df_valid,aggregate(data.frame(tpd,rpd,reaction_time),
                                            by=list(int_trial),mean,na.rm=T))

metrics_sd_trial<-with(df_valid,aggregate(data.frame(tpd,rpd), 
                                          by=list(int_trial),sd,na.rm=T))
#note: on trial level no variation of reaction time can be calculated

#-- state specific data (strState) ###
metrics_fore_mean_trial<-with(df_valid[df_valid$strState=='foreperiod',],aggregate(data.frame(tpd,rpd),
                                                                                   by=list(interaction(id.label,currentTrial,typeGame)),mean,na.rm=T))

metrics_stim_mean_trial<-with(df_valid[df_valid$strState=='stimuli',],aggregate(data.frame(tpd,rpd),
                                                                                by=list(interaction(id.label,currentTrial,typeGame)),mean,na.rm=T))

metrics_isi_mean_trial<-with(df_valid[df_valid$strState=='isi',],aggregate(data.frame(tpd,rpd),
                                                                           by=list(interaction(id.label,currentTrial,typeGame)),mean,na.rm=T))

#relabel
names(metrics_fore_mean_trial)[2:3]<-c('tpd_fore','rpd_fore')
names(metrics_stim_mean_trial)[2:3]<-c('tpd_stim','rpd_stim')
names(metrics_isi_mean_trial)[2:3]<-c('tpd_isi','rpd_isi')
names(metrics_mean_trial)[2:4]<-c('tpd_m','rpd_m','rt_m')
names(metrics_sd_trial)[2:3]<-c('tpd_sd','rpd_sd')

###MERGE several data.frame at once (Reduce, which is very fast)
df_trial<-Reduce(function(x,y){merge(x = x, y = y, by = "Group.1", all.x=T)}, 
                 list(factors_trial,metrics_mean_trial,metrics_sd_trial, 
                      metrics_fore_mean_trial, metrics_stim_mean_trial, metrics_isi_mean_trial))

df_trial$group3<-with(df_trial,ifelse(groupASD=='ASD','ASD',
                                      ifelse(groupADHD=='ADHD','ADHD','TD')))


##create phasic response variable on trial level
rpd_phasic<-with(df_trial,rpd_stim-rpd_fore)
df_trial<-data.frame(df_trial,rpd_phasic)


#NOTE: TPD - global median corrected PD values, RPD - first 200ms of trial corrected values
#TPD, RPD do not matter for the difference measure rpd_phasic, but TPD needs to be chosen
# for BPS as otherwise, it would also be a change measure

#with(df_trial,hist(rpd_phasic,100))
#with(df_trial,hist(rt_m,50))

#change typeGame levels to meaningful labels - for facet labelling in plotting
levels(df_trial$typeGame)<-c('low-utility','high-utility')


#DEFINE MODELS - that random variance is taken for power analysis ####

    #FAST TASK STUDY - model:
      #-->make sure, there is no missing data, otherwise simulation will give errors --> and POWER==
      df_power<-na.omit(df_trial)
      
      #NOTE: we removed random effect currentTrial due to low variance and excluded several interactions, as these seem to overload the model
      power.lmm<-lme4::lmer(scale(rt_m)~typeGame+group3+
                              rpd_fore+
                              rpd_phasic+
                              typeGame:group3+
                              (1|id.label),data=df_power,REML=F)
      
                                      #EXCLUDED EFFECTS
                                        # typeGame:scale(rpd_fore)+
                                        # typeGame:scale(rpd_phasic)+
                                        # groupASD:scale(rpd_fore)+
                                        # groupASD:scale(rpd_phasic)+
                                        # groupADHD:scale(rpd_fore)+
                                        # groupADHD:scale(rpd_phasic)+
                                        # typeGame:groupASD:scale(rpd_fore)+
                                        # typeGame:groupADHD:scale(rpd_fore)+
                                        # typeGame:groupASD:scale(rpd_phasic)+
                                        # typeGame:groupADHD:scale(rpd_phasic)+
                                        # sex+scale(age)+scale(pIQ)
                            
      
      #compares to model without the interested effect
      compared_model<-~typeGame+group3+
        rpd_fore+
        rpd_phasic
        
                                    
#-->DEFINE EFFECT SIZE AND EFFECT that should be tested (group difference) ####      
      summary(power.lmm)
      VarCorr(power.lmm) #random effect variance
      fixef(power.lmm)
      
      #simulation 1
      #fixef(power.lmm)["group3ASD"]<-0.15 #define effect size to be tested
      
      #simulation 2 (rerun model estimation before defining effect to reset observed fixed effects)
      fixef(power.lmm)["typeGamehigh-utility:group3ASD"]<-0.10 #define effect size to be tested
      
#-->SIMULATE POWER #### 
      #---> fcompare is used as INTERESTED EFFECT is involved in interaction - fcompare compares to model without the interested fixed effect
doTest(power.lmm, test = fcompare(compared_model)) #test whether comparison works
powerSim(power.lmm, test = fcompare(compared_model), alpha=0.05, nsim=1000) #nsim == number of simulation
#powerSim(power.lmm,fixed('typeGamehigh-utility:groupASDno_ASD',method='lr')) #nsim == number of simulation
