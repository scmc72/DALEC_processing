# this pro is correct the in situ spectrum using method proposed by Jiang et al., 2020, ISPRS.
# Rrs400-Rrs850 must be available to use this method.
# 
#
# the flag of the corrected spectra: 
# 0, good data;
# 1, Oxygen aborption at 760nm is too strong, may influence the Rrs shape, temporal threshold is 0.0005, can adjust this to a more strict value
# 2, out of the training range of the method;
# 3, extremely turbid water, algae bloom, floating materials, high TSM etc;
# 4, the number of Rrs between 400-700 nm for delta corrected spectrum is higher than 50 wavelenths, temporal threshold is N=50, can adjust this to a more strict value
# 
#
# the single band Rrs750 model is used in this pro, but when calculate RHW and delta, median750-780,775-785,805-815,835-845 were used
# the input of this pro is the 1nm spectrum
# 
# Dalin Jiang, University of Stirling, UK
# email: dalin.jiang@stir.ac.uk
#
# June 20, 2021
#
#************************************************************************************************

time_start<-Sys.time()
#
Rrs_file<-"~/Documents/StirlingUni/Project/DALEC_processing/R-code/data/2022.csv"
Rrs_dt<-read.csv(Rrs_file,header=TRUE)
#
WAVE_MIN <- 400
WAVE_MAX <- 1000
Rrs_begin<-which(names(Rrs_dt)=="Rrs400")
Rrs_end<-which(names(Rrs_dt)=="Rrs1000")
#
Ref_dt<-Rrs_dt[,Rrs_begin:Rrs_end]
#

cor_jiang<-function(site_Rrs,wave_min,wave_max,wave_int){
	site_Rrs<-as.numeric(site_Rrs)
	names(site_Rrs)<-paste("Rrs",seq(wave_min,wave_max,wave_int),sep="")
	#
	pos_400<-which(names(site_Rrs)=="Rrs400")
	pos_490<-which(names(site_Rrs)=="Rrs490")
	pos_700<-which(names(site_Rrs)=="Rrs700")
	pos_750<-which(names(site_Rrs)=="Rrs750")
	pos_765<-which(names(site_Rrs)=="Rrs765")
	pos_780<-which(names(site_Rrs)=="Rrs780")
	pos_810<-which(names(site_Rrs)=="Rrs810")
	pos_840<-which(names(site_Rrs)=="Rrs840")
	#
	# -------RHW correction------
	md_750_780<-median(site_Rrs[pos_750:pos_780],na.rm=TRUE)
	Rrs780<-median(site_Rrs[(pos_780-5):(pos_780+5)],na.rm=TRUE)
	Rrs810<-median(site_Rrs[(pos_810-5):(pos_810+5)],na.rm=TRUE)
	Rrs840<-median(site_Rrs[(pos_840-5):(pos_840+5)],na.rm=TRUE)
	RHW<-Rrs810-Rrs780-(Rrs840-Rrs780)*(810.0-780.0)/(840.0-780.0)
	#
	est_md_750_780<-18267.884*RHW^3-129.158*RHW^2+3.072*RHW
	
	if (RHW > 0){
		delta<-md_750_780-est_md_750_780
	}else{
		delta<-md_750_780
	}
	cor_Rrs<-site_Rrs-delta	
	#
	# -------flag the quality--------
	cor_490<-median(cor_Rrs[(pos_490-5):(pos_490+5)],na.rm=TRUE)
	cor_750<-median(cor_Rrs[(pos_750-5):(pos_750+5)],na.rm=TRUE)
	med_765<-median(site_Rrs[(pos_765-3):(pos_765)],na.rm=TRUE)     # only 4 nm
	med_750<-median(site_Rrs[(pos_750-5):(pos_750+5)],na.rm=TRUE)
	med_780<-median(site_Rrs[(pos_780-5):(pos_780+5)],na.rm=TRUE)
	
	OAI<-abs(med_765-med_750-(med_780-med_750)*(765.0-750.0)/(780.0-750.0))   # Oxygen absorption peak/drop at 760nm
	N_neg<-length(which(cor_Rrs[pos_400:pos_700] < 0))
	#print(OAI)
	if(RHW > 0.0075){
		tmp_flag <- 2          # out of the training range of the method
	} else if(cor_750 > cor_490 & cor_750 > 0.01){
		tmp_flag <- 3          # extremely turbid, algae bloom, floating materials, high TSM waters
	} else if(N_neg > 50){
		tmp_flag <- 4          # negative results between 400-700 nm after delta correction 
	} else if(OAI > 0.0005){
		tmp_flag <- 1          # Oxygen aborption at 760nm is too strong, may influence the Rrs shape
	}else {
		tmp_flag <- 0          # good data
	}
	
	cor_result<-data.frame(DeltaCor_method="Jiang_2020",DeltaCor_flag=tmp_flag,RHW=RHW,Delta=delta,t(cor_Rrs))
	return(cor_result)
}


# do the correction
tmp_cor<-apply(Ref_dt,1,cor_jiang,wave_min=WAVE_MIN,wave_max=WAVE_MAX,wave_int=1)    
# tranform data format
tmp_result<-do.call("rbind",tmp_cor)
# conbine data information
final_result<-cbind(Rrs_dt[,2:(Rrs_begin-1)],tmp_result)
#

# output result
#
out_name<-gsub(".csv","_DeltaCor_Jiang.csv",Rrs_file) 
write.csv(final_result,out_name)
#
print(Sys.time()-time_start)




#-------------plot one example-----------
m<-19
xx<-seq(WAVE_MIN,WAVE_MAX)
plot(xx,as.numeric(Ref_dt[m,]),type="l",ylim=c(-0.005,0.03))                                 #,xlim=c(700,850)
lines(xx,as.numeric(tmp_result[m,paste("Rrs",xx,sep="")]),col="red")
lines(xx,rep(0,length(xx)),col="green")


# free up workspace and clear memory
rm(list=ls())
gc()
