###takes gated or non gated populati0ns - expects a split population on the primary detector depicted by !#####!.
#runs on all FCS files in directory. 


# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#    BiocManager::install("flowCore")
#    
#load packages

library(flowCore)

#SETUP DETECTORS CHANNELS HERE
iterate_list <- c(7:26,28:41) #enter fluorescent detector dentifiers here, discard scatters (normally 1-6) and any other detectors #ver9 configuration
#iterate_list <- c(3:25,27:40) #enter fluorescent detector dentifiers here, discard scatters (normally 1-6) and any other detectors #ver10 configuration (no width/height)
#iterate_list <- c(31:35) #enter fluorescent detector dentifiers here, discard scatters (normally 1-6) and any other detectors #ver10 configuration only violet laser...

# #set the working directory from which the files will be read from
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)
#    
#make a list of all .fcs files in folder
filelist <- list.files(pattern = ".fcs")
 
#placeholder for information.
output_list <- c("filename","seperating_index", "sensitivity_index", "abs_seperating_index","abs_sensitivity_index")

for (a in 1:length(filelist)){
#read in FCS file (will be pregated lymphocytes)
data <- read.FCS(filelist[a])

#Primary detector should be enclosed in symbol TBD - !***!
#get primary detector name from filename and search detectors
pri_det <- substr(filelist[a],(unlist(regexec("!",filelist[a])))+1,(unlist(regexec("!",filelist[a]))+7))
pri_det <- paste(pri_det,"-A",sep="")
#split population into two on primary detector using kmeans filtering
pops <- eval(parse(text = paste('split(data, kmeansFilter("',pri_det,'"=c("Negative","Positive"), filterId="myKMeans"))',sep='')))

#draw out necessary statistics - set which parameters to use
#primary seperating index = Pm-Nm/2XSD
stain_index <- ((median(pops$Positive@exprs[,pri_det]))-(median(pops$Negative@exprs[,pri_det])))/(2*sd(pops$Negative@exprs[,pri_det]))

#calculate sensitivity index
#Pm-Nm/((84N-Nm)/0.995)
seperation_index <-  as.numeric(((median(pops$Positive@exprs[,pri_det]))-(median(pops$Negative@exprs[,pri_det])))/
  (((quantile(pops$Negative@exprs[,pri_det],0.84))-(median(pops$Negative@exprs[,pri_det])))/0.995))
  
#calculate absolute stain index - sum all stain indexes (example below)
stain_index_L <- 0 #reset to zero
seperation_index_L <- 0
for (b in 1:length(iterate_list)){
  detector <- iterate_list[b]
stain_index_in_list <- ((median(pops$Positive@exprs[,detector]))-(median(pops$Negative@exprs[,detector])))/(2*sd(pops$Negative@exprs[,detector]))
stain_index_L <- stain_index_L+stain_index_in_list

seperation_index_in_list <-  as.numeric(((median(pops$Positive@exprs[,detector]))-(median(pops$Negative@exprs[,detector])))/
                                  (((quantile(pops$Negative@exprs[,detector],0.84))-(median(pops$Negative@exprs[,detector])))/0.995))
seperation_index_L<- seperation_index_L+seperation_index_in_list

}
print(paste(filelist[a]))
print(paste("Primary stain index = ",stain_index,sep=""))
print(paste("Primary seperating index = ",seperation_index,sep=""))
print(paste("Absolute stain index = ",stain_index_L,sep=""))
print(paste("Absolute seperation index = ",seperation_index_L,sep=""))

row_entry <- c(filelist[a],stain_index, seperation_index, stain_index_L,seperation_index_L)
output_list <- rbind(output_list,row_entry)


}

write.csv(output_list,"completed.csv")

