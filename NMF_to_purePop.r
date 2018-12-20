# make the NMF calculation- using the ls-nmf- weight by population size
# find the pur population and make the ggplot
library(base)
library (NMF)
library(ggplot2)

###############################################
#          Part 1- NMF
##############################################

path3<-"C:/Users/Alexandra/Documents/GitHub/Haplotype-Admixture/Date/"
setwd(path3)
freq<-read.csv("freq_norm.csv", header=F,sep=',')
pop_name<-read.csv("Population_name.csv", header=F, sep=',')
haplo<-read.csv("Haplo_s.csv", header=F, sep=',')
haplo_Size<-as.matrix(read.csv("Population_size.csv", header=F, sep=','))

dev<-c(1.e5)
pop_size<-t(strtoi(haplo_Size))/1.e5 # get the populaton size and devide them by 1.e5
Z<-matrix(pop_size, nrow(freq),ncol(freq),byrow=TRUE )
Res<-NULL
w<-NULL
h<-NULL
pl=7 # number  rows in the h martrix of col in w matrix
Res <- nmf(freq,pl, 'ls-nmf', weight = Z)
#w<-basis(Res)
h<-coef(Res)
save(file='h_matrix.RDate',h)
###############################################
#          Part 2- Find the pure populations
##############################################

PopName_mat<-read.csv('Population_name_and_index_Sort.csv', header=F, sep=',' ) # sorted population with the index of the primer order
PopName<-as.character(PopName_mat[1:(nrow(PopName_mat)),2]) # get the sort pop name
PoP_Order<-PopName_mat[1:(nrow(PopName_mat)),1]


h_colSum<-colSums(h) 
norm_h<-matrix(nrow=nrow(h),ncol=ncol(h)) # mat a matrix  for norm resolts 
norm_h<-t((t(h))/h_colSum) # checked that the sum of each row is 1 -rows are the population and col is ranks
print(colSums(norm_h)) # check that summ s 1 
# change the col order of the h matrix
norm_h<-norm_h[,PoP_Order]

#--------------------------looking to the pure group------------------

# get the population that have the highest value in each rank  

# do only onse 
matrix_pop10_0.8<-as.matrix(norm_h,nrow=nrow(h),ncol=ncol(h))
index<-NULL
pop<-NULL
max_val<-NULL
for (j in 1:nrow(h)){
  index[j]<-which.max(norm_h[j,]) # find the max population from 66 in each pure col pop
  print(norm_h[j,index[j]]) 
  if  (norm_h[j,index[j]]>0.90)  #change the cut off 
  {
    pop[j]<-PopName[index[j]] # find the pur pop
  }
}
# define pur pop name 
# 1.E.Asian 1-6
# 2.S_Asian 7-13
# 3. Middle_East 34:35
# 4. Hispanic 51:56
# 5. Jews 14:33
# 6. African 47:50
# 7.Europe/other 36:46
# make a list of the population order of the population names
sub_group<-NULL
sub_group[[1]]<-c(1:8)# E_Asian 
sub_group[[2]]<-c(9:13)# S_Asian 
sub_group[[3]]<-c(14:31)# Jews 
sub_group[[4]]<-c(32:35)#  Middle_East 
sub_group[[5]]<-c(36:46)# Europe/other 
sub_group[[6]]<-c(47:50)# African 
sub_group[[7]]<-c(51:56)#   Hispanic

name_subG<-c('E_Asian','S_Asia','Jews','Middle_East','Europe','African','Hispanic')

# sort the lines in the norm_h matrix in the sub group order 
line_order<-order(index) #   the order of the line in the norm_h matrix
norm_h<-norm_h[line_order,]  # change the row order in matrix
pop<-pop[line_order] # chage the pure pop name order
index<-index[line_order] # change the index of the pur pop 

# Sort the pure populations by subgorup connection according to the index, counting the number of pure populations in each subgroup
num_pur_pop<-replicate(length(name_subG), list(0)) 
for (j in 1:length(sub_group)){
  for (jj in 1:nrow(h)){ 
    if (sum((index[jj]==sub_group[[j]])*1)){
      num_pur_pop[[j]]<-num_pur_pop[[j]]+1
    }
  }}
# convert  population index to population name
pop_order_all<-NULL
for (k in 1:length(sub_group)){
  pop_order_all[[k]]=PopName[sub_group[[k]]]}
#save( pop_order_all,file="pop_order_all.RData")
##===========================
# change the row name to genaric name to each pure population
row_name<-NULL
for (jj in 1:length(sub_group)){
  if (num_pur_pop[[jj]]>1){
    for (ii in 1:num_pur_pop[[jj]]){
      row_name<-c(row_name,paste0(name_subG[jj],'_',ii))
    }}
  else if(num_pur_pop[[jj]]==1){
    row_name<-c(row_name,name_subG[jj])
  }
}

colnames(norm_h)<-c(t(as.matrix(PopName))) # the col names must be a row vector 
row.names(norm_h)<-row_name

#-----------------------------------------
# convert the data to long for bar plot - preparing for ggplot using 
#-----------------------------------------

library (reshape2)
long_new_matrix<-melt(norm_h) # make a long format
colnames(long_new_matrix)<-c('rank','PopName','value')

# define the color that we whant to use
color_range<-c('#9966FF','#6633FF','#3300cc','#330099',#E.Asian
               '#FFFF00' ,'#FF9900','#FF6600','#FF3300', # S_Asian
               '#00FFFF','#33CCFF','#3399FF','#0066CC',# Jews
               '#cccccc','#999999','#666666','#333333',#,'#000000', # Middle_East
               '#99FF33','#00CC33','#009900','#006600',#Europe/other 
               '#CC6600','#993300','#663300','#330000', # African 
               '#ff6666','#ff0000','#CC0000','#990000')#Hispanic 
pie(rep(1,length(color_range)),col=color_range) # colore plots 


#------------------------------------------------------
# color to SQRT

# use the pop_Rank to change the colore in color_range _ we have more color than we need
# change the color vector  for the new levels/factor oreder 
#length(pop_Rank[[1]])# find the number of allemnts in each group

cr_subGroup<-NULL # define the color to each sub group otomticly
k=1
for (iii in 1:length(num_pur_pop)){
  if (num_pur_pop[[iii]]==0){
    k<-((k%/%4-k%%4%/%10+1)*4+1)
    next
  }
  else{
    for (ii in 1:num_pur_pop[[iii]]){
      cr_subGroup<-c(cr_subGroup,color_range[k])
      k=k+1}}
  if (k%%4==0){
    k<-k+1}
  else if (k>28){
    k<-k}
  else if (k%%4==1){
    k<-k}
  else {
    k<-((k%/%4-k%%4%/%10+1)*4+1)}
}
pie(rep(1,nrow(h)),col=cr_subGroup) # colore plots 
color_range1<-cr_subGroup

#  define the rank as a levels - by using factor 
levels(long_new_matrix$rank) 
factor(long_new_matrix$rank)
# change the rank order and the color ordet to-I want the color be in the same order in the bar
g<-levels(factor(long_new_matrix$rank)) # convert the rank to factor and then put in levels and this to a var g -vector


#bar plot
pdf(paste0("1B.pdf"),15, 8)
barplot<-ggplot(data=long_new_matrix,aes(x=PopName,y=value,fill=rank)) +
  geom_bar(colour='black',stat='identity')+  
  scale_fill_manual(values=color_range1)+ 
  ggtitle("56 Population sortingpure population")
barplot+ theme(axis.title.x = element_text( colour="black", size=14), 
               axis.text.x  = element_text(angle=45,hjust=1,color='black', size=14),
               axis.title.y=element_blank(),axis.text.y=element_text(size=14),
               legend.title=element_text(size=14),legend.text=element_text(size=14))
dev.off()


