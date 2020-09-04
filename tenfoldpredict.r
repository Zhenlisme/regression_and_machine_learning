library(ROCR)
data=read.table(file="C:/Users/luozh/Desktop/Gallus.feature.xls",header=T,sep="\t")
length(data$ID)
modelingdata=function(data,index){
  attach(data)
  vlabel = as.factor(Label[-index])
  vmrna = mRNA[-index]
  vorf = ORF[-index]
  vfickett = Fickett[-index]
  vhexamer = Hexamer[-index]
  mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
  test <- data.frame(vmrna = mRNA[index], vorf = ORF[index],vfickett = Fickett[index], vhexamer = Hexamer[index], vlabel=Label[index])
  test$prob <- predict(mylogit,newdata=test,type="response")
  output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
  rownames(output)=ID[index]
  detach(data)
  return(output)
}

tenfoldpredict=function(data){
  shuffled_data=c()
  for(i in sample(length(data$ID),length(data$ID))){
    shuffled_data=c(shuffled_data,as.vector(as.matrix(data[i,])))
  }
  shuffled_data=as.data.frame(matrix(shuffled_data,ncol = 6,byrow = T),stringsAsFactors = F)
  colnames(shuffled_data)=names(data)
  data=shuffled_data
  data$mRNA=as.numeric(data$mRNA)
  data$ORF=as.numeric(data$ORF)
  data$Fickett=as.numeric(data$Fickett)
  data$Hexamer=as.numeric(data$Hexamer)
  data$Label=as.numeric(data$Label)
  Response_list=list()
  Labls_list=list()
  for(i in seq(1,length(data$ID),length(data$ID)/10)){
    index=seq(i,i+length(data$ID)/10-1)
    mymodel=modelingdata(data=data,index = index)
    Response_list=c(Response_list,list(as.vector(mymodel[,"Prob"])))
    Labls_list=c(Labls_list,list(as.vector(mymodel[,"Label"])))
  }
  ROCR_data = list(predictions=Response_list,Labels=Labls_list)
  pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
  return(pred)
}

pred=tenfoldpredict(data = data)

perf <- performance(pred,"tpr","fpr")

par(mfrow=c(2,2),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)
plot(perf,col="blue",lty=3,xlab="1-Specificity",ylab="Sensitivity",ylim=c(0.7,1),xlim=c(0,0.3),main="",cex.axis=1.5,cex.label=1.5)	#AUC = 0.9927 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",xlab="1-specificity",ylab="sensitivity",main="",cex.axis=1.2,cex.label=1.2) 
abline(v=0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.05,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)

d=performance(pred,measure="prec", x.measure="rec")
plot(d,col="blue",lty=3,xlab="Recall (TPR)",ylab="Precision (PPV)",xlim=c(0.7,1),ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2)
plot(d,lwd=2,avg="vertical",col="red",xlab="Recall (TPR)",ylab="Precision (PPV)",add=T,cex.axis=1.2,cex.label=1.2)
abline(v=1.0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.95,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)

perf <- performance(pred,"acc")
plot(perf,col="blue",lty=3,xlab="Coding probability cutoff",ylab="Accuracy",ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2) 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",cex.axis=1.2,cex.label=1.2) 

abline(h=1,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)

S <- performance(pred,measure="sens")
P <- performance(pred,measure="spec")
plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(0.8,1),cex.axis=1.2,cex.label=1.2) 
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue") 
plot(P,col="red",lty=3, add=TRUE,) 
plot(P,lwd=2,avg="vertical",add=TRUE,col="red") 
legend(0.4,0.85,col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))
abline(v=0.55,lty="dashed",lwd=0.5)


dev.off()
