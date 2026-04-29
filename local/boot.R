library(ubc)
r<-c(rep(.4,20),rep(0,25))
names(r)<-paste0("x",1:45)
N<-200
N1<-round(N/2)
data<-simulate_correlated_y(200,r)

bigone<-function(data,datadis,vars) {

  discovery<-ubc_lm_estimates(formula = y~.,data = datadis,varnames = vars)
  confirmation<-ubc_lm_estimates(formula = y~.,data = data,varnames = vars)
 .data<-as.data.frame(cbind(y=confirmation$estimate,x=discovery$estimate))
 .weights<-1/confirmation$SE^2
 .loo_estimates(.data,loo,weights = .weights)

}

.loo_estimates<-function(data,loo=TRUE,weights=NULL) {
  if (loo) {
    est<-sapply(seq_len(nrow(data)), function(i) {
      mod<-stats::lm(y~x,data = data[-i,],weights = weights[-i])
      predict(mod,newdata = data)[i]
    })
  } else {
    mod<-stats::lm(y~x,data = data,weights = .weights)
    est<-predict(mod)
  }
  return(est)
}

bot.fun<-function(data,inds,model_weights,loo=TRUE) {
  .data<-data[inds,]
  .weights<-model_weights[inds]
  .loo_estimates(.data,loo=loo,weights = model_weights)
}
r<-boot::boot(.data,bot.fun,R=10,model_weights=.weights)
apply(r$t,2,sd)

data1<-data[1:(N1),]
data2<-data[(N1+1):N,]
discovery<-ubc::lm_candidates(y~.,data1)

bot.fun<-function(data,inds,datadis,vars) {
  .data1<-data[inds,]
  .data2<-datadis[inds,]
  bigone(.data2,.data1,vars)
}

r<-boot::boot(data2,bot.fun,R=100,datadis=data1,vars=discovery$var)
apply(r$t,2,sd,na.rm=T)

summary(ubc_lm(y~.,discovery,data2,se_method = "full"))
