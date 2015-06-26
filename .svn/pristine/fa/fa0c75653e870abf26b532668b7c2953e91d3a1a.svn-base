LogN <-
function(Input, InputSP, EmpiricalR, EmpiricalRSP, NumOfEachGroup, AlphaIn, BetaIn,  PIn, NoneZeroLength)
{
    #2 condition case (skip the loop then maybe run faster? Code multi condition cases later)

        #For each gene (m rows of Input---m genes)
        #Save each gene's F0, F1 for further likelihood calculation. 
  
        #Get F0 for EE
        F0Log=f0(Input,  AlphaIn, BetaIn, EmpiricalR, NumOfEachGroup, log=T)
        #Get F1 for DE
        F1Log=f1(InputSP[[1]], InputSP[[2]], AlphaIn, BetaIn, EmpiricalRSP[[1]],EmpiricalRSP[[2]], NumOfEachGroup, log=T)

        #Get z
		#Use data.list in logfunction
		F0LogMdf=F0Log+600
		F1LogMdf=F1Log+600
		F0Mdf=exp(F0LogMdf)
		F1Mdf=exp(F1LogMdf)

		z.list=PIn*F1Mdf/(PIn*F1Mdf+(1-PIn)*F0Mdf)
		zNaNName=names(z.list)[is.na(z.list)]
		zGood=which(!is.na(z.list))
		if(length(zGood)==0){
		#Min=min(min(F0Log[which(F0Log!=-Inf)]), 
		#	min(F1Log[which(F1Log!=-Inf)]))
		tmpMat=cbind(F0Log,F1Log)
		tmpMean=apply(tmpMat,1,mean)
		F0LogMdf=F0Log-tmpMean
		F1LogMdf=F1Log-tmpMean
		F0Mdf=exp(F0LogMdf)
		F1Mdf=exp(F1LogMdf)

		z.list=PIn*F1Mdf/(PIn*F1Mdf+(1-PIn)*F0Mdf)
		zNaNName=names(z.list)[is.na(z.list)]
		zGood=which(!is.na(z.list))
		
		}
		###Update P
        #PFromZ=sapply(1:NoneZeroLength,function(i) sum(z.list[[i]])/length(z.list[[i]]))
        PFromZ=sum(z.list[zGood])/length(z.list[zGood])
        F0Good=F0Log[zGood]
		F1Good=F1Log[zGood]
		### MLE Part ####
        # Since we dont wanna update p and Z in this step
        # Each Ng for one row
		
		NumGroupVector=rep(c(1:NoneZeroLength),NumOfEachGroup)
		
		NumGroupVector.zGood=NumGroupVector[zGood]
		NumOfEachGroup.zGood=tapply(NumGroupVector.zGood,NumGroupVector.zGood,length)

        StartValue=c(AlphaIn, BetaIn,PIn)
   		     
		Result<-optim(StartValue,Likefun,InputPool=list(InputSP[[1]][zGood,],InputSP[[2]][zGood,],Input[zGood,],z.list[zGood], NoneZeroLength,EmpiricalR[zGood, ],EmpiricalRSP[[1]][zGood,], EmpiricalRSP[[2]][zGood,], NumOfEachGroup.zGood))
        #LikeOutput=Likelihood( StartValue, Input , InputSP , PNEW.list, z.list)
		AlphaNew= Result$par[1]
		BetaNew=Result$par[2:(1+NoneZeroLength)]
        PNew=Result$par[2+NoneZeroLength]
		##
        Output=list(AlphaNew=AlphaNew,BetaNew=BetaNew,PNew=PNew,ZNew.list=z.list,PFromZ=PFromZ, zGood=zGood, zNaNName=zNaNName,F0Out=F0Good, F1Out=F1Good)
        Output
    }

