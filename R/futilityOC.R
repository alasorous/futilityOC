futilityOC <- function (txtype = "norm", pow.type = "PREDICTIVE", p.stop = 0.15, ind = 1, epts.r=1, inF = 0.5, H1 = TV, TV = H1,
	                      benef = NA, ctr.r = 0.5, NB.k = 0.5, over.dp = 1.5, alpha = 0.1, power = 0.8, n = 0, mat = 1, sd = 1,
	                      LRV = "", delta = "", ran.r = 1, pr.tv = 0.1, pr.lrv = 0.8, print=T) {

  #' Calculates the operating characteristics of futility rules
  #'
  #' Args:
  #'   txtype: The type of treatment (tx) effect "norm", HR","NB","PR","QP","RD", "RR". Default="norm"
  #'   pow.type: The type of power for futility rule "PREDICTIVE" or "CONDITIONAL" power. Default="PREDICTIVE"
  #'   p.stop: power boundary to stop for futility. Default = 0.15
  #'   ind: value 1, 2 or 3. 1 = power for statistical significance, 2 = power to meet Go criteria,
  #'                       3 = power to meet better than NoGo boundary. Default = 1
  #'   epts.r: Assumed correlation between the endpoints used at the futility analysis and the final. Default = 1
  #'   inF: information fraction at interim. Default = 0.5
  #'   H1: The treatment effect size assumed under the alternative hypothesis. Default = TV
  #'   TV: Target value for Go No Go decision. Default = H1
  #'   benef: A logical value. TRUE indicates a larger effect is beneficial. FALSE indicates a small effect is beneficial.
  #'        Default values depend on the type of treatment effect:
  #'              = TRUE for effect types of norm and RD
  #'              = FALSE for effect types of HR, NB, PR, QP and RR
  #'   ctr.r: Assumed event rate for the control arm, which is necessary to be specified for effect types
  #'        of NB, PR, QP, RD and RR. Default = 0.5
  #'   NB.k: The shape parameter of negative binomial distribution, which is necessary when effect
  #'       type is NB. Default = 0.5.
  #'   over.dp: The over-dispersion estimate for treatment effect type QP. Default = 1.5
  #'   alpha: Alpha value that is used for determining sample size n, if n is not specified. Default = 0.1 one-sided.
  #'   power: Statistical power that is used for determining sample size if n is not specified. Default = 0.8.
  #'   n: The sample size of the trial. When it is not specified (=0), its value is estimated. Note: n represents the number of
  #'      events when treatment effect is in HR with time-to-event endpoint and the default mat = 1 is not altered.
  #'   mat: Data maturity rate of time-to-event endpoint. Default = 1
  #'   sd: The standard deviation specified for endpoints following a normal distribution. Default = 1
  #'   LRV: The lower reference value specified for its use in Go/NoGo decision rules. Default = 2/3 * TV
  #'   delta: A vector of assumed/hypothesized treatment effects for which the operating characteristics
  #'        of futility analysis are to be evaluated. Default = a vector of 4 with values = (TV, LRV, TV/4, 0).
  #'   ran.r: The randomization ratio. Default = 1
  #'   pr.tv: Standard GNG rules. Probability that the treatment effect is equal to or greater than TV. Default = 0.1
  #'   pr.lrv: Standard GNG rules. Probability that the treatment effect is equal to or greater than LRV. Default = 0.8
  #'
  #' Returns:
  #'   Design.specifications
  #'   Go.NoGo.parameters
  #'   Prob.stop.under.Hypotheses
  #'   Futility.stopping.boundary
  #'   Conditional.power.at.stopping.boundary
  #'   Predictive.power.at.stopping.boundary
  #'   Critical.values.for.Go.NoGo.based.on.Interim.data

  L.pr <- qnorm(c(pr.tv,pr.lrv))
  z.a <- qnorm(1-alpha)

	Type <- is.element(txtype, c("HR","NB","PR","QP","RR"))

	if (Type) {
	  i.b <- ifelse(is.na(benef), -1, 2*benef-1)
		H1a <- i.b*log(H1)
		TVa <- i.b*log(TV)
		if (LRV=="") {
		  LRVa <- 2*TVa/3
		  LRV <- exp(i.b*LRVa)
		} else
		  LRVa <- i.b*log(LRV)
		if (delta[1] == "")
		  delta <- exp(i.b*c(TVa,LRVa,TVa/4,0))
		  Del <- i.b*log(delta)
	}

	if (!Type) {
	  i.b <- ifelse(is.na(benef), 1, 2*benef-1)
		H1a <- i.b*H1
		TVa <- i.b*TV
		if (LRV=="") {
		  LRVa <- 2*TVa/3
		  LRV <- i.b*LRVa
		} else
		  LRVa <- i.b*LRV
		if (delta[1]=="")
		  delta <- i.b*c(TVa,LRVa,TVa/4,0)
		  Del <- i.b*delta
	}

	m.D <- length(Del)
	tv.lrv <- c(TVa,LRVa)

	if (txtype=="norm") {
	  ss <- (ran.r+1)/ran.r*sd^2*rep(1,m.D)
	  sst <- ss[1]*c(1,1)
		D.scale <- "norm: mean diff"
	}
	if (txtype=="HR") {
	  ss <- (ran.r+1)/ran.r*1/mat*rep(1,m.D)
	  sst <- ss[1]*c(1,1)
		D.scale <- "HR: hazard ratio"
	}
	if (txtype=="NB") {
	  ss <- (1/(delta*ctr.r)+NB.k)/ran.r+(1/ctr.r)+NB.k
		sst <- (1/(c(1,H1)*ctr.r)+NB.k)/ran.r+(1/ctr.r)+NB.k
		D.scale <- "NB: rate ratio"
	}
	if (txtype=="PR") {
	  ss <- 1/(delta*ctr.r)/ran.r+1/ctr.r
		sst <- 1/(c(1,H1)*ctr.r)/ran.r+1/ctr.r
		D.scale <- "PR: rate ratio"
	}
	if (txtype=="QP") {
	  ss <- over.dp*(1/(delta*ctr.r)/ran.r+1/ctr.r)
		sst <- over.dp*(1/(c(1,H1)*ctr.r)/ran.r+1/ctr.r)
		D.scale <- "QP: rate ratio"
	}
	if (txtype=="RC") {
	  ss <- (1-(delta+ctr.r))*(delta+ctr.r)/ran.r+(1-ctr.r)*ctr.r
		sst <- (1-(c(0,H1)+ctr.r))*(c(0,H1)+ctr.r)/ran.r+(1-ctr.r)*ctr.r
		ave.r <- ctr.r+H1*ran.r/(1+ran.r)
		sst[1] <- (ran.r+1)*(1-ave.r)*ave.r
		D.scale <- "RC: rate diff"
	}
	if (txtype=="RR") {
	  ss <- (1-delta*ctr.r)/(delta*ctr.r)/ran.r+(1-ctr.r)/ctr.r
		sst <- (1-c(1,H1)*ctr.r)/(c(1,H1)*ctr.r)/ran.r+(1-ctr.r)/ctr.r
		D.scale="RR: rate ratio"
	}

	ss0 <- sst[1]
	ss1 <- sst[2]

	m2 <- (z.a*sqrt(ss0)+qnorm(power)*sqrt(ss1))^2/H1a^2
	if (n==0) {
	  n <- 2*ceiling(m2*(ran.r+1)/2)                                  # calculate sample size if n not given
	} else {
	  m2 <- n/(ran.r+1);
	  power <- pnorm(sqrt(m2*H1a^2/ss1)-z.a)                           # calculate power if n given
	}

	se <- sqrt(ss/m2)
	se1 <- se/sqrt(inF)
	se.h1 <- sqrt(ss1/m2/inF)

	sigm <- array(c(se1,se*sqrt(epts.r),se*sqrt(epts.r),se)^2,c(m.D,2,2))

	gng.c <- rep(1,m.D)%*%t(tv.lrv)+se%*%t(L.pr)
	gng.c[,2] <- apply(gng.c,1,max)
	gng.F <- rep(1,m.D)%*%t(tv.lrv)+se1%*%t(L.pr)
	gng.F[,2] <- apply(gng.F,1,max)

	tv.lrv2 <- tv.lrv+sqrt(ss1/m2)*L.pr
	tv.lrv2[2] <- max(tv.lrv2)
	tv.lrv2 <- c(tv.lrv2,z.a*sqrt(ss0/m2))

	# Calculate futility stopping rule given choice of predictive or conditional power

	tc0 <- ifelse(pow.type=="PREDICTIVE", qnorm(p.stop,tv.lrv2,sqrt((ss1/m2)*(1-inF)/inF))[4-ind],
	                                    qnorm(p.stop,tv.lrv2,sqrt(ss1/m2*(1-inF)))[4-ind])

	Pr.stop.H1 <- pnorm(tc0,H1a,se.h1)
	Pr.false.stop <- pnorm(tc0,c(0,tv.lrv),se.h1)
	names(Pr.false.stop) <- c("H0","TV","LRV")

	tc <- qnorm(Pr.stop.H1,H1a,se.h1)
	Tab0 <- Tab1 <- Tab2 <- Tab5 <- Tab6 <- array(0,c(m.D,5))

	# Conditional power calculation at futility stopping boundary

	c.pw <- pnorm(inF*tc+(1-inF)*c(tc,Del),rep(1,m.D+1)%*%t(tv.lrv2),sqrt(ss1/m2*(1-inF)))
	c.pw <- cbind(c(ifelse(is.element(txtype,c("norm","RC")),tc,exp(i.b*tc)),delta),c.pw)
	c.pw[,2] <- 1-c.pw[,2]

	# Predictive power calculation at futility stopping boundary

	p.pw <- pnorm(tc,tv.lrv2,sqrt(ss1/m2*(1-inF)/inF))
	p.pw[1] <- 1-p.pw[1]
	p.pw <- matrix(p.pw,1,3)

	Tab1[,2] <- pnorm(gng.c[,1],Del,se )
	Tab1[,5] <- 1-Tab1[,2]
	Tab0[,2] <- pnorm(gng.F[,1],Del,se1)
	Tab0[,5] <- 1-Tab0[,2]
	Tab1[,4] <- 1-pnorm(gng.c[,2],Del,se )
	Tab1[,3] <- Tab1[,5]-Tab1[,4]
	Tab0[,4] <- 1-pnorm(gng.F[,2],Del,se1)
	Tab0[,3] <- Tab0[,5]-Tab0[,4]

	for (j in 1:m.D)	{
		Tab2[j,3] <- mvtnorm::pmvnorm(c(tc,gng.c[j,1]),Inf,Del[j],sigma=sigm[j,,])
		Tab2[j,4] <- mvtnorm::pmvnorm(c(tc,gng.c[j,2]),Inf,Del[j],sigma=sigm[j,,])
		Tab5[j,2] <- mvtnorm::pmvnorm(c(gng.F[j,2],gng.c[j,2]),Inf,Del[j],sigma=sigm[j,,])
		Tab6[j,2] <- mvtnorm::pmvnorm(c(tc,gng.c[j,2]),c(Inf,Inf),Del[j],sigma=sigm[j,,])
		Tab6[j,3] <- mvtnorm::pmvnorm(c(tc,  -Inf),c(Inf,gng.c[j,1]),Del[j],sigma=sigm[j,,])
		Tab6[j,4] <- mvtnorm::pmvnorm(c(-Inf,gng.c[j,2]),c(tc,Inf),Del[j],sigma=sigm[j,,])
		Tab6[j,5] <- mvtnorm::pmvnorm(c(-Inf,-Inf),c(tc,gng.c[j,1]),Del[j],sigma=sigm[j,,])
	}

	Tab0[,1] <- Tab1[,1] <- Tab2[,1] <- Tab5[,1] <- Tab6[,1] <- delta

	Tab2[,2] <- 1-Tab2[,3]
	Tab2[,3] <- Tab2[,3]-Tab2[,4]
	Tab5[,4] <- Tab0[,4]-Tab5[,2]
	Tab5[,5] <- Tab2[,4]-Tab5[,2]
	Tab5[,3] <- 1-apply(Tab5[,c(2,4,5)],1,sum)
	Tab1 <- round(Tab1,3)
	Tab2 <- round(Tab2,3)
	Tab0 <- round(Tab0,3)
	i.stp <- pnorm(tc,Del,se1)
	Tab3 <- round(cbind(delta,i.stp,1-i.stp),3)

	if (Type) {
	  gng.c <- exp(i.b*gng.c)
	  gng.F <- exp(i.b*gng.F)
	} else {
	  gng.c <- i.b*gng.c
	  gng.F <- i.b*gng.F
	}

	T.tit <- c(D.scale," Pr(NoGo)","Pr(Pause)","   Pr(Go)")
	colnames(Tab0) <- colnames(Tab1) <- colnames(Tab2) <- c(T.tit," ")
	T.tit <- c(D.scale,"  IA & F","     None","  IA only"," Final only")
	colnames(Tab5) <- T.tit
	T.tit <- c(D.scale,"Cont'd|Go","Cont'd|NG","Disc'd|Go","Disc'd|NG")
	colnames(Tab6) <- T.tit
	colnames(gng.F) <- colnames(gng.c) <- c("Pause|NoGO","Go|Pause")
	colnames(Tab3) <- c(D.scale,"Pr(Discont'n)","Pr(Cont'n)")
	rownames(c.pw) <- c("Current trend:","Assumed trend:",rep(" ",m.D-1))
	colnames(c.pw) <- c("future.trend","Pr(NoGo|*)\u2267","Pr(Go|*)\u2266","Cond.Power")
	colnames(p.pw) <- c(" Pred.Pr(NoGo|*)\u2267"," Pred.Pr(Go|*)\u2266"," Pred.Power")
	Tab4 <- Tab2
	Tab4[,2:5] <- Tab2[,2:5]-Tab1[,2:5]
	if (print == T) {
	  cat(paste("Design.specifications:\n   H1 =", H1, ", alpha =", alpha, ", power =", round(power, 3), ", n =", n, ", inF =", inF,
		        "\nGo.NoGo.parameters:\n   TV =", TV, ", LRV =", round(LRV, 2), ", pr.tv =", pr.tv, ", pr.lrv =", pr.lrv,
       ##	Specified.Futility.Criterion=cbind(Pr.stop.H1),
		        "\nProb.stop.under.Hypotheses:\n   H0 =", round(Pr.false.stop[1],3), ", TV =", round(Pr.false.stop[2],3), ", LRV =", round(Pr.false.stop[3],3),
            "\nFutility.stopping.boundary =", round(c.pw[1,1],3),
		     #   "\nConditional.power.at.stopping.boundary:\n   Conditional Power =", round(c.pw[1,c(4)],3), ", Pr(Go|*) =", round(c.pw[1,c(3)],3),
         #              ", Pr(No Go|*) =", round(c.pw[1,c(2)],3),
		        "\nPredictive.power.at.stopping.boundary:\n   Predictive Power =", round(p.pw[,c(3)],3), ", Pr(Go|*) =", round(p.pw[,c(2)],3),
		                  ", Pr(No Go|*) =", round(p.pw[,c(1)],3),
       ##	Prob.of.termination.at.FA=Tab3[,c(1,3,2)],
		        "\nCritical.values.for.Go.NoGo.based.on.Interim.data=\n   Go.Pause =", round(gng.F[1,c(2)],4), ", Pause.NoGo =", round(gng.F[1,c(1)],4)	))
	}
       ##	Decision.prob.based.on.Interim.data=Tab0[,c(1,4,3,2)],
       ## Critical.values.for.Go.NoGo.based.on.Final.data=round(gng.c[1,c(2,1)],4),
       ## Decision.prob.based.on.Final.without.FA=Tab1[,c(1,4,3,2)],
       ## Decision.prob.based.on.overall.trial.with.FA=Tab2[,c(1,4,3,2)],
       ## Decision.prob.change.due.to.FA=Tab4[,c(1,4,3,2)],
       ## Potential.final.GNG.prob.by.FA.decision=round(Tab6,3),
       ## Go.prob.across.IA.and.Final.analyses=round(Tab5[,c(1,2,4,5,3)],3)	)
	d.out <- data.frame(H1, alpha, power, n, inF, TV, LRV, pr.tv, pr.lrv, t(Pr.false.stop), stop.boundary = c.pw[1,1], #t(c.pw[1,c(4,3,2)]),
	                    t(p.pw[,c(3,2,1)])) #, t(gng.F[1,c(2,1)]))
	names(d.out) <- c("H1", "alpha", "power", "n", "inF", "TV", "LRV", "pr.tv", "pr.lrv",
	                  "pr.stop.if.H0", "pr.stop.if.TV", "pr.stop.if.LRV",
	                  "stop.boundary", "pred.power.at.boundary", "pr.GO.at.boundary", "pr.NOGO.at.boundary")
	return(d.out)
}
