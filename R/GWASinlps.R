GWASinlps = function( y, x, family=c("normal","binomial"), method=c("rigorous","quick"), cor_xy=NULL, mmle_xy=NULL, prior=c("mom", "imom", "zellner", "horseshoe"), tau, priorDelta = modelbbprior(1,1), k0, m, rxx, nskip = 3, niter = 2000, verbose = FALSE, seed = NULL, tau.hs.method = "halfCauchy", sigma.hs.method = "Jeffreys" )
{
	if(family == "normal")
	{
		if(!exists("time", mode="integer")) time = Sys.time()
		if(!is.null(seed)) set.seed(seed)

		# Compute corr if not provided
		if(is.null(cor_xy)) cor_xy = arma_cor(x,y) #	find corr of all x's with y, if not ptovided
		names(cor_xy) = colnames(x)	

		# Run GWASinlps iteration
		varsleft = colnames(x)
		varsselected = "" 
		selected_iterwise = list()

		max_nocollect = 0   # number of times no variables show up in the hppm 
		iter = 0
		while(looprun(varsselected,varsleft,max_nocollect,m,nskip))
		{			 
			# enter ith iteration
			iter = iter + 1
			if(verbose) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "") 
						
			# run nlpsLM
			run = nlpsLM( y = y, x = x[ , varsleft], cor_xy=cor_xy, prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose, tau.hs.method = tau.hs.method, sigma.hs.method = sigma.hs.method )
			lastcollect = run$hppm 

			# add the selected vars to final output & regress them out. if no var is selected, note the instance.
			if(length(lastcollect) > 0)
			{
				selected_iterwise = c(selected_iterwise, list(lastcollect))
				varsselected = c(varsselected, lastcollect)  # include HPPM vars in selected set
				y = fastLm( y ~ x[ , lastcollect] )$ residuals  # regress out HPPM vars
			} else
			{
				selected_iterwise = c(selected_iterwise, '')
				max_nocollect = max_nocollect + 1
				if(verbose) {cat("***","nskip=",max_nocollect,"***","\n")}
			}
			
			# exclude unselected vars from further consideration
			remove = run$not.selected	
			varsleft = setdiff1( varsleft, c(lastcollect,remove) )	

			# compute contributions (MLEs with previously selected vars in model) for the rest of the vars
			cor_xy = arma_cor(x[,varsleft],y)
			names(cor_xy) = varsleft
		}

		# Final output vars
		varsfinal = varsselected[ 2 : min(length(varsselected),m) ]

		# Print result
		cat("=================================", "\n","Number of selected variables: ", length(varsfinal), "\n", "Time taken: ", round(difftime(Sys.time(), time, units = "mins"),2), " min", "\n",  "=================================", "\n", sep = "")  

		# Return
		return( list(selected=varsfinal, selected_iterwise=selected_iterwise) )

	} else
	#
	#
	#
	#
	#
	#
	#
	if(family == "binomial")
	{
		if(method == "rigorous")
		{
			if(!exists("time", mode="integer")) time = Sys.time()
			if(!is.null(seed)) set.seed(seed)

			# Check if inputs are sufficient and consistent
			if(prior=="horseshoe") stop("horseshoe prior is not available for binomial family.")

			# Compute mmle if not provided
			mode(x) = "double"	#fastglm needs design matrix in double mode
			if(is.null(mmle_xy)) mmle_xy = apply( x, 2, function(z) coef( fastglm(y=y, x=cbind(1,matrix(z,ncol=1)), family = binomial(link = "logit")) )[2] )
			
			# Run GWASinlps iterations
			varsleft = colnames(x)
			varsselected = "" 
			selected_iterwise = list()

			max_nocollect = 0   # number of times no variables show up in the hppm 
			iter = 0
			contrib_xy = mmle_xy
			while(looprun(varsselected,varsleft,max_nocollect,m,nskip))
			{			 
				# enter ith iteration
				iter = iter + 1
				if(verbose) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "") 

				# run nlpsGLM
				run = nlpsGLM( y = y, x = x[ , varsleft], mmle_xy = contrib_xy, prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose )
				lastcollect = run$hppm 

				# add the selected vars to final output. if no var is selected, note the instance.
				if(length(lastcollect) > 0)
				{
					selected_iterwise = c(selected_iterwise, list(lastcollect))
					varsselected = c(varsselected, lastcollect)  # include HPPM vars in selected set
				} else
				{
					selected_iterwise = c(selected_iterwise, '')
					max_nocollect = max_nocollect + 1
					if(verbose) {cat("***","nskip=",max_nocollect,"***","\n")}
				}
				
				# exclude unselected vars from further consideration
				remove = run$not.selected	
				varsleft = setdiff1( varsleft, c(lastcollect,remove) )		

				# compute contributions (MLEs with previously selected vars in model) for the rest of the vars
				contrib_xy = apply( x[,varsleft], 2, function(z) coef( fastglm(y=y, x=cbind(1,x[,varsselected[-1]],z), family = binomial(link = "logit")) )[length(lastcollect)+2] )
			}

			# Final output vars

			varsfinal = varsselected[ 2 : min(length(varsselected),m) ] 

			# Print result
			cat("=================================", "\n","Number of selected variables: ", length(varsfinal), "\n", "Time taken: ", round(difftime(Sys.time(), time, units = "mins"),2), " min", "\n",  "=================================", "\n", sep = "")  

			# Return
			return( list(selected=varsfinal, selected_iterwise=selected_iterwise) )
		} else
		#
		#
		if(method == "quick")
		{
			if(!exists("time", mode="integer")) time = Sys.time()
			if(!is.null(seed)) set.seed(seed)

			# Check if inputs are sufficient and consistent
			if(prior=="horseshoe") stop("horseshoe prior is not available for binomial family.")

			# Compute mmle if not provided
			mode(x) = "double"	#fastglm needs design matrix in double mode
			if(is.null(mmle_xy)) mmle_xy = apply( x, 2, function(z) coef( fastglm(y=y, x=cbind(1,matrix(z,ncol=1)), family = binomial(link = "logit")) )[2] )
			
			# Run GWASinlps iterations
			varsleft = colnames(x)
			varsselected = ""
			selected_iterwise = list()

			max_nocollect0 = 0   # number of times no variables show up in the hppm 
			iter = 0
			
			while(looprun(varsselected,varsleft,max_nocollect0,1,nskip))  #exit if at least one var is selected 
			{
				# enter ith iteration
				iter = iter + 1
				if(verbose) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "")

				# run nlpsGLM
				run = nlpsGLM( y = y, x = x[ , varsleft], mmle_xy = mmle_xy[varsleft], prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose )
				lastcollect = run$hppm
				remove = run$not.selected	
				varsleft = setdiff1( varsleft, c(lastcollect,remove) )

				# add the selected vars to final output & regress them out. if no var is selected, note the instance.
				if(length(lastcollect) > 0)
				{
					selected_iterwise = c(selected_iterwise, list(lastcollect))
					varsselected = c(varsselected, lastcollect)  # include HPPM vars in selected set
					y = residuals( fastglm(y=y, x=x[ , lastcollect, drop=F], family = binomial(link = "logit") ), type = "deviance")  # collect deviance residuals
				} else
				{
					selected_iterwise = c(selected_iterwise, '')
					max_nocollect0 = max_nocollect0 + 1
					if(verbose) {cat("***","nskip=",max_nocollect0,"***","\n")}
				}

				# exclude unselected vars from further consideration
				remove = run$not.selected	
				varsleft = setdiff1( varsleft, c(lastcollect,remove) )	

				# compute contributions (MLEs with previously selected vars in model) for the rest of the vars
				if(length(lastcollect) > 0)
				{
					cor_xy = arma_cor(x[,varsleft],y)
					names(cor_xy) = varsleft
				}					
			}

			max_nocollect = 0
			while(looprun(varsselected,varsleft,max_nocollect,m-length(lastcollect),nskip-max_nocollect0))
			{			 
				# enter ith iteration
				iter = iter + 1
				if(verbose) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "") 
							
				# run nlpsLM
				run = nlpsLM( y = y, x = x[ , varsleft], cor_xy=cor_xy, prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose, tau.hs.method = tau.hs.method, sigma.hs.method = sigma.hs.method )
				lastcollect = run$hppm 

				# add the selected vars to final output & regress them out. if no var is selected, note the instance.
				if(length(lastcollect) > 0)
				{
					selected_iterwise = c(selected_iterwise, list(lastcollect))
					varsselected = c(varsselected, lastcollect)  # include HPPM vars in selected set
					y = fastLm( y ~ x[ , lastcollect] )$ residuals  # regress out HPPM vars
				} else
				{
					selected_iterwise = c(selected_iterwise, '')
					max_nocollect = max_nocollect + 1
					if(verbose) {cat("***","nskip=",max_nocollect,"***","\n")}
				}
				
				# exclude unselected vars from further consideration
				remove = run$not.selected	
				varsleft = setdiff1( varsleft, c(varsselected,remove) )	

				# compute contributions (MLEs with previously selected vars in model) for the rest of the vars
				cor_xy = arma_cor(x[,varsleft],y)
				names(cor_xy) = varsleft
			}

			# Final output vars
			varsfinal = varsselected[ 2 : min(length(varsselected),m) ]

			# Print result
			cat("=================================", "\n","Number of selected variables: ", length(varsfinal), "\n", "Time taken: ", round(difftime(Sys.time(), time, units = "mins"),2), " min", "\n",  "=================================", "\n", sep = "")  

			# Return
			return( list(selected=varsfinal, selected_iterwise=selected_iterwise) )
		}
	}      
}

nlpsLM = function( y, x, cor_xy, prior = c("mom", "imom", "zellner", "horseshoe"), tau, priorDelta = modelbbprior(1,1), k0, rxx, niter = 2000, verbose = F, tau.hs.method = "halfCauchy", sigma.hs.method = "Jeffreys" )
{ 	
	k0 = min(k0,ncol(x)) #if x has only 1 snp, but k0=2, then just reset k0=1
	names_sorted_cor_xy = names( sort( abs(cor_xy), decreasing = T ) [1:k0] )  # find x's with top k0 cors

	hppm = list()

	names_xx_input_set = NULL

	for(i in 1:k0)  #take the i'th of top k0 x's
	{
		corr_xx = c(abs( arma_cor(x, x[,names_sorted_cor_xy[i]]) ))   # find cor of all x vars with xi and give names
		names(corr_xx) = colnames(x) 

		names_xx_thres = na.omit( names(corr_xx[which(corr_xx >= rxx)]) ) # take x's with corr > rxx with xi

		names_xx_input =  setdiff1(names_xx_thres, unlist(hppm))  # exclude x's previously chosen in hppm and set the rest as the input x's 

		names_xx_input_set = c( names_xx_input_set, names_xx_input )  # store the input x's

		if(verbose) cat( "j =",i, "\ninput :", names_xx_input, "\n" )

		if(length(names_xx_input) != 0  )  # if there is some input x
		{
			if(prior == "mom") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], priorCoef = momprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input[ which(bb $ postMode == 1) ] # collect the HPPM vars
			} else
			#
			#
			if(prior == "imom") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], priorCoef = imomprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input [ which(bb $ postMode == 1) ] # collect the HPPM vars
			} else
			#
			#
			if(prior == "zellner") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], priorCoef = zellnerprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input [ which(bb $ postMode == 1) ] # collect the HPPM vars
			} else
			# emom not yet implemented in modelSelection. Also, emomLM won't run with only one x variable.
			#
			#
			if(prior == "horseshoe") 
			{
				fit = horseshoe(y, x[, names_xx_input, drop=F], method.tau = tau.hs.method, method.sigma = sigma.hs.method)
				fitsel = HS.var.select(fit, y, "intervals")
				hppm[[i]] = names_xx_input [ which( fitsel == 1 ) ]
			}
			if(verbose) cat( "selected :", hppm[[i]], "\n")  # print the HPPM vars
		}
	}

	names_xx_not_selected = setdiff1(names_xx_input_set, unlist(hppm) )

	return( list( hppm = unlist(hppm), not.selected = names_xx_not_selected ) )
}

nlpsGLM = function( y, x, mmle_xy, prior = c("mom", "imom", "zellner"), tau, priorDelta = modelbbprior(1,1), k0, rxx, niter = 2000, verbose = F )
{ 
	k0 = min(k0,ncol(x)) #if x has only 1 snp, but k0=2, then just reset k0=1
	names_sorted_mmle_xy = names( sort( abs(mmle_xy), decreasing = T ) [1:k0] )  # find x's with top k0 cors

	hppm = list()

	names_xx_input_set = NULL

	for(i in 1:k0)  #take the i'th of top k0 x's
	{
		corr_xx = c(abs( arma_cor(x, x[,names_sorted_mmle_xy[i]]) ))   # find cor of all x vars with xi and give names
		names(corr_xx) = colnames(x) 

		names_xx_thres = na.omit( names(corr_xx[which(corr_xx >= rxx)]) ) # take x's with corr > rxx with xi

		names_xx_input =  setdiff1(names_xx_thres, unlist(hppm))  # exclude x's previously chosen in hppm and set the rest as the input x's 

		names_xx_input_set = c( names_xx_input_set, names_xx_input )  # store the input x's

		if(verbose) cat( "j =",i, "\ninput :", names_xx_input, "\n" )

		if(length(names_xx_input) != 0  )  # if there is some input x
		{
			if(prior == "mom") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], family="binomial", priorCoef = momprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input[ which(bb $ postMode == 1) ] # collect the HPPM vars
			} else
			#
			#
			if(prior == "imom") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], family="binomial", priorCoef = imomprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input [ which(bb $ postMode == 1) ] # collect the HPPM vars
			} else
			#
			#
			if(prior == "zellner") 
			{
				bb = modelSelection( y, x = x[, names_xx_input, drop=F], family="binomial", priorCoef = zellnerprior(tau=tau), priorDelta = priorDelta, niter = niter, center=T, scale=T, verbose=F )  # NLP-MCMC with those vars only
				hppm[[i]] = names_xx_input [ which(bb $ postMode == 1) ] # collect the HPPM vars
			} 
			# emom not yet implemented in modelSelection. Also, emomLM won't run with only one x variable.
			if(verbose) cat( "selected :", hppm[[i]], "\n")  # print the HPPM vars
		}
	}

	names_xx_not_selected = setdiff1(names_xx_input_set, unlist(hppm) )

	return( list( hppm = unlist(hppm), not.selected = names_xx_not_selected ) )
}

setdiff1 <- function (a, b, no.dup.guaranteed = TRUE) {
  if (no.dup.guaranteed) {
    au <- a
    bu <- b
  } else {
    au <- unique(a)
    bu <- unique(b)
  }
  ind <- match(bu, au, nomatch = 0)
  DIFF <- au[-c(ind, length(au) + 1)]  ## https://stackoverflow.com/a/52772380
  DIFF
}
