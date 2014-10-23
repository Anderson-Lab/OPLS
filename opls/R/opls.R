library(pracma)
library(caret)

#' opls
#'
#' This allows you to create an OPLS model if you know the number of orthogonal components.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_OPLS_fact - Integer specifying the number of OPLS orthogonal components.
#' @return the found opls model
#' @export
#' @examples
#' model <- opls(rand(10,10),rand(10,1),1)
opls <-
  function(X,Y,num_OPLS_fact){
    w_ortho = matrix(nrow=ncol(X),ncol=num_OPLS_fact)
    t_ortho = matrix(nrow=nrow(X),ncol=num_OPLS_fact)
    p_ortho = matrix(nrow=ncol(X),ncol=num_OPLS_fact)
    
    Xres = scale(X,center=TRUE,scale=FALSE)
    Yres = scale(Y,center=TRUE,scale=FALSE)
    SS_Y=sum(sum(Yres^2))
    SS_X=sum(sum(Xres^2))
    
    if (num_OPLS_fact > 0) {
      for (iter in 1:num_OPLS_fact) {
        #find PLS component
        w = t(t(Yres)%*%Xres / (t(Yres)%*%Yres)[1])
        w = w / norm(w)
        t = Xres%*%w / (t(w)%*%w)[1]
        p = t(t(t)%*%Xres / (t(t)%*%t)[1])
        
        #run OSC filter on Xres
        w_ortho[,iter] = p - (t(w)%*%p / (t(w)%*%w)[1])[1] * w
        w_ortho[,iter] = w_ortho[,iter] / norm(matrix(w_ortho[,iter]))
        t_ortho[,iter] = Xres%*%w_ortho[,iter] / (t(w_ortho[,iter])%*%w_ortho[,iter])[1]
        p_ortho[,iter] = t(t(t_ortho[,iter])%*%Xres / (t(t_ortho[,iter])%*%t_ortho[,iter])[1])
        Xres = Xres - t_ortho[,iter]%*%t(p_ortho[,iter])
      }
    }
    
    
    ############
    # PLS on full data
    ############
    #find PLS component
    w = t(t(Yres)%*%Xres / (t(Yres)%*%Yres)[1])
    w = w / norm(w)
    t = Xres%*%w / (t(w)%*%w)[1]
    c = t(t(t)%*%Yres / (t(t)%*%t)[1])
    u = Yres%*%c / (t(c)%*%c)[1]
    p = t(t(t)%*%Xres / (t(t)%*%t)[1])
    # b coef
    b_l=((t(t)%*%t)[1]^(-1))*(t(u)%*%t)[1]
    
    #save model params
    model = list()
    model$b=b_l
    model$b_l=b_l
    model$c=c
    model$p=p
    model$w=w
    model$t=t
    model$u = u
    model$t_ortho = t_ortho
    model$p_ortho = p_ortho
    model$w_ortho = w_ortho
    model$num_OPLS_fact = num_OPLS_fact
    # Original space
    Wstar=w%*%solve(t(p)%*%w)
    B_pls=Wstar%*%b_l%*%t(c)
    Xres = scale(X,center=TRUE,scale=FALSE)
    z=Xres;
    # filter out OPLS components
    if (num_OPLS_fact > 0) {
      for (filter in 1:num_OPLS_fact) {
        z = (z - (z%*%w_ortho[,filter]/(t(w_ortho[,filter])%*%w_ortho[,filter])[1])%*%t(p_ortho[,filter]))
      }
    }
    #predict
    model$Y_pred = z%*%B_pls + mean(Y)
    
    model$R2_X=(t(model$t)%*%model$t)[1]*(t(model$p)%*%model$p)[1]/SS_X
    model$R2_Y=(t(model$t)%*%model$t)[1]*(model$b^2)[1]*(t(model$c)*model$c)[1]/SS_Y
    
    return(model)
  }

#' opls_CV
#'
#' This allows you to create an OPLS model if you know the number of orthogonal components.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_OPLS_fact - Integer specifying the number of OPLS orthogonal components.
#' @param folds - Number of k-fold cross-validation groups or -1 for leave one out cross-validation.
#' @return List containing \item{Q^2}{Cross-validated R^2} \item{Q2s}{One for each iteration} \item{press}{residual calculation used in Q^2 calculation} \item{Y_pred_validation}{Predicted values from the validation data} \item{accuracy}{standard accuracy that is only relevant if this is a classification problem}
#' @export
#' @examples
#' res <- opls_CV(rand(10,10),rand(10,1),1,-1)
opls_CV <-
  function(X,Y,num_OPLS_fact,folds) {
    press = 0
    y_sum = 0
    errors = 0
    Q2s = matrix(nrow=1,ncol=length(folds))
    Y_pred_validation = list()
    
    for (CV_count in 1:length(folds)) {
      # set up CV data
      Xtemp = X[-folds[[CV_count]],]
      Ytemp = Y[-folds[[CV_count]],]
      
      m = colMeans(Xtemp)
      m_Y = mean(Ytemp)
      Xres = sweep(Xtemp,2,m)
      Yres = Ytemp - m_Y    
      
      model = opls(Xres,Yres,num_OPLS_fact)
      
      # calc partial press
      X_leftOut = as.matrix(X[folds[[CV_count]],])
      if (length(folds[[CV_count]]) == 1) {
        X_leftOut = t(X_leftOut)
      }
      X_leftOut = sweep(X_leftOut,2,m)
      Y_leftOut = Y[folds[[CV_count]],]
      Y_leftOut = Y_leftOut -m_Y    
      temp_press = 0
      temp_y_sum = 0
      Y_preds = c()
      for (cpp in 1:length(Y_leftOut)) {
        Wstar=model$w%*%solve(t(model$p)%*%model$w)
        B_pls=Wstar%*%model$b_l%*%t(model$c)
        z=t(data.matrix(X_leftOut[cpp,]))
        # filter out OPLS components 
        if (num_OPLS_fact > 0) {
          for (filter in 1:num_OPLS_fact) {
            w_ortho = data.matrix(model$w_ortho[,filter])
            p_ortho = data.matrix(model$p_ortho[,filter])
            z = (z - (z%*%w_ortho/(t(w_ortho)%*%w_ortho)[1])%*%t(p_ortho))
          }
        }
        # predict
        Y_pred = z%*%B_pls
        Y_preds[cpp] = Y_pred + m_Y
        temp_press = temp_press + (Y_pred - Y_leftOut[cpp])^2
        temp_y_sum = temp_y_sum + (Y_leftOut[cpp])^2
        correct_Y = Y_pred - (Y_leftOut[cpp])
        for (k in 1:length(Y)) {
          if (abs(Y_pred - (Y[k])) < abs(correct_Y))
            errors = errors+1
          break
        }
      }
      Y_pred_validation[[CV_count]] = Y_preds
      Q2s[CV_count] = 1 - temp_press/temp_y_sum
      press = press + temp_press
      y_sum = y_sum + temp_y_sum
    }
    
    Q2 = 1 - press/y_sum
    accuracy = (length(Y)-errors) / length(Y)
    
    res = list(Q2=Q2,Q2s=Q2s,press=press,accuracy=accuracy,Y_pred_validation=Y_pred_validation)
    return(res)
  }

#' run_opls
#'
#' This will perform OPLS analysis including cross-validation and a permutation test.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_permutations - Number of permutation iterations.
#' @param CV - Number of k-fold cross-validation groups or -1 for leave one out cross-validation.
#' @param min_num_OPLS_fact - minimum number of OPLS orthogonal components. Default 0.
#' @return List containing \item{model}{Final OPLS model} \item{Q2}{Cross-validated R^2} \item{accuracy}{standard accuracy that is only relevant if this is a classification problem} \item{num_OPLS_fact}{number of orthogonal components} \item{permutation_Q2s}{Q2 from each permutation test} \item{pvalue}{p-value using the permutation test}
#' @export
#' @examples
#' res <- run_opls(rand(10,10),rand(10,1),100,-1)
run_opls <-
  function(X,Y,num_permutations,CV,min_num_OPLS_fact=0) {
    library(caret)
    
    if (CV == -1) {
      CV = length(Y) # Leave one out
    }
    folds = createFolds(Y,CV)
    
    num_OPLS_fact = min_num_OPLS_fact
    res = opls_CV(X,Y,num_OPLS_fact,folds)
    while (num_OPLS_fact < length(X[1,])) {
      next_res = opls_CV(X,Y,num_OPLS_fact+1,folds)
      R = next_res$press/res$press
      if (R >= 1)         
        break
      else {
        res = next_res
        num_OPLS_fact = num_OPLS_fact + 1
      }
    }
    
    model = opls(X,Y,num_OPLS_fact)
    cv_res = opls_CV(X,Y,num_OPLS_fact,folds)
    res = list(model=model,Q2=res$Q2,Q2s=res$Q2s,accuracy=res$accuracy,num_OPLS_fact=num_OPLS_fact)
    res$folds = folds
    res$Y_pred_validation = cv_res$Y_pred_validation
    
    res$permutation_Q2s = permutation_test(X,Y,num_OPLS_fact,num_permutations,folds)
    
    pvalue<-1 - length(which(res$Q2[1] >= res$permutation_Q2s))/length(res$permutation_Q2s)
    if (pvalue == 0) {
      res$pvalue = paste("<",1/length(res$permutation_Q2s))
    } else {
      res$pvalue = paste("<",pvalue)
    }
    
    return(res)
  }

#' determine_significant_features
#'
#' Permute the values of each variable and recompute its loadings. Compare this distribution of loadings with the original loading for each variable.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param orig_model - original model.
#' @param num_permutations - Number of permutation iterations.
#' @param talpha - test alpha to use as cutoff
#' @param inside_num_permutations - the smaller number of permutations used to efficiently rule out insignificant features. Should be less than num_permutations.
#' @param inside_talpha - corresponding test alpha for ruling out loadings that are not close to being significant. This is purely for efficient calculations.
#' @return List containing \item{sig_inxs}{Significant indices} \item{not_sig_inxs}{Indices that are not significant} \item{significant}{0/1 bit vector for significant variables} \item{p_permuted}{the permuted loading values}
#' @export
#' @examples
#' X <- rand(10,10)
#' Y <- rand(10,1)
#' model <- opls(X,Y,1)
#' res <- determine_significant_features(X,Y,model,500,0.05,100,0.2)
determine_significant_features <- function(X,Y,orig_model,num_permutations,talpha,inside_num_permutations,inside_talpha) {
  num_samples = nrow(X)
  num_variables = ncol(X)
  p_permuted = list()
  
  ## Determine the significant features for a model
  P_original = orig_model$p # Grab the original P
  # for each variable, permute the labels N times and recalculate P
  max_num_permutations = factorial(num_samples)
  #print('Maximum number of permutations is\n',max_num_permutations)
  # fprintf('The number of permutations is #d\n',N)
  # fprintf('The test alpha is #f\n',talpha)
  
  ## Inner permutation test
  variables = 1:num_variables
  N = min(c(inside_num_permutations,max_num_permutations))
  inner_P_permuted = run_det_sig(X,Y,orig_model,N,variables)
  
  # Determine the outright not significant bins
  significant = vector(length=num_variables)
  significant[1:length(significant)] = NaN
  sig_inxs = list()
  not_sig_inxs = list()
  for (v in 1:num_variables) {
    sorted = sort(inner_P_permuted[v,],decreasing=TRUE)
    inside_ix = max(c(1,round(N*inside_talpha/2))) # Two tailed
    inside_thres1 = sorted[inside_ix]
    inside_thres2 = sorted[length(sorted)-inside_ix+1]
    if (inside_thres1 >= P_original[v] && P_original[v] >= inside_thres2) { # Not close to the boundary
      p_permuted[[length(p_permuted)+1]] = inner_P_permuted[v,]
      not_sig_inxs[[length(not_sig_inxs)+1]] = v
      significant[v] = FALSE
    }
  }
  
  ## Border permutation test
  variables = which(is.nan(significant))
  N = min(c(num_permutations,max_num_permutations))
  border_P_permuted = run_det_sig(X,Y,orig_model,N,variables)
  
  # Determine the outright (not) significant bins
  sig_inxs = list()
  not_sig_inxs = list()
  for (v in variables) {
    sorted = sort(border_P_permuted[v,],decreasing=TRUE)
    ix = max(c(1,round(N*talpha/2))) # Two tailed
    thres1 = sorted[ix]
    thres2 = sorted[length(sorted)-ix+1]
    if (P_original[v] >= thres1) {
      p_permuted[[v]] = border_P_permuted[v,]
      significant[v] = TRUE
      sig_inxs[[length(sig_inxs)+1]] = v
    } else if (P_original[v] <= thres2) {
      p_permuted[[v]] = border_P_permuted[v,]
      significant[v] = TRUE
      sig_inxs[[length(sig_inxs)+1]] = v
    } else {
      p_permuted[[v]] = border_P_permuted[v,]
      not_sig_inxs[[length(not_sig_inxs)+1]] = v
      significant[v] = FALSE
    }
  }
  return(list(sig_inxs=unlist(sig_inxs),not_sig_inxs=unlist(not_sig_inxs),significant=significant,p_permuted=p_permuted))
}

#' run_det_sig
#'
#' This is an internal helper function. Do not call directly.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param orig_model - original model.
#' @param N - Number of permutation iterations.
#' @param variables - list of variables to test
#' @export
run_det_sig <- function(X,Y,orig_model,N,variables) {
  num_samples = nrow(X)
  num_variables = ncol(X)
  P_original = orig_model$p#/sqrt(sum(orig_model$p^2)) # Grab the original P
  
  num_OPLS_fact = orig_model$num_OPLS_fact
  P_permuted = matrix(NaN,nrow=N,ncol=num_variables)
  for (v in variables) {
    v_P_permuted = matrix(NaN,nrow=N,ncol=1)
    for (n in 1:N) {
#       #Newer version
#       remaining = 1:length(Y)
#       # Perform length(Y) swaps, but we need to make sure it is swapped with the other class
#       for (j in 1:(length(Y)-1)) {
#         while (T) {
#           six = round(1+(length(remaining)-1)*rand(1))
#           if (Y[six] == Y[j]) # Select another
#             next
#           X_permuted[j,v] = X[six,v]
#           #remaining = setdiff(remaining,six) # Remove the one that we swapped with
#           #remaining = setdiff(remaining,c(j+1)) # Remove the next one because we are getting ready to swap it
#           break
#         }
#       }
      X_permuted = X
      inxs = randperm(num_samples,num_samples)
      X_permuted[,v] = X[inxs,v]
      model = opls(X_permuted,Y,num_OPLS_fact)
      #model$p = model$p/sqrt(sum(orig_model$p^2)) # Grab the original P
      v_P_permuted[n] = model$p[v]
      # Make sure the direction of the vector is the same
      model$p[v] = 0 # Remove from calculation
      P_test = P_original
      #ixs = which(is.na(P_test))
      #P_test[ixs] = 0
      P_test[v] = 0      
      err1 = sum((P_test - model$p)^2)      
      err2 = sum(P_test - (-1*model$p)^2)
      if (err2 < err1) {
        v_P_permuted[n] = -1 * v_P_permuted[n]
      }
    }
    P_permuted[,v] = v_P_permuted
    #fprintf('Finished #d\n',v)
  }
  P_permuted = t(P_permuted)
  return (P_permuted)
}

#' apply_opls_model
#'
#' Apply a model to new (or old) data. Computes the t scores and the predicted-y value for each sampel in new_X. 
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param opls_model - opls model.
#' @param new_X - m x p matrix, where m is the number of samples and p is the number of variables.
#' @return List containing \item{t}{t-score} \item{t_ortho}{t-orthogonal scores} \item{Y_pred}{the predicted-y values for the samples}
#' @export
#' @examples
#' X <- rand(10,10)
#' new_X <- rand(5,10)
#' Y <- rand(10,1)
#' model <- opls(X,Y,1)
#' res <- apply_opls_model(X,Y,model,new_X)
apply_opls_model <- function(X,Y,opls_results,new_X) {
  model = opls_results$model
  t_ortho = matrix(nrow=nrow(new_X),ncol=opls_results$num_OPLS_fact)
  
  Xres = scale(new_X,center=colMeans(X),scale=FALSE)
  #Xres = bsxfun(@minus,new_X,mean(X));
  if (opls_results$num_OPLS_fact > 0) {
    
    for (iter in 1:opls_results$num_OPLS_fact) {
      #run OSC filter on Xres
      t_ortho[,iter] = Xres%*%model$w_ortho[,iter] / (t(model$w_ortho[,iter])%*%model$w_ortho[,iter])[1]
      Xres = Xres - t_ortho[,iter]%*%t(model$p_ortho[,iter])
    }
  }
  t = Xres%*%model$w / (t(model$w)%*%model$w)[1]
  
  # Original space
  Wstar=model$w%*%solve(t(model$p)%*%model$w)
  B_pls=Wstar%*%model$b_l%*%t(model$c)
  
  Xres = scale(new_X,center=colMeans(X),scale=FALSE)
  z=Xres
  # filter out OPLS components
  if (opls_results$num_OPLS_fact > 0) {
    for (filter in 1:opls_results$num_OPLS_fact) {
      z = (z - (z%*%model$w_ortho[,filter]/(t(model$w_ortho[,filter])%*%model$w_ortho[,filter])[1])%*%t(model$p_ortho[,filter]))
    }
  }
  #predict
  Y_pred = z%*%B_pls + mean(Y)
  
  return(list(t=t,t_ortho=t_ortho,Y_pred=Y_pred))
}

#' n.group.opls.helper
#'
#' An internal helper function. Do not call this directly.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_permutations - number of permutation for the randomization test.
#' @param CV - number of folds for k-fold cross-validation or -1 for leave one out.
#' @param min_num_OPLS_fact - minimum number of OPLS factors to consider.
#' @export
n.group.opls.helper = function(X,Y,num_permutations,CV,min_num_OPLS_fact=0) {
  uniqueY = unique(Y)
  pairs = list()
  Q2s = c()
  opls.history = list()
  c = 1
  for (i1 in 1:length(uniqueY)) {
    if (i1+1 > length(uniqueY))
      break
    for (i2 in (i1+1):length(uniqueY)) {
      class1 = uniqueY[i1]
      class2 = uniqueY[i2]
      ixs = which(Y == class1 | Y == class2)
      opls<-run_opls(X[ixs,],as.matrix(Y[ixs,]),0,CV,min_num_OPLS_fact=min_num_OPLS_fact) #runs OPLS-DA
      pairs[[c]] = c(i1,i2)
      Q2s[c] = opls$Q2
      c = c + 1
    }
  }
  mx = max(Q2s)
  classes = uniqueY[c(pairs[[which(Q2s==mx)]])]
  classesIx = c(pairs[[which(Q2s==mx)]])
  remainingY = uniqueY[-classesIx]
  modelIxs = which(Y == classes[1] | Y == classes[2])
  adjustedY = as.matrix(Y[modelIxs,])
  ixs = which(adjustedY == classes[1])
  adjustedY[ixs] = 1
  ixs = which(adjustedY == classes[2])
  adjustedY[ixs] = 2
  origUniqueY = classes
  newUniqueY = c(1,2)
  opls<-run_opls(X[modelIxs,],adjustedY,num_permutations,CV,min_num_OPLS_fact=min_num_OPLS_fact) #runs OPLS-DA
  opls.history[[length(opls.history)+1]] = opls
  
  # We now have our first pair, so the next goals is to continue to add classes to the
  # list. One at a time. Each time picking the one that gives you highest Q2. We will
  # want to report how the Q2 degrades over the course of adding the classes.
  while (length(remainingY) > 0) {
    Q2s = c()
    for (i1 in 1:length(remainingY)) {
      class1 = remainingY[i1]
      ixs = which(Y == class1)
      opls_results_new<-apply_opls_model(X[modelIxs,],adjustedY,opls,X[ixs,]) #projects new groups onto first model
      newY = median(opls_results_new$Y_pred)
      opls<-run_opls(X[c(modelIxs,ixs),],as.matrix(c(adjustedY[,1],rep(newY,length(ixs)))),num_permutations,CV,min_num_OPLS_fact=min_num_OPLS_fact)
      Q2s[i1] = opls$Q2
    }
    mx = max(Q2s)
    c = which(Q2s == mx)
    class1 = remainingY[c]
    ixs = which(Y == class1)
    origUniqueY[length(origUniqueY)+1] = class1
    opls_results_new<-apply_opls_model(X[modelIxs,],adjustedY,opls,X[ixs,]) #projects new groups onto first model
    newY = median(opls_results_new$Y_pred)
    newUniqueY[length(newUniqueY)+1] = newY
    opls<-run_opls(X[c(modelIxs,ixs),],as.matrix(c(adjustedY[,1],rep(newY,length(ixs)))),num_permutations,CV,min_num_OPLS_fact=min_num_OPLS_fact)
    opls.history[[length(opls.history)+1]] = opls
    adjustedY = as.matrix(c(adjustedY[,1],rep(newY,length(ixs))))
    modelIxs = c(modelIxs,ixs)
    remainingY = remainingY[-c]
  }
  return(list(opls=opls,opls.history=opls.history,origUniqueY=origUniqueY,newUniqueY=newUniqueY,adjustedY=adjustedY))
}

#' n.group.opls 
#'
#' Create an OPLS-DA model when there are more than 2 classes.  This function is not applicable if you have a continuous response variable. It iteratively adds classes to the model. New classes are added such that they maximize the Q2. The class label is determine by the previous model (i.e., new data is projected into the previous model).
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_permutations - number of permutation for the randomization test.
#' @param CV - Parameter for internal cross-validation. -1 for leave-one-out cross-validation. The value of k in k-fold cross-validation otherwise.
#' @param nIterations - number of iterations for external validation. One of each sample is held out each iteration.
#' @param min_num_OPLS_fact - minimum number of OPLS factors. Default 0.
#' @return List containing \item{Q2}{External Q2 value} \item{helper.results}{Results from running helper function, including an opls model, opls model history, original unique Y values, new unique Y values, and adjusted Y values}
#' @export
#' @examples
#' X <- rand(12,10)
#' new_X <- rand(5,10)
#' Y <- rand(12,1)
#' Y[1:4,] = 1
#' Y[5:8,] = 2 
#' Y[9:12,] = 3 
#' res <- n.group.opls(X,Y,100,-1)
n.group.opls = function(X,Y,num_permutations,CV,nIterations=100,min_num_OPLS_fact=0) {
  uniqueY = unique(Y)
  # We don't need to determine the group labels
  if (length(uniqueY) == 2) {
    results = run_opls(X,Y,num_permutations,CV,min_num_OPLS_fact=min_num_OPLS_fact)
    return(list(Q2=results$Q2,helper.results=results))    
  }
  
  squaredResiduals = c()
  squaredResidualsMean = c()
  q = 1
  for (iter in 1:nIterations) {
    leaveOutIxs = c()
    c = 1
    for (cls in uniqueY) {
      ixs = which(Y == cls)
      leaveOutIxs[c] = sample(ixs,1)
      c = c + 1
    }
    leaveOutX = X[leaveOutIxs,]
    trainingX = X[-leaveOutIxs,]
    leaveOutY = as.matrix(Y[leaveOutIxs,])
    trainingY = as.matrix(Y[-leaveOutIxs,])
    results = n.group.opls.helper(trainingX,trainingY,0,CV,min_num_OPLS_fact=min_num_OPLS_fact)
    
    opls_results_holdOut<-apply_opls_model(trainingX,results$adjustedY,results$opls,leaveOutX)
    newLeaveOutY = c()
    for (i in 1:length(leaveOutY)) {
      ix = which(results$origUniqueY == leaveOutY[i])
      newLeaveOutY[i] = results$newUniqueY[ix]
      squaredResiduals[q] = (opls_results_holdOut$Y_pred[i]-newLeaveOutY[i])^2
      squaredResidualsMean[q] = (mean(results$adjustedY)-newLeaveOutY[i])^2
      q = q + 1
    }
  }
  Q2 = 1 - sum(squaredResiduals)/sum(squaredResidualsMean)      
  
  # Now create a model with the full data
  results = n.group.opls.helper(X,Y,num_permutations,CV,min_num_OPLS_fact=min_num_OPLS_fact)
    
  #pvalue<-1 - length(which(results$opls$Q2[1] >= results$opls$permutation_Q2s))/length(results$opls$permutation_Q2s)
  
  return(list(Q2=Q2,helper.results=results))
}

#' permutation_test
#'
#' Internal helper function. Do not use.
#' @param X - n x p matrix, where n is the number of samples and p is the number of variables.
#' @param Y - n x 1 matrix. Must be numeric 
#' @param num_OPLS_fact - number of orthogonal OPLS factors
#' @param num_permutations - number of permutation for the randomization test.
#' @param folds - number of folds for CV
#' @export
permutation_test <-
  function(X,Y,num_OPLS_fact,num_permutations,folds) {
    #permutation
    permutation_Q2s = vector(length=num_permutations)
    for (i in 1:num_permutations) {
      #permute labels
      Y_rand_idx = randperm(length(Y))
      Y_rand = matrix(Y[Y_rand_idx])
      
      #run OPLS on permuted data
      res = opls_CV(X,Y_rand,num_OPLS_fact,folds)
      permutation_Q2s[i] = res$Q2
    }
    return(permutation_Q2s)
  }


# Example usage
if (F) {
  Y<-as.matrix(as.numeric(as.factor(NMRData.m$Species)))
  labels = as.factor(NMRData.m$Species)
  uLabels = unique(labels)
  results<-n.group.opls(X,Y, 100, -1, 100)
  uYs = unique(as.character(results$helper.results$adjustedY))
  cols = as.character(results$helper.results$adjustedY)
  cols[cols==uYs[1]] = 'royalblue1'
  cols[cols==uYs[2]] = 'darkgray'
  cols[cols==uYs[3]] = 'mediumseagreen'
  cols[cols==uYs[4]] = 'goldenrod'
  cols[cols==uYs[5]] = 'darkred'
  model = results$helper.results$opls$model
  plot(model$t, model$t_ortho[,1], col=cols, pch=(cols), cex=2, xlab="t",ylab="t-orthogonal")
  legend('topleft', legend=uLabels, col=c('royalblue1','darkgray','mediumseagreen','goldenrod','darkred'), pch=c(6,5), cex=.5, bty='n')
  
  num_permutations <- 500; talpha <- 0.01; inside_num_permutations = 100; inside_talpha = 0.2
  # Now if you look at results$helper.results$adjustedY, there are really two groups that are the same, 1.39 and 1.354
  adjustedYmanually = round(results$helper.results$adjustedY*10)/10
  adjustedOPLSresults = run_opls(X,adjustedYmanually,num_permutations,CV,min_num_OPLS_fact=0)
  sig_features_results <- determine_significant_features(X,adjustedYmanually,results$helper.results$opls$model,num_permutations,
                                                         talpha,inside_num_permutations,inside_talpha)  
  sig_features_results$sig_inxs
  
  Y=results$helper.results$adjustedY
  orig_model = results$helper.results$opls$model
  N = 100
  variables = 82
  test.results = run_det_sig(X,Y,orig_model,1000,variables)
}
