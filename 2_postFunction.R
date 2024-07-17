# posterior mean of parameters -----------

post.mean <- function(df_of_draws, dat, nBasis){
  nSub <- length(unique(dat$ID))
  b0 <- mean(df_of_draws[,'beta0'])
  bMatrix <- matrix(
    apply(df_of_draws[,sprintf('beta[%d,%d]',rep(seq(nSub), each=nBasis),seq(nBasis))],
          2, mean),
    nrow=nSub,byrow = T)
  G <-  matrix(
    apply(df_of_draws[,sprintf('VarBeta[%d,%d]',rep(seq(nBasis), each=nBasis),seq(nBasis))],
          2, mean),
    nrow=nBasis,byrow = T)
  Gamma <-
    mean(df_of_draws[,'sigmaepsilon'])^2

  list(b0=b0,bMatrix=bMatrix,G=G,Gamma=Gamma)
}

# fitted values based on posterior mean ----------------------

pc.fit <- function(ls_par, dat, ls_idxA,seed){
  mat_fitted <- data.frame(
    matrix(NA,nrow=nrow(dat),ncol=length(ls_idxA),dimnames=list(c(),sprintf('y_hatA%d',seq(length(ls_idxA))))))
  id <-  as.numeric(dat$ID)

  b0 <- ls_par$b0
  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma


  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    set.seed(seed)

    cat(sprintf('Set %d: calculating fitted curve at MAP\n',t))

    bAMatrix <- bMatrix[,idxA]

    mat_fitted_t <- my.fit(idxA, idxB, dat, G,bAMatrix,b0, id)
    list(mat_fitted_t)
  }->temp
  list(
    mat_fitted=matrix(unlist(lapply(temp,'[[',1)),ncol=length(ls_idxA),byrow=F,dimnames=list(c(),sprintf('y_hatA%d',seq(length(ls_idxA)))))
  )
}

# cluster number by KL method --------------------------------

pc.KL <- function(ls_par, dat, ls_idxA, nIter, thKL, regQ, seed){
  id <-  as.numeric(dat$ID)
  nSub <- length(unique(id))

  bMatrix <- ls_par$bMatrix
  G <- ls_par$G
  Gamma <- ls_par$Gamma


  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    set.seed(seed)

    cat(sprintf('Set %d: optimizing cluster number by k means\n',t))

    bAMatrix <- bMatrix[,idxA]
    Q <- Q.func(idxA, idxB, dat,Gamma, G,id)
    Kcluster.obj <- Kcluster.pick(nIter, bAMatrix,Q, regQ, do.plot=F, seed=seed)
    threshold <- Kcluster.obj$KL[1]*thKL # find ncluster such that KL<threshold
    nClusters_t <- which(Kcluster.obj$KL<threshold)[1]
    KLs_t <- Kcluster.obj$KL
    cluster_t <- Kcluster.obj$clusterMatrix

    list(nClusters_t, KLs_t, cluster_t)
  }->temp

  list(
    nClusters=unlist(lapply(temp,'[[',1)),
    KLs=lapply(temp,'[[',2),
    cluster0=lapply(temp,'[[',3)
  )
}

# draw from posterior with chosen nCluster --------------------

pc.pair <- function(df_of_draws, dat, nBasis, ls_idxA, nIter, nDraw,nClusters, regQ, seed){
  id <-  as.numeric(dat$ID)
  nSub <- length(unique(id))
  if(length(nClusters)==1) nClusters <- rep(nClusters,length(ls_idxA))

  foreach(t = seq(1,length(ls_idxA))) %dopar% {
    idxA <- ls_idxA[[t]]
    idxB <- setdiff(ls_idxA[[1]],idxA)
    nCluster <- nClusters[t]

    set.seed(seed)

    cat(sprintf('Set %d: drawing samples from posterior\n',t))
    post.obj <-
      cluster.posterior(df_of_draws,nDraw, dat, nBasis,idxA, idxB, id, nCluster,nIter, regQ, seed)

    # # posterior of KL
    # hist(post.obj$KL,main=sprintf('Posterior of KL with %d clusters',nCluster))
    #
    # # posterior of centroid
    # plot(post.obj$dKArray[,1,],post.obj$dKArray[,2,],pch='.', main=sprintf('Posterior of centroids with %d clusters',nCluster))

    # posterior probability of 2 subjects in 1 cluster
    tb <- matrix(NA, nrow= nSub, ncol=nSub)
    for(i in seq(2,nSub)){
      for(j in seq(i-1)){
        tb[i,j] <- mean(post.obj$clusterMatrix[i,]==post.obj$clusterMatrix[j,])
      }
    }
    ls_prob_t <- tb
    list(ls_prob_t,
         post.obj$clusterMatrix, post.obj$clusterVec)
  }->temp

  list(
    ls_prob=lapply(temp,'[[',1),
    ls_clust=lapply(temp,'[[',2),
    ls_clust0=lapply(temp,'[[',3)
  )
}


# cluster number by bootstrap on the fitted curve -----------------

nCluster.boot <- function(mat_fitted, dat, nB, nCl, r=5, threshold=0.5, scheme_2 = TRUE, seed){
  set.seed(seed)
  nSet <- ncol(mat_fitted)
  sBtb <- matrix(NA, nrow=nCl-1, ncol=nSet)
  nClusters <- numeric(nSet)
  
  for(j in seq(nSet)){
    xx <- dat %>% # formatting observations for distance
      select(ID,t) %>%
      bind_cols(data.frame(y=mat_fitted[,j])) %>%
      pivot_wider(names_from = t, values_from=y) %>%
      select(-ID) %>%
      as.matrix
    nSub <- nrow(xx)
    
    sB <- sapply(seq(2,nCl), function(cl){
      print(cl)
      valid_result <- FALSE
      while(!valid_result){
        try({
          min.agr <- numeric(length=nB)
          clst.mat <- matrix(NA_real_, nrow=nB+1, ncol=nSub)
          
          km <- kmeansCBI(data=xx, k=cl, scaling=FALSE, runs=r)
          center <- km$result$centers
          clst <- km$result$cluster
          clst.mat[1,] <- clst
          
          stab.matrix <- matrix(NA_real_, nrow=nSub, ncol=nB)
          
          for(b in 1:nB){
            ind <- sample(nSub,nSub,replace=T)
            xx.star <- xx[ind,]
            km.star <- kmeansCBI(data=xx.star, k=cl, scaling=FALSE, runs=r)
            center.star <- km.star$result$centers
            class.star <- mapping.Euclidean(center.star, xx, 1:cl)
            clst.mat[(b+1),] <- class.star
          }
          
          nB1 <- nB+1
          agree.mat <- matrix(NA, nrow=nB1, ncol=nB1)
          diag(agree.mat) <- 1
          
          for(i in 1:(nB1-1)){
            for(j in (i+1):nB1){
              agree.mat[i,j] <- mean(agreement(clst.mat[i,], clst.mat[j,]))
              agree.mat[j,i] <- agree.mat[i,j]
            }
          }
          
          mean.agr <- rowMeans(agree.mat)
          ref <- which.max(mean.agr)
          clst.ref <- clst.mat[ref, , drop = FALSE]
          
          stab.mat <- matrix(NA, nrow = nB1, ncol = nSub)
          
          for(i in 1:nB1){
            stab.mat[i,] <- agreement(clst.ref, clst.mat[i,])
          }
          
          if(!scheme_2){ ref <- 1 }
          
          ref.scheme1 <- clst.mat[1,] #
          c.mat.scheme1 <- clst.mat[-1,] #
          min.agr.scheme1 <- c()
          
          ref.cl <- clst.mat[ref,]
          c.mat.1 <- clst.mat[-ref,]
          for(i in 1:nB){
            min.agr[i] <- min.agreement(ref.cl, agreement(ref.cl, c.mat.1[i,]))
          }
          temp <- mean(min.agr)
          valid_result <- !is.na(temp)
        }, silent = TRUE)
      }
      temp
    })
    sBtb[,j] <- sB
    if(any(sB > threshold)){
      nClusters[j] <- max(which(sB > threshold))+1
    }else{
      nClusters[j] <- 1
    }
  }
  ls_sB <- list(sBtb=sBtb, nClusters=nClusters)
  ls_sB
}


jaccard <- function(set1, set2){
  jaccard <- length(intersect(set1, set2))/length(union(set1, set2))
  return(jaccard)
}


mapping.Euclidean <- function(center.ori, center.map, label.ori){
  dist.mat <- dist2(center.map, center.ori)
  mapping.pre <- apply(dist.mat, 1, which.min)
  mapping <- label.ori[mapping.pre]
  return(mapping)
}


agreement <- function(clst1, clst2){
  n1 <- length(clst1)
  n1_occur <- data.frame(table(clst1))
  nk1 <- as.numeric(n1_occur[n1_occur$Freq>1,]$clst1)
  
  n2 <- length(clst2)
  n2_occur <- data.frame(table(clst2))
  nk2 <- as.numeric(n2_occur[n2_occur$Freq>1,]$clst2)
  if(n1!=n2)warning('sample size is not equal')
  n <- n1
  if(length(nk1)==0|length(nk2)==0){
    stab.vec <- rep(0,n)
  }else{
    cluster.sets.1 <- list()
    cluster.sets.2 <- list()
    for(i in nk1){
      cluster.sets.1[[i]] <- which(clst1==i)
    }
    
    for(i in nk2){
      cluster.sets.2[[i]] <- which(clst2==i)
    }
    
    jaccard.matrix <- matrix(0, nrow=length(unique(clst1)), ncol=length(unique(clst2)))
    for(i in 1:length(nk1)){
      for(j in 1:length(nk2)){
        jaccard.matrix[i,j] <- jaccard(cluster.sets.1[[i]], cluster.sets.2[[j]])
      }
    }
    
    stab.vec <- c()
    for(i in 1:n){
      memb <- which(levels(n1_occur$clst1)==clst1[i])
      memb.star <- which(levels(n2_occur$clst2)==clst2[i])
      stab.vec[i] <- jaccard.matrix[memb, memb.star]
    }
    
  }
  return(stab.vec)
}


min.agreement <- function(clst, agrmt){
  clst.list <- unique(clst)
  clst.sta <- c()
  n <- length(clst.list)
  for(i in 1:n) {
    clst.sta[i] <- mean(agrmt[clst==clst.list[i]])
  }
  min.agrmt <- min(clst.sta)
  return(min.agrmt)
}
