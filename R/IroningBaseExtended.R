### SMOOTHING AND IRONING ###################
#this code smooths existing likelihood profiles, dealing with both spiked values and steep cliffs
testtestest <- function(){
####STUFF YOU NEED TO SPECIFY####

#Fit the model to the experimental data
all.names <- c("sigmaN.CM",
               "sigmaN.EMS",
               "rhoTCM",
               "rhoTCMbase",
               "sigmaCM.EMS",
               "sigmaCM.EMSbase",
               "sigmaEMS.CM",
               "sigmaEMS.CMbase",
               "rhoTEMS",
               "rhoTEMSbase",
               "sigmaEMS.EML",
               "sigmaEMS.EMLbase",
               "sigmaEML.EMS",
               "sigmaEML.EMSbase",
               "rhoTEML",
               "rhoTEMLbase",
               "sigmaEML.RM",
               "sigmaEML.RMbase",
               "sigmaRM.EML",
               "sigmaRM.EMLbase",
               "rhoRM",
               "rhoRMbase",
               "sigmaEMS.RM",
               "sigmaEMS.RMbase",
               "sigmaRM.EMS",
               "sigmaRM.EMSbase",
               "Length2",
               "Length3")

#data variables
data <- rbind(c(1, 'sigmaN.CM', c(0.1,10,0.1)),
              c(2, 'rhoTCMbase', c(0,0.4,0.005)),
              c(3, 'sigmaCM.EMS', c(0,10,0.1)),
              c(4, 'sigmaCM.EMSbase', c(0,0.3,0.005)),
              c(5, 'sigmaEMS.CM', c(0,10,0.1)),
              c(6, 'rhoTEMS', c(1,5,0.1)),
              c(7, 'sigmaEMS.EMLbase', c(0.1,0.5,0.01)),
              c(8, 'rhoTEML', c(2,5,0.1)),
              c(9, 'rhoTEMLbase', c(-4,0, 0.1)),
              c(10, 'sigmaEML.RM', c(0,0.001, 0.00001)),
              c(11, 'rhoRM', c(0,3, 0.1)),
              c(12, 'Length2', c(0.1,4,0.05)),
              c(13, 'Length3', c(0.0,3,0.05))
)

#define initial values####
init.cond <- c(N = 100 , TCM =  0, TEMS = 0, TEML = 0, RM = 0, Ag = 1)


#folder to save files in
main.folder <- "PLAMSRestricted"
folder <- paste0("/net/home.isilon/ag-graw/Michael/Malaria/",main.folder,"/PLof",data[rowx,2],".rds")

##if this function is not called via AutoSubmission specify
#rowx = the parameter index of the profile
#name.par = the name of the respective parameter


#specify at which difference in the MLE values the program should try to refit it
#this takes care of steep drops in your profile
trigger <- 0.3

#specify the minimal difference between values that will be taken into account when searching for spikes
no.diff <- 0.01

#specify the number of limit optimisation
#after testing the left and right parameter combinations for spiked values, the program continues by taking random values from intervals based upon the left and right parameters
lim.number <- 1

#afterwards the program samples from the intervals based on the highest on lowest value of the WHOLE data set. Specify the number of tries
###NOTE: Does not work very well!
all.number <- 0


#read in previously generated data
#data is of form c(fixed parameter (value), corresponding MLE, remaining fitted parameters)

if(file.exists(paste0("/net/home.isilon/ag-graw/Michael/Malaria/",main.folder,"/PLof",data[rowx,2],".rds"))){
  print("load previously smoothed data")
  pl <- readRDS(paste0("/net/home.isilon/ag-graw/Michael/Malaria/",main.folder,"/PLof",data[rowx,2],".rds"))
}else{
  pl <- as.data.frame(read.table(paste0("/net/home.isilon/ag-graw/Michael/Malaria/",main.folder,"/PLof",
                                        data[rowx,2],".csv"),header = F))
  #name the data frame
  colnames(pl) <- c(name.par, "MLE", data[-rowx,2])
}


#specify columns of fitted parameters in data frame
par.range <- c(3,ncol(pl))


####LIBRARIES AND FUNCTIONS#####

#specify range of the fixed parameter
range <- pl[,1]
#index of the parameter
pl.par = rowx

#print first two columns for comparisons in .Rout
print(pl[,1:2])

#redefine data frame
prs.two <- pl

#set count variables for first run
improvement <- 1
badfit <- 1

#run until refitting of cliffs does not yield any improvement and no more spikes are present
###MAIN LOOP####
while(improvement > 0 || length(badfit) > 0){
  ###SMOOTHING####
  #find steep cliffs and specify the type: 1 - cliff to the left, 2- cliff to the right, 12 - spiked value
  sel.steepcliff <- rep(0, length(range))
  #cycle through generated profile
  for(i in 1:length(range)){
    #check left border
    if(i==1){
      if((prs.two[i,2] - prs.two[i+1,2]) > trigger){
        sel.steepcliff[i] <- 2
      }
    }else if(i==length(range)){#check right border
      if((prs.two[i,2] - prs.two[i-1,2]) > trigger){
        sel.steepcliff[i] <- 1
      }
    }else{#check all values inbetween
      if((prs.two[i,2] - max(prs.two[c(i-1,i+1),2]))> trigger){#if value is a cliff peak
        sel.steepcliff[i] <- 12
      }else if ((prs.two[i,2] - prs.two[i-1,2]) > trigger){#if cliff is to the left
        sel.steepcliff[i] <- 1
      }else if ((prs.two[i,2] - prs.two[i+1,2]) > trigger){#if cliff is to the right
        sel.steepcliff[i] <- 2
      }

    }

  }
  #get problematic entries
  steepcliff <- which(sel.steepcliff > 0)
  print(paste0("Cliff higher than ", trigger))
  print(steepcliff)

  #define improvement variable
  improvement <- 0

  #fit only if cliffs are present
  if(length(steepcliff)>0){
    #cycle through cliffs
    for(i in 1:length(steepcliff)){
      k <- steepcliff[i]
      cliff <- sel.steepcliff[k]

      #set seed for random fitting
      set.seed(unclass(Sys.time()))

      #use the parameter values of the better side of the cliff
      if(cliff == 1){
        params.test <- as.numeric(prs.two[k - 1,par.range[1]:par.range[2]])
        names(params.test) <- names(prs.two[k - 1,par.range[1]:par.range[2]])
      }else if(cliff == 2){
        params.test <- as.numeric(prs.two[k + 1,par.range[1]:par.range[2]])
        names(params.test) <- names(prs.two[k + 1,par.range[1]:par.range[2]])
      }else{#cliff == 12
        params.test <- as.numeric(prs.two[k - 1,par.range[1]:par.range[2]])
        names(params.test) <- names(prs.two[k - 1,par.range[1]:par.range[2]])

        params.test2 <- as.numeric(prs.two[k + 1,par.range[1]:par.range[2]])
        names(params.test2) <- names(prs.two[k + 1,par.range[1]:par.range[2]])
      }

      print(paste0("Index ",k))
      print(paste0("Use side ",cliff))
      print(paste0("Best value is currently ", prs.two[k,2]))
      #try and see if parameter combination work
      works <- is.finite(ls2(par = params.test,
                             init = init.cond,
                             vals = range[k],
                             namez = name.par))
      #if parameter combination failed, try 100 times to create a working one by jittering
      if(works == FALSE){
        #save old parameter set
        params.base <- params.test
        #define count variable
        count.test <- 0
        while(works == FALSE & count.test < 100){
          #create new parameter set for testing
          params.test <- runif(n = length(params.base),
                               min = params.base - 0.05*sign(params.base)*params.base,
                               max = params.base + 0.05*sign(params.base)*params.base)
          names(params.test) <- names(params.base)
          #test the set
          works <-  is.finite(ls2(par = params.test,
                                  init = init.cond,
                                  vals = range[k],
                                  namez = name.par))
          #update count variable
          count.test <- count.test + 1
        }
      }
      #run optimization only if parameter combination works
      if(works == TRUE){

        #optimise
        p.scale <- par.scales(par = abs(params.test),
                              scale = abs(params.test),
                              fix = length(params.test))

        opt <- optms(params = params.test, range = range[k])

        #jump to next value if fitting succeded,
        if(opt$value < prs.two[k,2]){
          #overwrite data frame
          prs.two[k,] <- c(range[k], opt$value, opt$par)
          #update improvement
          improvement <- improvement + 1
          #save results
          saveRDS(prs.two, folder)
          next
        }

      }
      #if the value is spiked, try parameter combination of the right side
      if(cliff == 12){

        #try and see if parameter combination work
        works <- is.finite(ls2(par = params.test2,
                               init = init.cond,
                               vals = range[k],
                               namez = name.par))
        #if parameter combination failed, try 100 times to create a working one by jittering
        if(works == FALSE){
          #save old parameter set
          params.base <- params.test2
          #define count variable
          count.test <- 0
          while(works == FALSE & count.test < 100){
            #create new parameter set for testing
            params.test2 <- runif(n = length(params.base),
                                  min = params.base - 0.05*sign(params.base)*params.base,
                                  max = params.base + 0.05*sign(params.base)*params.base)
            names(params.test2) <- names(params.base)
            #test the set
            works <-  is.finite(ls2(par = params.test2,
                                    init = init.cond,
                                    vals = range[k],
                                    namez = name.par))
            #update count variable
            count.test <- count.test + 1
          }
        }
        #run optimization only if parameter combination works
        if(works == TRUE){

          #optimise
          p.scale <- par.scales(par = abs(params.test2),
                                scale = abs(params.test2),
                                fix = length(params.test2))

          opt <- optms(params = params.test2, range = range[k])

          #jump to next value if fitting succeded,
          if(opt$value < prs.two[k,2]){
            #overwrite data frame
            prs.two[k,] <- c(range[k], opt$value, opt$par)
            #update improvement
            improvement <- improvement + 1
            #save results
            saveRDS(prs.two, folder)
            next
          }

        }
      }


    }
  }
  print(paste0("Improved ", improvement, " out of ", length(steepcliff), " cliffs"))

  #IRONING####
  #find spiked values
  badfit <- rep(0, length(range))
  for(i in 2:(length(range)-1)){
    if((prs.two[i,2] > (prs.two[i-1,2] + no.diff) && #check if higher than left value
        prs.two[i,2] > (prs.two[i+1,2] + no.diff)) #|| #check if higher than right value
    ){
      badfit[i] <- 1
    }
  }
  #display spiked values
  badfit <- which(badfit==1)
  print(paste0("Spiked values are found in:"))
  print(badfit)

  #define a limit when ironing should stop to try smoothing for cliffs again
  limit1 <- 0
  #repeat until there are no spiked values or limit is reached
  #while(length(badfit) >= 1 && limit1 < 5){
  if(length(badfit) >= 1){
  #cycle through bad fits
  for(i in 1:length(badfit)){
    k <- badfit[i]
    print(paste0("Has to get better than ",max(prs.two[c(k-1, k+1),2]) + no.diff))

    #try the smaller value first
    if(prs.two[k - 1,2] < prs.two[k + 1,2]){
      first <- k - 1
      second <- k + 1
      print("Try left, then right")
    }else{
      first <- k + 1
      second <- k - 1
      print("Try right, then left")
    }

    #set seed for random fitting
    set.seed(unclass(Sys.time()))

    #first run####
    #get parameter combination from better side
    params.first <- as.numeric(prs.two[first,par.range[1]:par.range[2]])
    names(params.first) <- names(prs.two[first,par.range[1]:par.range[2]])

    #try and see if parameter combination work
    #try and see if parameter combination work
    works <- is.finite(ls2(par = params.first,
                           init = init.cond,
                           vals = range[k],
                           namez = name.par))
    #if parameter combination failed, try 100 times to create a working one by jittering
    if(works == FALSE){
      #save old parameter set
      params.base <- params.first
      #define count variable
      count.test <- 0
      while(works == FALSE & count.test < 100){
        #create new parameter set for testing
        params.first <- runif(n = length(params.base),
                              min = params.base - 0.05*sign(params.base)*params.base,
                              max = params.base + 0.05*sign(params.base)*params.base)
        names(params.first) <- names(params.base)
        #test the set
        works <-  is.finite(ls2(par = params.first,
                                init = init.cond,
                                vals = range[k],
                                namez = name.par))
        #update count variable
        count.test <- count.test + 1
      }
    }
    #run optimization only if parameter combination works
    if(works == TRUE){

      #optimise
      p.scale <- par.scales(par = abs(params.first),
                            scale = abs(params.first),
                            fix = length(params.first))

      opt <- optms(params = params.first, range = range[k])

      #jump to next value if fitting succeded,
      if(opt$value < (max(prs.two[c(k-1, k+1),2]) + no.diff)){
        prs.two[k,] <- c(range[k], opt$value, opt$par)
        saveRDS(prs.two, folder)
        next
      }

    }
    #second run####
    #use worse parameter combination
    print("second")
    params.second <- as.numeric(prs.two[second,par.range[1]:par.range[2]])
    names(params.second) <- names(prs.two[second,par.range[1]:par.range[2]])
    #optimise

    #try and see if parameter combination work
    works <- is.finite(ls2(par = params.second,
                           init = init.cond,
                           vals = range[k],
                           namez = name.par))
    #if parameter combination failed, try 100 times to create a working one by jittering
    if(works == FALSE){
      #save old parameter set
      params.base <- params.second
      #define count variable
      count.test <- 0
      while(works == FALSE & count.test < 100){
        #create new parameter set for testing
        params.second <- runif(n = length(params.base),
                               min = params.base - 0.05*sign(params.base)*params.base,
                               max = params.base + 0.05*sign(params.base)*params.base)
        names(params.second) <- names(params.base)
        #test the set
        works <-  is.finite(ls2(par = params.second,
                                init = init.cond,
                                vals = range[k],
                                namez = name.par))
        #update count variable
        count.test <- count.test + 1
      }
    }
    #run optimization only if parameter combination works
    if(works == TRUE){

      opt <- optms(params = params.second, range = range[k])

      #jump to next value if fitting succeded,
      if(opt$value < (max(prs.two[c(k-1, k+1),2]) + no.diff)){
        prs.two[k,] <- c(range[k], opt$value, opt$par)
        saveRDS(prs.two, folder)
        next
      }

    }

    #Same value####
    #use worse parameter combination
    print("same")
    params.same <- as.numeric(prs.two[k,par.range[1]:par.range[2]])
    names(params.same) <- names(prs.two[k,par.range[1]:par.range[2]])
    #optimise

    #try and see if parameter combination works
    #try and see if parameter combination work
    works <- is.finite(ls2(par = params.same,
                           init = init.cond,
                           vals = range[k],
                           namez = name.par))
    #if parameter combination failed, try 100 times to create a working one by jittering
    if(works == FALSE){
      #save old parameter set
      params.base <- params.same
      #define count variable
      count.test <- 0
      while(works == FALSE & count.test < 100){
        #create new parameter set for testing
        params.same <- runif(n = length(params.base),
                             min = params.base - 0.05*sign(params.base)*params.base,
                             max = params.base + 0.05*sign(params.base)*params.base)
        names(params.same) <- names(params.base)
        #test the set
        works <-  is.finite(ls2(par = params.same,
                                init = init.cond,
                                vals = range[k],
                                namez = name.par))
        #update count variable
        count.test <- count.test + 1
      }
    }
    #run optimization only if parameter combination works
    if(works == TRUE){

      opt <- optms(params = params.same, range = range[k])

      #jump to next value if fitting succeded,
      if(opt$value < (max(prs.two[c(k-1, k+1),2]) + no.diff)){
        prs.two[k,] <- c(range[k], opt$value, opt$par)
        saveRDS(prs.two, folder)
        next
      }

    }

    #Lim fit####
    if(lim.number > 0){
      #use left and right parameter sets as intervals and sample from them
      for(r in 1:lim.number){
        print("lim.run")
        limits <- rbind(params.first, params.second)
        params.lim <- runif(length(params.first),
                            min = apply(limits, 2, min),
                            max = apply(limits, 2, max))
        names(params.lim) <- names(prs.two[k-1,par.range[1]:par.range[2]])

        #try and see if parameter combination work
        if(is.finite(ls2(par = params.lim,
                         init = init.cond,
                         vals = range[k],
                         namez = name.par
        ))){

          #optimise
          opt <- optms(params = params.lim, range = range[k])

          #jump to next value if fitting succeded,
          if(opt$value < (max(prs.two[c(k-1, k+1),2]) + no.diff)){
            prs.two[k,] <- c(range[k], opt$value, opt$par)
            saveRDS(prs.two, folder)
            next
          }

        }

      }#for 1:lim.number

    }# if lim.number > 0

    #random fit####
    if(all.number > 0){
      #use min and max of the whole data set as intervals and sample from them
      for(r in 1:all.number){
        print("all.run")
        params.all <- runif(length(params.first),
                            min = apply(prs.two[,par.range[1]:par.range[2]], 2, min),
                            max = apply(prs.two[,par.range[1]:par.range[2]], 2, max))
        names(params.all) <- names(prs.two[k-1,par.range[1]:par.range[2]])

        #try and see if parameter combination work
        if(is.finite(ls2(par = params.all,
                         init = init.cond,
                         vals = range[k],
                         namez = name.par
        ))){

          #optimise
          opt <- optms(params = params.all, range = range[k])

          #jump to next value if fitting succeded,
          if(opt$value < (max(prs.two[c(k-1, k+1),2]) + no.diff)){
            prs.two[k,] <- c(range[k], opt$value, opt$par)
            saveRDS(prs.two, folder)
            next
          }

        }

      }#for r in all.number

    }#if all.number > 0

  }

  #update badfit vector to continue while loop if needed
  badfit <- rep(0, length(range))
  for(i in 2:(length(range)-1)){
    if((prs.two[i,2] > (prs.two[i-1,2] + no.diff) && #check if higher than left value
        prs.two[i,2] > (prs.two[i+1,2] + no.diff)) #|| #check if higher than right value
    ){
      badfit[i] <- 1
    }
  }
  badfit <- which(badfit==1)
  print(paste0("Spiked values are found in:"))
  print(badfit)
  limit1 <- limit1 + 1
  }

}#end of while(improvement > 0 || length(badfit) > 0){

}
