
###############################################
#
#  Stutt et al. 2020
#  "A modelling framework to assess the likely effectiveness of facemasks 
# in combination with 'lock-down' in managing the COVID-19 pandemic"
#  Code for Model 2
#
###############################################


library(deSolve)

options(stringsAsFactors = FALSE)

#Note that all tX parameters are in fact defined and used as rate parameters (i.e. 1/T) and are not timescales
#Model time units are days
createInoculumPool = function() {
  
  inoculum = list(
    D=0,
    F=0
  )
  
  parameters = list(
    tD = 1.0/(10/(24*3600)), #10 seconds
    tF = 1.0/2.0
  )
  
  return(list(inoculum = inoculum,
              parameters = parameters))
}

createPopParameters = function(bD = 1.0, bF = 1.0, mI = 1.0, mD = 1.0, mF = 1.0) {
  
  epi = list(
    tE = 1.0 / 3.8,
    tC = 1.0 / 1.2,
    tI = 1.0 / 3.2,
    bC = 2.71,#cryptic infectiousness multiple of bI
    bI = 1.0,#w.l.o.g.
    bD = bD,
    bF = bF,
    mI = mI,
    mC = mI,#Intentional
    mD = mD,
    mF = mF
  )
  
  return(epi)
}

createPopulation = function(name, popState, pop_bD, pop_bF, pop_mI = 1.0, pop_mD = 1.0, pop_mF = 1.0) {
  
  state = list(S = 0, 
               E = 0, 
               C = 0, 
               I = 0, 
               R = 0)
  
  for(cName in names(popState)) {
    state[[cName]] = popState[[cName]]
  }
  
  parameters = createPopParameters(bD = pop_bD, bF = pop_bF, mI = pop_mI, mD = pop_mD, mF = pop_mF)
  
  return(list(name = name,
              state = state, 
              parameters = parameters))
}

createPandemicManagementParameters = function(tFirstQuarantine = 0.0, durationQuarantineCycle, durationBetweenQuarantineCycles, nQuarantineCycles = 1, quarantineReductionInGoingOut = 1.0) {
  return(list(
    tFirstQuarantine = tFirstQuarantine, 
    durationQuarantineCycle = durationQuarantineCycle,
    durationBetweenQuarantineCycles = durationBetweenQuarantineCycles,
    nQuarantineCycles = nQuarantineCycles, 
    quarantineReductionInGoingOut = quarantineReductionInGoingOut
  ))
}

noQuarantine = function() {
  return(createPandemicManagementParameters(nQuarantineCycles = 0, durationQuarantineCycle = 0, durationBetweenQuarantineCycles = 0))
}

isInQuarantine = function(t, pandemicManagementParams) {
  
  if(t < pandemicManagementParams$tFirstQuarantine) {
    return(FALSE)  
  }
  
  quarantineSinglePeriod = pandemicManagementParams$durationQuarantineCycle + pandemicManagementParams$durationBetweenQuarantineCycles
  
  quarantineEndTIme = pandemicManagementParams$tFirstQuarantine + quarantineSinglePeriod * pandemicManagementParams$nQuarantineCycles
  
  if(t >= quarantineEndTIme) {
    return(FALSE)  
  }
  
  #fmod type arithmetic
  nQuarantinePeriod = floor((t - pandemicManagementParams$tFirstQuarantine) / quarantineSinglePeriod)
  
  tInPeriod = (t - pandemicManagementParams$tFirstQuarantine) - nQuarantinePeriod * quarantineSinglePeriod
  
  if(tInPeriod < pandemicManagementParams$durationQuarantineCycle) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}

plotQuarantinePeriods = function(pandemicManagementParams, shadingColour = rgb(red = 0,green = 0,blue = 0,alpha = 0.1)) {
  yrange = c(0, 1e9)#Just absurdly large to make sure it runs the height of the graph
  
  qPeriodTotal = pandemicManagementParams$durationQuarantineCycle + pandemicManagementParams$durationBetweenQuarantineCycles
  
  for(iQPeriod in seq_len(pandemicManagementParams$nQuarantineCycles)) {
    tStart = pandemicManagementParams$tFirstQuarantine + (iQPeriod - 1) * qPeriodTotal
    tEnd = tStart + pandemicManagementParams$durationQuarantineCycle
    
    rect(xleft = tStart, xright = tEnd, ybottom = yrange[1], ytop = yrange[2], col = shadingColour)
  }
}

runModel = function(initialModelState, timeSpan, pandemicManagement = createPandemicManagementParameters(nQuarantineCycles = 0, durationQuarantineCycle = 0, durationBetweenQuarantineCycles = 0)) {
  
  #Assemble state vector
  modelState = c()
  
  #Populations:
  for(population in initialModelState$populations) {
    for(cName in names(population$state)) {
      globalClassName = paste0(population$name, cName)
      modelState[[globalClassName]] = population$state[[cName]]
    }
  }
  
  #Inoculum
  for(iName in names(initialModelState$inoculumPool$inoculum)) {
    modelState[[iName]] = initialModelState$inoculumPool$inoculum[[iName]]
  }
  
  
  #Assemble parameters vector
  modelParams = c()
  
  #Populations:
  for(population in initialModelState$populations) {
    for(pName in names(population$parameters)) {
      globalParamName = paste0(population$name, pName)
      modelParams[[globalParamName]] = population$parameters[[pName]]
    }
  }
  
  #Inoculum
  for(pName in names(initialModelState$inoculumPool$parameters)) {
    modelParams[[pName]] = initialModelState$inoculumPool$parameters[[pName]]
  }
  
  #Create derivative function
  #This isn't too slow, but could probably be enhanced with some index packing to avoid excessive string manipulation if needed
  #We are capturing initialModelState and pandemicManagement in this (budget) lambda
  dModel = function(t, state, parameters) {
    
    #Obtain a named vector guaranteed to be in the same order as the state vector
    dF = state * 0
    
    isQuarantine = isInQuarantine(t = t, pandemicManagementParams = pandemicManagement)
    
    qF = 1.0
    if(isQuarantine) {
      qF = 1.0 - pandemicManagement$quarantineReductionInGoingOut
    }
    
    #Populations:
    totalDropletProductionRate = 0.0
    for(population in initialModelState$populations) {
      #S
      gS = paste0(population$name, "S")
      
      gmD = paste0(population$name, "mD")
      gbD = paste0(population$name, "bD")
      sLossD = parameters[[gbD]] * parameters[[gmD]] * qF * state[["D"]] * state[[gS]]
      
      gmF = paste0(population$name, "mF")
      gbF = paste0(population$name, "bF")
      sLossF = parameters[[gbF]] * parameters[[gmF]] * qF * state[["F"]] * state[[gS]]
      
      sLoss = sLossD + sLossF
      
      dF[[gS]] = -sLoss
      
      #E
      gE = paste0(population$name, "E")
      
      gtE = paste0(population$name, "tE")
      eLoss = parameters[[gtE]] * state[[gE]]
      
      dF[[gE]] = sLoss - eLoss
      
      #C
      gC = paste0(population$name, "C")
      
      gtC = paste0(population$name, "tC")
      cLoss = parameters[[gtC]] * state[[gC]]
      
      gbC = paste0(population$name, "bC")
      gmC = paste0(population$name, "mC")
      dProdByC = parameters[[gbC]] * parameters[[gmC]] * qF * state[[gC]]
      
      totalDropletProductionRate = totalDropletProductionRate + dProdByC
      
      dF[[gC]] = eLoss - cLoss
      
      #I
      gI = paste0(population$name, "I")
      
      gtI = paste0(population$name, "tI")
      iLoss = parameters[[gtI]] * state[[gI]]
      
      gbI = paste0(population$name, "bI")
      gmI = paste0(population$name, "mI")
      dProdByI = parameters[[gbI]] * parameters[[gmI]] * qF * state[[gI]]
      
      totalDropletProductionRate = totalDropletProductionRate + dProdByI
      
      dF[[gI]] = cLoss - iLoss
      
      #R
      gR = paste0(population$name, "R")
      
      dF[[gR]] = iLoss
      
    }
    
    #Inoculum
    dLoss = parameters[["tD"]] * state[["D"]]
    dF[["D"]] = totalDropletProductionRate - dLoss
    dF[["F"]] = dLoss - parameters[["tF"]] * state[["F"]]
    
    return(list(dF))
  }
  
  
  #Output timepoints
  tStart = 0.0
  tEnd = timeSpan
  if(length(timeSpan) == 2) {
    tStart = timeSpan[1]
    tEnd = timeSpan[2]
  }
  
  modelTimes = seq(from = tStart, to = tEnd, by = 0.1)#Note timestep here doesn't affect integration timestep, this is purely for the output
  
  #Actually solve:
  odeOut = ode(y = modelState, times = modelTimes, func = dModel, parms = modelParams)
  
  
  #Package output back up into same structure as input
  stateOut = list(time = odeOut[, "time"])
  
  #Populations:
  for(population in initialModelState$populations) {
    for(cName in names(population$state)) {
      globalClassName = paste0(population$name, cName)
      stateOut[[population$name]][[cName]] = odeOut[, globalClassName]
    }
  }
  
  #Inoculum
  for(iName in names(initialModelState$inoculumPool$inoculum)) {
    stateOut[["inoculum"]][[iName]] = odeOut[, iName]
  }
  
  return(stateOut)
}

fitModelParamsFromAnalytic = function(assumedDropletInfectionProportion = 0.5, assumedR0 = 2.5, populationSize = 6e7) {
  #Note that all the tX parameters are actually _rates_, so need to use 1/tX wherever it appears in an analytic expression
  
  if(assumedDropletInfectionProportion > 1.0 | assumedDropletInfectionProportion < 0.0) {
    stop("assumedDropletInfectionProportion must be in range [0.0, 1.0]")
  }
  
  baselineEpiParams = createPopParameters()
  
  baselineInoculumParams = createInoculumPool()$parameters
  
  #This quantity appears in all the rate parameters
  coreRatio = assumedR0 / (((baselineEpiParams$bC / baselineEpiParams$tC) + (baselineEpiParams$bI / baselineEpiParams$tI)) * populationSize)
  
  params = list(
    bD = coreRatio *        assumedDropletInfectionProportion  / (1.0 / baselineInoculumParams$tD),
    bF = coreRatio * (1.0 - assumedDropletInfectionProportion) / (1.0 / baselineInoculumParams$tF)
  )
  
  return(params)
}

runScenario = function(modelParameters, maskParameters, quarantineParameters, initialExposed, populationSize, timespan) {
  
  #Checking the mask wearing start time - don't want it to be longer than the simulation length, or things will go wonky
  if(maskParameters$startTime > timespan) {
    stop(paste0("Specified start time for mask wearing (", maskParameters$startTime, ") was later than the end of the simulation (", timespan, ")"))
  }
  
  #Pre mask
  modelPopulationsPreMask = list(createPopulation(name = "NoMask", 
                                           popState = list(S = populationSize, E = initialExposed), 
                                           pop_bD = modelParameters$bD, pop_bF = modelParameters$bF)
  )
  
  modelStatePreMask = list(populations = modelPopulationsPreMask,
                           inoculumPool = createInoculumPool())
  
  if(maskParameters$startTime > 0.0) {
    modelResultsPreMask = runModel(initialModelState = modelStatePreMask,
                                   timeSpan = maskParameters$startTime, 
                                   pandemicManagement = quarantineParameters)
  } else {
    modelResultsPreMask = list(time = 0.0,
                               NoMask = modelPopulationsPreMask[[1]]$state,
                               inoculum = modelStatePreMask$inoculumPool$inoculum)
  }
  
  iFinalTime = length(modelResultsPreMask$time)
  #Allocate people randomly as mask and non-mask wearers
  #TODO: Could have people preferentially wearing masks if they have symptoms, 
  #but this is unlikely to make a difference as the effect is very transient due to only affecting individuals at the time masks are first adopted
  popNoMask = list()
  popMask = list()
  for(cName in names(modelResultsPreMask$NoMask)) {
    classFinalValue = modelResultsPreMask$NoMask[[cName]][iFinalTime]
    
    popNoMask[[cName]] = (1.0 - maskParameters$proportion) * classFinalValue
    popMask[[cName]] = maskParameters$proportion * classFinalValue
  }
  
  #Post mask
  modelPopulationsPostMask = list(createPopulation(name = "NoMask", 
                                           popState = popNoMask, 
                                           pop_bD = modelParameters$bD, pop_bF = modelParameters$bF),
                          createPopulation(name = "Mask", 
                                           popState = popMask, 
                                           pop_bD = modelParameters$bD, pop_bF = modelParameters$bF, 
                                           pop_mI = 1.0 - maskParameters$exhaleCapture, pop_mD = 1.0 - maskParameters$inhaleCapture, pop_mF = 1.0 - maskParameters$fomiteCapture)
  )
  
  modelInoculumPoolPostMask = createInoculumPool()
  for(iName in names(modelResultsPreMask$inoculum)) {
    modelInoculumPoolPostMask$inoculum[[iName]] = modelResultsPreMask$inoculum[[iName]][iFinalTime]
  }
  
  modelStatePostMask = list(populations = modelPopulationsPostMask,
                            inoculumPool = modelInoculumPoolPostMask)
  
  modelResultsPostMask = runModel(initialModelState = modelStatePostMask, 
                                  timeSpan = c(maskParameters$startTime, timespan), 
                                  pandemicManagement = quarantineParameters)
  
  modelTime = c(modelResultsPreMask$time, modelResultsPostMask$time)
  
  noMaskResults = list()
  maskResults = list()
  for(cName in names(modelResultsPreMask$NoMask)) {
    noMaskResults[[cName]] = c(      modelResultsPreMask$NoMask[[cName]], modelResultsPostMask$NoMask[[cName]])
    maskResults[[cName]]   = c(0.0 * modelResultsPreMask$NoMask[[cName]], modelResultsPostMask$Mask[[cName]])#Mask wearer stats are zero before they exist
  }
  
  inoculumResults = list(D=c(modelResultsPreMask$inoculum$D, modelResultsPostMask$inoculum$D), 
                         F=c(modelResultsPreMask$inoculum$F, modelResultsPostMask$inoculum$F))
  
  runResults = list(time = modelTime,
                    noMaskResults = noMaskResults,
                    maskResults = maskResults,
                    inoculumResults = inoculumResults)
  
  scenarioData = list(runResults = runResults,
                      maskParameters = maskParameters,
                      quarantineParameters = quarantineParameters,
                      initialExposed = initialExposed, 
                      populationSize = populationSize)
}

getMaskString = function(maskParams, omitProportion=FALSE) {
  
  maskString = "No Masks"
  
  if(maskParams$proportion > 0.0 | omitProportion == TRUE) {#omitProportion probably wants us to ignore No Masks output?
    maskString = paste0("Masks: ")
    if(maskParams$startTime > 0) {
      maskString = paste0(maskString, "Start: T+", maskParams$startTime, " days, ")
    }
    if(omitProportion != TRUE) {
      maskString = paste0(maskString, sprintf("%.0f%%", 100.0 * maskParams$proportion), " wearing, ")
    }
    maskString = paste0(maskString, "inoculum capture: ", sprintf("%.0f%%", 100.0 * maskParams$exhaleCapture), " exhale, ", sprintf("%.0f%%", 100.0 * maskParams$inhaleCapture), " inhale, fomite ", sprintf("%.0f%%", 100.0 * maskParams$fomiteCapture))
  }
  
  return(maskString)
}

getQuarantineString = function(quarantineParams) {
  
  quarantineString = "No Lock-down"
  
  if(quarantineParams$nQuarantineCycles > 0) {
    quarantineString = paste0("Lock-down: Start: T+", quarantineParams$tFirstQuarantine, " days, duration: ", quarantineParams$durationQuarantineCycle, " days, ", sprintf("%.0f%%", 100.0 * quarantineParams$quarantineReductionInGoingOut), " contact reduction, ", quarantineParams$durationBetweenQuarantineCycles, " days off, ", quarantineParams$nQuarantineCycles, " cycles")
  }
  
  return(quarantineString)
}

plotLine = function(x, y, col, lty, lwd, legendString) {
  lines(x = x, y = y, col = col, lty = lty, lwd = lwd)
  
  legendData = list(string = legendString,
                    col = col,
                    lty = lty,
                    lwd = lwd)
  
  return(legendData)
}

plotScenario = function(scenarioData) {
  
  maskString = getMaskString(scenarioData$maskParameters)
  
  quarantineString = getQuarantineString(scenarioData$quarantineParameters)
  
  #Show y axis as cases per million, not absolute nummbers
  caseScaling = 1e-6
  
  totalResults = list()
  for(cName in names(scenarioData$runResults$noMaskResults)) {
    scenarioData$runResults$noMaskResults[[cName]] = scenarioData$runResults$noMaskResults[[cName]] * caseScaling
    scenarioData$runResults$maskResults[[cName]]   = scenarioData$runResults$maskResults[[cName]]   * caseScaling
    totalResults[[cName]] = scenarioData$runResults$noMaskResults[[cName]] + scenarioData$runResults$maskResults[[cName]]
  }
  
  yMax = scenarioData$populationSize * caseScaling
  while(yMax > 20.0*max(totalResults$R)) {
    yMax = yMax / 10.0
  }
  
  #Probably would have been better to use ggplot
  plot(type = "n", 
       x = -1,
       xlim = range(scenarioData$runResults$time),
       ylim = c(0, yMax), 
       xlab = paste0("Time since outbreak first reached ", scenarioData$initialExposed, " infected individuals (days)"), 
       ylab = "Cases (millions)", 
       main = paste0(maskString, "\n", quarantineString), 
       cex.main = 2, cex.lab = 1.9, cex.axis = 2)#TODO: Larger axis labels
  
  legendData = plotLine(x = scenarioData$runResults$time, y = totalResults$R, col = "black", lty = 1, lwd = 2, legendString = "Total Removed")
  
  ld = plotLine(x = scenarioData$runResults$time, y = totalResults$I, col = "red", lty = 1, lwd = 2, legendString = "Total active symptomatic")
  for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
  
  if(scenarioData$maskParameters$proportion > 0.0 & scenarioData$maskParameters$proportion < 1.0) {
    
    ld = plotLine(x = scenarioData$runResults$time, y = scenarioData$runResults$maskResults$R, col = "grey", lty = 2, lwd = 2, legendString = "Mask wearers removed")
    for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
    
    ld = plotLine(x = scenarioData$runResults$time, y = scenarioData$runResults$noMaskResults$R, col = "black", lty = 3, lwd = 2, legendString = "Non-mask wearers removed")
    for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
  }
  
  #Inoculum:
  #scaledInoc = scenarioData$runResults$inoculumResults$F * yMax / max(scenarioData$runResults$inoculumResults$F)
  #ld = plotLine(x = scenarioData$runResults$time, y = scaledInoc, col = "blue", lty = 1, lwd = 2, legendString = "Scaled fomite inoculum")
  #for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
  
  plotQuarantinePeriods(pandemicManagementParams = scenarioData$quarantineParameters)
  
  #Date mask start usage:
  if(scenarioData$maskParameters$proportion > 0.0) {
    abline(v = scenarioData$maskParameters$startTime, col = "blue", lty = 3, lwd = 3)
  }
  
  legend("topleft", 
         legend = legendData$string, 
         col = legendData$col, 
         lty = legendData$lty,
         lwd = legendData$lwd,
         cex = 2)
}


#Simulation settings:
individualPlotSize = c(1920*0.5, 1080*0.5)

popSize = 6e7
fittedModelParameters = fitModelParamsFromAnalytic(assumedR0 = 4.0, assumedDropletInfectionProportion = 0.5, populationSize = popSize)

modelTimespan = 540
modelInitialExposed = 100

#Baseline parameters:
maskEffectivenessParameters = list(exhaleCapture = 0.5, inhaleCapture = 0.5, fomiteCapture = 0.0)
noMasks = list(
  exhaleCapture = 0.0,
  inhaleCapture = 0.0,
  fomiteCapture = 0.0,
  startTime = 0.0,
  proportion = 0.0)

quarantineParameters = createPandemicManagementParameters(tFirstQuarantine = 45,
                                                          durationQuarantineCycle = 90,
                                                          durationBetweenQuarantineCycles = 30,
                                                          nQuarantineCycles = 4,
                                                          quarantineReductionInGoingOut = 0.5)


#Panel: quarantine + various levels of mask wearing
masks = list(noMasks)

maskProportions = c(0.25, 0.50, 1.00)

for(maskProportion in maskProportions) {
  maskParams = list(
    exhaleCapture = maskEffectivenessParameters$exhaleCapture,
    inhaleCapture = maskEffectivenessParameters$inhaleCapture,
    fomiteCapture = maskEffectivenessParameters$fomiteCapture,
    startTime = 0.0,
    proportion = maskProportion)
  
  masks[[length(masks)+1]] = maskParams
}

quarantines = list()

qStarts = c(45, 45, 45, 45)
#qStarts = c(45, 50, 55, 100)

for(qStart in qStarts) {
  
  qParams = quarantineParameters
  
  qParams$tFirstQuarantine = qStart
  
  quarantines[[length(quarantines)+1]] = qParams
}

mqs = list()#TODO: function called zip does this?
for(iMask in seq_len(length(masks))){
  mqs[[length(mqs)+1]] = list(m = masks[[iMask]], q = quarantines[[iMask]])
}

nRows = 2
nCols = ceiling(length(mqs) / nRows)

png(filename = "Paper_5_MaskProportionPanel.png", width = individualPlotSize[1]*nCols, height = individualPlotSize[2]*nRows)
oPar = par(mfrow = c(nRows, nCols),
           oma = c(5,4,0,0) + 0.1,
           mar = c(4,5,4,1) + 0.1)

iFig = 1

for(mq in mqs) {
  
  scenarioResults = runScenario(modelParameters = fittedModelParameters, 
                                maskParameters = mq$m, 
                                quarantineParameters = mq$q, 
                                initialExposed = modelInitialExposed,
                                populationSize = popSize, 
                                timespan = modelTimespan)
  
  plotScenario(scenarioData = scenarioResults)
  
  #Label
  fb = par("usr")
  text(labels = paste0("(", letters[iFig], ")"), x = fb[1] + 0.95 * (fb[2] - fb[1]), y = fb[3] + 0.9 * (fb[4] - fb[3]), cex = 5)
  
  iFig = iFig + 1
}

#Restore old par settings
par(oPar)
dev.off()


#Panel: quarantine + timing of mask wearing before first Q, around peak of I, at end of first Q
masks = list()

maskStartTimes = c(30, 60, 90, 120)

for(maskStartTime in maskStartTimes) {
  maskParams = list(
    exhaleCapture = maskEffectivenessParameters$exhaleCapture,
    inhaleCapture = maskEffectivenessParameters$inhaleCapture,
    fomiteCapture = maskEffectivenessParameters$fomiteCapture,
    startTime = maskStartTime,
    proportion = 1.0)
  
  masks[[length(masks)+1]] = maskParams
}


nRows = 2
nCols = ceiling(length(masks) / nRows)

png(filename = "Paper_6_MaskStartTimePanel.png", width = individualPlotSize[1]*nCols, height = individualPlotSize[2]*nRows)
oPar = par(mfrow = c(nRows, nCols),
           oma = c(5,4,0,0) + 0.1,
           mar = c(4,5,4,1) + 0.1)

iFig = 1

for(mask in masks) {
  scenarioResults = runScenario(modelParameters = fittedModelParameters, 
                                maskParameters = mask,
                                quarantineParameters = quarantineParameters,
                                initialExposed = modelInitialExposed,
                                populationSize = popSize, 
                                timespan = modelTimespan)
  
  plotScenario(scenarioData = scenarioResults)
  
  #Label
  fb = par("usr")
  text(labels = paste0("(", letters[iFig], ")"), x = fb[1] + 0.95 * (fb[2] - fb[1]), y = fb[3] + 0.9 * (fb[4] - fb[3]), cex = 5)
  
  iFig = iFig + 1
}

#Restore old par settings
par(oPar)
dev.off()


#Supp: 2x2 no quarantine + various levels of mask wearing
masks = list(noMasks)

maskProportions = c(0.25, 0.50, 1.00)

for(maskProportion in maskProportions) {
  maskParams = list(
    exhaleCapture = maskEffectivenessParameters$exhaleCapture,
    inhaleCapture = maskEffectivenessParameters$inhaleCapture,
    fomiteCapture = maskEffectivenessParameters$fomiteCapture,
    startTime = 0.0,
    proportion = maskProportion)
  
  masks[[length(masks)+1]] = maskParams
}


nRows = 2
nCols = ceiling(length(masks) / nRows)

png(filename = "Paper_7_MaskProportionPanelNoQuarantine.png", width = individualPlotSize[1]*nCols, height = individualPlotSize[2]*nRows)
oPar = par(mfrow = c(nRows, nCols),
           oma = c(5,4,0,0) + 0.1,
           mar = c(4,5,4,1) + 0.1)

iFig = 1

for(mask in masks) {
  scenarioResults = runScenario(modelParameters = fittedModelParameters, 
                                maskParameters = mask,
                                quarantineParameters = noQuarantine(),
                                initialExposed = modelInitialExposed,
                                populationSize = popSize,
                                timespan = modelTimespan)
  
  plotScenario(scenarioData = scenarioResults)
  
  #Label
  fb = par("usr")
  text(labels = paste0("(", letters[iFig], ")"), x = fb[1] + 0.95 * (fb[2] - fb[1]), y = fb[3] + 0.9 * (fb[4] - fb[3]), cex = 5)
  
  iFig = iFig + 1
}

#Restore old par settings
par(oPar)
dev.off()


#Reduction in individual chance of being infected as a mask or non-mask wearer, given X% of population wearing masks
maskSweep = list(list(maskID = "Positive", maskParams = list(exhaleCapture = 0.5,
                                                             inhaleCapture = 0.5,
                                                             fomiteCapture = 0.0,
                                                             startTime = 0.0,
                                                             proportion = 0.0)),
                 list(maskID = "Neutral", maskParams = list(exhaleCapture = 0.5,
                                                            inhaleCapture = 0.0,
                                                            fomiteCapture = 1.0 - 1/(1-0.5),#To exactly cancel out the exhaleCapture (used as mF = 1 - fomiteCapture, i.e. mF = 1/(1-maskExhale) = 1/mI)
                                                            startTime = 0.0,
                                                            proportion = 0.0)),
                 list(maskID = "Negative", maskParams = list(exhaleCapture = 0.5,
                                                             inhaleCapture = 0.0,
                                                             fomiteCapture = 1.0 - 1/(1-0.75),#To more than cancel out the exhaleCapture (used as mF = 1 - fomiteCapture, i.e. mF = 1/(1-maskExhale) = 1/mI)
                                                             startTime = 0.0,
                                                             proportion = 0.0))
           )

nRows = length(maskSweep)
nCols = 1

png(filename = "Paper_8_MaskBenefitAllocation.png", width = individualPlotSize[1] * nCols, height = individualPlotSize[2] * nRows)
oPar = par(mfrow = c(nRows, nCols),
           oma = c(5,4,0,0) + 0.1,
           mar = c(5,5,7,1) + 0.1)

iFig = 1

for(sweepMask in maskSweep) {

  masks = list()
  
  #Generate masks:
  templateMask = sweepMask$maskParams
  
  maskProportions = seq(from = 0.0, to = 1.00, by = 0.01)
  
  for(maskProportion in maskProportions) {
    maskParams = list(
      exhaleCapture = templateMask$exhaleCapture,
      inhaleCapture = templateMask$inhaleCapture,
      fomiteCapture = templateMask$fomiteCapture,
      startTime = 0.0,
      proportion = maskProportion)
    
    masks[[length(masks)+1]] = maskParams
  }
  
  
  pRNoMasks = c()
  pRMasks = c()
  pRs = c()
  
  for(mask in masks) {
    scenarioResults = runScenario(modelParameters = fittedModelParameters, 
                                  maskParameters = mask,
                                  quarantineParameters = noQuarantine(),
                                  initialExposed = modelInitialExposed,
                                  populationSize = popSize,
                                  timespan = modelTimespan)
    
    iFinalTime = length(scenarioResults$runResults$time)
    
    #Could just do by multiplying popSize * maskProportion...
    totalNoMask = 0
    totalMask = 0
    for(cName in names(scenarioResults$runResults$noMaskResults)) {
      totalNoMask = totalNoMask + scenarioResults$runResults$noMaskResults[[cName]][iFinalTime]
      totalMask = totalMask + scenarioResults$runResults$maskResults[[cName]][iFinalTime]
    }
    
    total = totalNoMask + totalMask
    
    finalRNoMask = scenarioResults$runResults$noMaskResults$R[iFinalTime]
    finalRMask = scenarioResults$runResults$maskResults$R[iFinalTime]
    
    finalR = finalRNoMask + finalRMask
    
    pRNoMask = finalRNoMask / totalNoMask
    pRMask = finalRMask / totalMask
    
    pR = finalR / total
    
    pRs = c(pRs, pR)
    
    #Values don't exist when we have none of a class, NA is more usable than nan
    pRNoMask[is.nan(pRNoMask)] = NA
    pRMask[is.nan(pRMask)] = NA
    
    pRNoMasks = c(pRNoMasks, pRNoMask)
    pRMasks = c(pRMasks, pRMask)
  }
  
  plot(x = -1, type = "n",
       xlim = range(maskProportions),
       ylim = c(min(pRNoMasks, pRMasks, na.rm = TRUE), 1), 
       xlab = "Proportion of population wearing masks", ylab = "Per-capita probability of infection", 
       main = paste0("Allocation of benefits from wearing\n", getMaskString(sweepMask$maskParams, omitProportion=TRUE)),
       cex.main = 3, cex.axis = 2, cex.lab = 3)
  
  legendData = plotLine(x = maskProportions, y = pRNoMasks, col = "black", lty = 1, lwd = 2, legendString = "Non-Mask wearers")
  
  ld = plotLine(x = maskProportions, y = pRMasks, col = "black", lty = 2, lwd = 2, legendString = "Mask Wearers")
  for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
  
  ld = plotLine(x = maskProportions, y = pRs, col = "black", lty = 3, lwd = 2, legendString = "Whole population")
  for(lName in names(legendData)) {legendData[[lName]] = c(legendData[[lName]], ld[[lName]])}
  
  abline(h = pRNoMasks[1], lty = 2, lwd = 2, col = "blue")
  
  legend("bottomleft", 
         legend = legendData$string, 
         col = legendData$col, 
         lty = legendData$lty,
         lwd = legendData$lwd,
         cex = 3)
  
  #Label
  fb = par("usr")
  text(labels = paste0("(", letters[iFig], ")"), x = fb[1] + 0.95 * (fb[2] - fb[1]), y = fb[3] + 0.9 * (fb[4] - fb[3]), cex = 5)
  
  iFig = iFig + 1
}

par(oPar)
dev.off()


