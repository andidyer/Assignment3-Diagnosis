IndependentVariableBoolean = function(datacolumn){
  pEvent = length(grep(1, datacolumn)) / length(datacolumn)
  qEvent = 1 - pEvent
  distribution = c(pEvent, qEvent); names(distribution) = c('p','q')
  
  out = matrix(distribution,1,2); colnames(out) = c('p','q')
  return(out)
}

IndependentVariable = function(datacolumn){
  outlist = list(
    pIndices = grep(1, datacolumn),
    qIndices = grep(0, datacolumn),
    probabilities = matrix(0, 1, 2)
  )
  colnames(outlist$probabilities) = c('p','q')
  
  #'This is the actual probability table, which will be used in the network
  #'if this is an independent variable (in this case, Sm, VTB and Pn are IV)
  outlist$probabilities[,'p'] = length(outlist$pIndices) / length(datacolumn)
  outlist$probabilities[,'q'] = length(outlist$qIndices) / length(datacolumn)
  return(outlist)
}

learn = function(hist){
  #Find indices for each
  Pn = IndependentVariable(hist$Pn)
  VTB = IndependentVariable(hist$VTB)
  TB = IndependentVariable(hist$TB)
  Sm = IndependentVariable(hist$Sm)
  LC = IndependentVariable(hist$LC)
  Br = IndependentVariable(hist$Br)
  XR = IndependentVariable(hist$XR)
  Dy = IndependentVariable(hist$Dy)
  
  #Br|Sm
  Br.Sm = matrix(0, 2, 3); colnames(Br.Sm) = c('Sm','p','q')
  Br.Sm[,'Sm'] = c(0,1)
  Br.Sm[1,'p'] = length(grep(1, hist[c(Sm$qIndices),'Br'])) / length(Sm$qIndices)
  Br.Sm[2,'p'] = length(grep(1, hist[c(Sm$pIndices),'Br'])) / length(Sm$pIndices)
  #CAN THE BELOW JUST BE 1 - p?
  Br.Sm[1,'q'] = length(grep(0, hist[c(Sm$qIndices),'Br'])) / length(Sm$qIndices)
  Br.Sm[2,'q'] = length(grep(0, hist[c(Sm$pIndices),'Br'])) / length(Sm$pIndices)
  
  #LC|Sm
  LC.Sm = matrix(0, 2, 3); colnames(LC.Sm) = c('Sm','p','q')
  LC.Sm[,'Sm'] = c(0,1)
  LC.Sm[1,'p'] = length(grep(1, hist[c(Br$qIndices),'LC'])) / length(Br$qIndices)
  LC.Sm[2,'p'] = length(grep(1, hist[c(Br$pIndices),'LC'])) / length(Br$pIndices)
  LC.Sm[1,'q'] = length(grep(0, hist[c(Br$qIndices),'LC'])) / length(Br$qIndices)
  LC.Sm[2,'q'] = length(grep(0, hist[c(Br$pIndices),'LC'])) / length(Br$pIndices)
  
  
  #TB|VTB
  TB.VTB = matrix(0, 2, 3); colnames(TB.VTB) = c('VTB','p','q')
  TB.VTB[,'VTB'] = c(0,1)
  TB.VTB[1,'p'] = length(grep(1, hist[c(VTB$qIndices),'TB'])) / length(VTB$qIndices)
  TB.VTB[2,'p'] = length(grep(1, hist[c(VTB$pIndices),'TB'])) / length(VTB$pIndices)
  TB.VTB[1,'q'] = length(grep(0, hist[c(VTB$qIndices),'TB'])) / length(VTB$qIndices)
  TB.VTB[2,'q'] = length(grep(0, hist[c(VTB$pIndices),'TB'])) / length(VTB$pIndices)
  
  #'NEXT ONES ALL HAVE MORE THAN ONE CAUSE
  #'Will need to use set operations
  
  #Dy|Br,LC
  Dy.Br_LC = matrix(0, 4, 4); colnames(Dy.Br_LC) = c('Br','LC','p','q')
  Dy.Br_LC[,'LC'] = c(0,1,0,1)
  Dy.Br_LC[,'Br'] = c(0,0,1,1)
  #NOT A OR B [1,]
  jointIndices = intersect(Br$qIndices,LC$qIndices)
  Dy.Br_LC[1,'p'] = length(grep(1, hist[jointIndices,'Dy'])) / length(jointIndices)
  Dy.Br_LC[1,'q'] = 1 - Dy.Br_LC[1,'p']
  #NOT A AND B [2,]
  jointIndices = setdiff(LC$pIndices,Br$pIndices)
  Dy.Br_LC[2,'p'] = length(grep(1, hist[jointIndices,'Dy'])) / length(jointIndices)
  Dy.Br_LC[2,'q'] = 1 - Dy.Br_LC[2,'p']
  #A AND NOT B [3,]
  jointIndices = setdiff(Br$pIndices,LC$pIndices)
  Dy.Br_LC[3,'p'] = length(grep(1, hist[jointIndices,'Dy'])) / length(jointIndices)
  Dy.Br_LC[3,'q'] = 1 - Dy.Br_LC[3,'p']
  #A AND B [4,]
  jointIndices = intersect(Br$pIndices,LC$pIndices)
  Dy.Br_LC[4,'p'] = length(grep(1, hist[jointIndices,'Dy'])) / length(jointIndices)
  Dy.Br_LC[4,'q'] = 1 - Dy.Br_LC[4,'p']
  
  #XR|Pn,LC,TB
  XR.Pn_LC_TB = matrix(0, 8, 5); colnames(XR.Pn_LC_TB) = c('Pn','LC','TB','p','q')
  XR.Pn_LC_TB[,'TB'] = c(0,1,0,0,1,0,1,1)
  XR.Pn_LC_TB[,'LC'] = c(0,0,1,0,1,1,0,1)
  XR.Pn_LC_TB[,'Pn'] = c(0,0,0,1,0,1,1,1)
  #-Pn -LC -TB [1,]
  jointIndices = union(union(TB$pIndices,LC$pIndices),Pn$pIndices)
  XR.Pn_LC_TB[1,'p'] = length(grep(1, hist[-jointIndices,'XR'])) / (nrow(hist) -length(jointIndices))
  XR.Pn_LC_TB[1,'q'] = 1 - XR.Pn_LC_TB[1,'p']
  #-Pn -LC TB [2,]
  jointIndices = setdiff(TB$pIndices,LC$pIndices)
  jointIndices = setdiff(jointIndices,Pn$pIndices)
  XR.Pn_LC_TB[2,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[2,'q'] = 1 - XR.Pn_LC_TB[2,'p']
  #-Pn LC -TB [3,]
  jointIndices = setdiff(LC$pIndices,TB$pIndices)
  jointIndices = setdiff(jointIndices,Pn$pIndices)
  XR.Pn_LC_TB[3,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[3,'q'] = 1 - XR.Pn_LC_TB[3,'p']
  #Pn -LC -TB [4,]
  jointIndices = setdiff(Pn$pIndices,TB$pIndices)
  jointIndices = setdiff(jointIndices,LC$pIndices)
  XR.Pn_LC_TB[4,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[4,'q'] = 1 - XR.Pn_LC_TB[4,'p']
  #-Pn LC TB [5,]
  jointIndices = intersect(LC$pIndices,TB$pIndices)
  jointIndices = setdiff(jointIndices,Pn$pIndices)
  XR.Pn_LC_TB[5,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[5,'q'] = 1 - XR.Pn_LC_TB[5,'p']
  #Pn LC -TB [6,]
  jointIndices = intersect(LC$pIndices,Pn$pIndices)
  jointIndices = setdiff(jointIndices,TB$pIndices)
  XR.Pn_LC_TB[6,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[6,'q'] = 1 - XR.Pn_LC_TB[6,'p']
  #Pn -LC TB [7,]
  jointIndices = intersect(TB$pIndices,Pn$pIndices)
  jointIndices = setdiff(jointIndices,LC$pIndices)
  XR.Pn_LC_TB[7,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[7,'q'] = 1 - XR.Pn_LC_TB[6,'p']
  #Pn LC TB [8,]
  jointIndices = intersect(intersect(TB$pIndices,LC$pIndices),Pn$pIndices)
  XR.Pn_LC_TB[8,'p'] = length(grep(1, hist[jointIndices,'XR'])) / length(jointIndices)
  XR.Pn_LC_TB[8,'q'] = 1 - XR.Pn_LC_TB[8,'p']
  
  
  #Finally, work out what to do about the temperature distribution
  
  
  #Build network, initialise as empty
  
  return(list(BR = Br.Sm,LC = LC.Sm,TB = TB.VTB,Dy = Dy.Br_LC, XR = XR.Pn_LC_TB))
  
  print('boo')
}