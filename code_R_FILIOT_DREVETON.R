###############################################
########## Importation des packages ###########
###############################################
library(zoo)
library(tseries)
library(ggplot2)
library(scales)
library(fUnitRoots)
library(lm)
library(tseries)
library(urca)
library(forecast)
library(Matrix)
library(xtable)
###############################################
#### Importation et traitement des donnees ####
###############################################


setwd(dir="/Users/bfiliot/Desktop/ENSAE/S2/stprojet")
table = read.csv("valeurs_mensuelles.csv", head=T, sep=";")

names(table)[2] = "Frequentation de passagers vols internationaux - Paris"

table = table[-2,]  # suppression des deux premieres lignes
table = table[-1,]  # elles ne correspondent pas a des valeurs

table$year = substring(table[,1],1,4)  # date au format annee / mois
table$month = substring(table[,1],6,7)

data = table[order(table[4],table[5]),]   # ordonnencement du plus vieux au plus recent

xm = as.zoo(ts(data[[2]],
          start     = c(1994,1),
          end       = c(2018,1),
          frequency = 12)) # transformation en objet zoo

dates = seq(as.POSIXct("1994-01-01"), by = "month", length.out = 289)

###############################################
####### Affichage graphique de la serie #######
###############################################

# Utilisation de ggplot (uniquement ici)

p = ggplot(data = xm, aes(x = dates, y = xm)) +
  geom_line(colour = 'blue') + 
  scale_x_datetime(labels = date_format("%Y"), breaks = date_breaks("years")) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Temps") +
  labs(y = "Fréquentation des passagers (millions)") +
  theme(text = element_text(size=20))
  #labs(title = "Serie temporelle d'origine") +
  
p # graphe de la serie temporelle

p2 = ggplot(data = xm[1:36], aes(x = seq(as.POSIXct("1994-01-01"), by = "month", length.out = 36), y = xm[1:36])) +
  geom_line(colour = 'blue') + 
  scale_x_datetime(labels = date_format("%Y-%m"), breaks = date_breaks("months")) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Temps") +
  labs(y = "Fréquentation des passagers (millions)") 

p2 # graphe de la serie temporelle les 3 premieres annees

###############################################
######## Stationnarisation de la serie ########
###############################################

# Pour ce genre de serie temporelle il est utile
# de transformer en log pour reduire l'heteroscedasticite
# des residus et/ou problemes de non-linearite.

xmlog = log(xm)

# Comme nous l'avons vu, la serie presente une tendance haussiere 
# ainsi qu'une saisonnalite annuelle. La tendance haussiere laisse
# presager d'une serie stationnaire. 

# Tout d'abord verifions la nature de la tendance deterministe.

T = length(xmlog)
time = 1:T
fit1 = lm(xmlog ~ poly(time,1))
fit2 = lm(xmlog ~ poly(time,2))
pred1 = predict(fit1)
pred2 = predict(fit2)
plot(dates, 
     xmlog, 
     type='l', 
     xlab='Temps', 
     ylab='Frequentation des passagers (millions)')
lines(dates, pred1, type = "l", col = "red", lty = 1, lwd = 2)
lines(dates, pred2, type = "l", col = "green", lty = 1, lwd = 2)
legend("bottomright",
       legend = c("Série", "Approx. linéaire", "Approx. quadratique"),
       col = c("black","red", "green"),
       text.col = c("black", "red", "green"),
       lty = c(1,1,1),
       ncol = 1,
       cex = 0.5)

# Conclusion : tendance lineaire en a + bt. On se place donc dans le cas "ct". 
# Notre serie log est-elle stationnaire ? 

adf = adfTest(xmlog, lags=13, type="ct")
adf@test$p.value # < 0.01 
# On rejette fortement la racine unite a 5%, la serie est donc stationnaire
# a 5%. Ce qui finalement n'est pas si contre-intuitif etant donne la tendance
# haussiere deterministe. Le modele est-il valide ? 
# Regardons la validite du modele avec la persistance des residus. 

Qtest = function(x, kmax=24){
  t(apply(matrix(1:kmax), 1, FUN=function(k){
    pv = Box.test(x, lag=k, type="Ljung-Box")$p.value
    return(c("lag"=k,"pval"=pv, "Non autocorreles a 5%" = (pv>0.05)))
  }))
}

resid = adf@test$lm$residuals
plot(resid,type="l")
Qtest(adf@test$lm$residuals)  
# Ici on voit que le test de blancheur des residus est rejete
# a 5% a tous les ordres. Les residus sont correles jusqu'a l'ordre 24.
# Le modele n'est donc pas valide si l'on considere un lag = 0 dans 
# la formule de l'ADF.
# Mais nous pourrions augmenter ce lag de telle sorte a voir a partir duquel
# les residus ne sont pas correles (a 5%) et ceci pour tous les 24 ordres.
# Ceci est verifie pour le lag 23. Mais alors, le test ADF est rejete ! 
# Ce que nous decidons de faire est de garder un test adf avec un lag = 0
# et tester la blancheur des residus a l'ordre 1. Des qu'une transformation
# de la serie xmlog sera telle que 
# le residu est non correle a l'ordre 1 a 5%, on choisira cette transformation
# et on regardera le lag minimum a partir duquel l'ADF/le modele est valide. 

# On prend tout d'abord en compte la saisonnalite.
xdesaison = diff(xmlog,12)
adf = adfTest(xdesaison, lags=0, type="ct")
adf@test$p.value
# On rejette aussi la racine unite a 5%, la serie differenciee est stationnaire
# ce qui est logique.
# On test la blancheur des residus a l'ordre 1 : 
# ici on remarque que la p-valeur est superieure a 5% mais tres proche (0.05417)
resid = adf@test$lm$residuals
plot(resid,type="l")
Box.test(resid, lag=1, type="Ljung-Box")$p.value 
# On decide donc de differencier egalement la serie pour avoir un rejet 
# de la correlation des residus a l'ordre 1 plus fort.

# Appliquons donc une desaisonnalisation et une differenciation.
xdd = diff(desaison,1)
adf = adfTest(xdd, lags=0, type="nc") #lag 11 pour 12 premiers ; lag 15 pour 24 premiers
adf@test$p.value
# On rejette la racine unite a 5%, la serie differenciee et desaisonnalisee
# est stationnaire ce qui est logique.
resid = adf@test$lm$residuals
plot(resid,type="l")
Box.test(resid, lag=1, type="Ljung-Box")$p.value
# Et cette fois l'absence d'autocorrelation n'est pas rejetee a l'ordre 1 a 5%
# car elle est egale a 0.39516. Nous gardons donc cette transformation.
# Verifions a partir de quel lag (dans la formule de l'ADF) le modele est valide,
# i.e tq les residus sont non-autocorreles jusqu'a l'ordre 24 a 5%.
# Au passage on peut verifier qu'a lag 0 le modele n'est pas valide meme
# si pour certains ordres les residus ne sont pas autocorreles a 5%.
Qtest(adf@test$lm$residuals)

lagmin = function(serie){
  i = 0
  if (serie == xdd){
    type = "nc"} # la serie desaisonnalite et differenciee est centree 
  if (serie == xdesaison){
    type = "c"} # la serie desaisonnalite est non centree
  if (serie == xmlog){
    type = "ct"} # la serie xmlog presente une tendance lineaire deterministe
  adf = adfTest(serie, lags=i, type=type)
  bbxtest = Qtest(adf@test$lm$residuals)
  while (sum(bbxtest[,3]) < 24) {
    i = i+1
    adf = adfTest(serie, lags=i, type=type)
    bbxtest = Qtest(adf@test$lm$residuals)
  }
  return("lagminADF"=i)
}

lagmin(xdd)
# On trouve un lag min de 13 ce qui est bien mieux que le 23 pour
# la serie non transformee.
# On a donc un modele ADF valide lorsque celui-ci prend en compte
# 13 retards. 
adfTest(xdd, lags=13, type="nc")
pp.test(xdd) # L'hypothese de racine unitaire est tres fortement rejetee.
kpss.test(xdd) # L'hypothese de stationnarite n'est pas rejetee a 5%. 

# On se contentera donc de notre transformation dans la mesure ou
# elle garantit la non-autocorrelation des residus jusqu'a l'ordre 24
# pour un lagADF egal a 13, tout ceci pour des seuils de 5%. 
# En soit nous aurions juste desaisonnaliser la serie... 
# Car on remarque que le lag minimal de l'ADF pour lequel les residus
# sont non-autocorreles a tous les ordres est...
lagmin(xdesaison)
adfTest(xdesaison, lags=13, type="c")
pp.test(xdesaison) # L'hypothese de racine unitaire est tres fortement rejetee.
kpss.test(xdesaison)
#... 13, comme pour la serie desaisonnalisee. On retiendra donc qu'a 
# lagADF = 0, les residus de la serie desaisonnalisee sont "presque" 
# auto-correles a l'ordre 1 au niveau 5% (le terme "presque" est 
# peu rigoureux) tandis que non pour la serie xdd desaisonnalisee et 
# differenciee.

# Difference avant / apres
par(mfrow=c(2,1))
plot(dates,
     xm, 
     type='l',
     xlab="Temps", 
     ylab = "",
     #ylab='Logarithme de la fréquentation des passagers (millions)',
     col='blue',
     main='Série originale')
T = length(xm)
time = 1:T
fitxm = lm(xm ~ poly(time,1))
predxm = predict(fitxm)
lines(dates, predxm, type = "l", col = "black", lty = 1, lwd = 1)

plot(dates[14:289],
     xdd, 
     type='l',
     xlab="Temps", 
     ylab = "",
     #ylab='Logarithme de la fréquentation des passagers (millions)',
     col='blue',
     main='Série transformée')

abline(a=0, b=0)

###############################################
########## Calibrage d'un modèle ARMA #########
###############################################

########### Fonctions prealables ##############

# Fonction de test de l'absence d'autocorrelation des residus
Qtest = function(x, kmax=24){
  t(apply(matrix(1:kmax), 1, FUN=function(k){
    pv = Box.test(x, lag=k, type="Ljung-Box")$p.value
    return(c("lag"=k,"pval"=pv, "Non autocorreles a 5%" = (pv>0.05)))
  }))
}

# Fonction de test des significativites individuelles des coefficients
signif = function(estim){ 
  coef = estim$coef
  se = sqrt(diag(estim$var.coef))
  t = coef/se
  pval = (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

# Fonction annexe
tf_binaire = function(signif){
  if (is.na(signif)) return(1)
  else return(signif)
}

# Selection des valeurs pmax, qmax, Pmax, Qmax #
# Pour cela nous sommes alles voir comment sont calcules les 
# limites de l'intervalle de confiance a 95% pour l'acf (et donc pacf)
# via la fonction 
getS3method("plot", "acf")
# et nous avons retrouve la variable clim egale a la borne 
par(mfrow=c(2,1))
clim95 = qnorm((1 + 0.95)/2)/sqrt(length(xdd))
clim99 = qnorm((1 + 0.99)/2)/sqrt(length(xdd))
acf = acf(xdd, lag.max = 300, main = "ACF série transformée")
abline(clim99,0,col="red", lty=2)
abline(-clim99,0,col="red", lty=2)
pacf = pacf(xdd, lag.max = 300, main = "PACF série transformée")
abline(clim99,0,col="red", lty=2)
abline(-clim99,0,col="red", lty=2)
# Les 2 valeurs coincident.

which(abs(acf$acf) > clim95)
which(abs(pacf$acf) > clim95)

# On remarque que les valeurs des autocorrelations pour des lag
# entre q=2 et q=11 sont dans l'intervalle de confiance a 95%  donc on peut
# supposer legitiment que qmax = 3.
# Pour la PACF, c'est un peu plus delicat. En effet les autocorrelations
# partielles de lag 1,2,4,5,7,11,12,13,21 ne sont pas dans l'intervalle de 
# confiance a 95%. On peut considerer que l'ecart entre lag 13 et lag 21
# permet de prendre pmax = 14.
# Au niveau des ordres Pmax et Qmax, qu'en est-il ? 
# Etant donne que les mutliples de 12 sont : 12,24,36,48,60,72,84,96,108,120,132,144
# on peut prendre Pmax = 4 et on pourrait aller jusque Qmax = 12...
# Si l'on regarde avec un intervalle de confiance a 99%...

clim99 = qnorm((1 + 0.99)/2)/sqrt(length(xdd))
par(mfrow=c(2,1))
acf = acf(xdd, lag.max = 300, main = "ACF série transformée")
abline(clim99,0,col="red", lty=2)
abline(-clim99,0,col="red", lty=2)
pacf = pacf(xdd, lag.max = 300, main = "PACF série transformée")
abline(clim99,0,col="red", lty=2)
abline(-clim99,0,col="red", lty=2)

which(abs(acf$acf) > clim99)
which(abs(pacf$acf) > clim99)
#... on peut prendre Qmax = 5. 

# On a donc les ordres maximums vraisemblables 
# (p*,d*,q*,P*,D*,Q*,s) = (14,1,3,4,1,5,12)

# La fonction suivante automatise le processus 
# de significativite et validation puis calcul
# des AIC, BIC. 
modelchoice = function(p,q,P,Q,serie,k=24){
  print(paste0("SARIMA (p=",p,",d=",1,",q=",q,",P=",P,",D=",0,",Q=",Q,",s=",12,")"))
  estim = try(arima(serie, order=c(p,0,q), seasonal=list(order = c(P,0,Q), period=12),
                     method="ML",
                     optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"d"=0,"q"=q,"P"=P,"D"=0,"Q"=Q,"s"=12,
                                          "ARsignif"=NA,"MAsignif"=NA,
                                          "ARSAISONsignif"=NA,"MASAISONsignif"=NA,
                                          "ModeleSignif"=NA,
                                          "ModeleValide"=NA,
                                          "Selectionne"=NA,
                                          "AIC"=NA, "BIC"=NA))
  arsignif = if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif = if (q==0) NA else signif(estim)[3,p+q]<=0.05
  ARsaisonsignif = if (P==0) NA else signif(estim)[3,p+q+P]<=0.05
  MAsaisonsignif = if (Q==0) NA else signif(estim)[3,p+q+P+Q]<=0.05
  ModeleSignif   = prod(apply(matrix(c(arsignif,masignif,ARsaisonsignif,MAsaisonsignif)), 1, tf_binaire))
  resnocorr = sum(Qtest(estim$residuals,k)[,3])==k
  aic = AIC(estim)
  bic = BIC(estim)
  return(c("p"=p,"d"=0,"q"=q,"P"=P,"D"=0,"Q"=Q,"s"=12,
           "ARsignif"=arsignif,"MAsignif"=masignif,
           "ARSAISONsignif"=ARsaisonsignif,"MASAISONsignif"=MAsaisonsignif,
           "ModeleSignif"=ModeleSignif,"ModeleValide"=resnocorr,
           "Selectionne"=ModeleSignif*resnocorr,
           "AIC"=aic, "BIC"=bic))
}

# On verifie si des modeles AR sont selectionnes...
pmax = 30 #par extreme securite 
resultAR = data.frame(t(mapply(function(p) modelchoice(p,0,0,0,xdd), c(0:pmax))))
resultAR = resultAR[with(resultAR, order(resultAR$Selectionne, decreasing = TRUE)),]
View(resultAR)
ARselec = resultAR[resultAR$Selectionne == 1,]
ARselec = ARselec[with(ARselec, order(ARselec$AIC, decreasing = FALSE)),]
View(ARselec)
# On selectionne AR(29), AR(25), AR(28), AR(24), AR(21)... dans l'ordre des
# AIC les plus faibles. Les ordres p sont beaucuop trop eleves pour 
# ce type de modelisation. On ecarte donc une modelisation de la serie
# xdd par un AR. D'autant plus que les BIC sont tres grands.

# Des modeles MA ? ...
qmax = 30 #par extreme securite
resultMA = data.frame(t(mapply(function(q) modelchoice(0,q,0,0,xdd), c(0:qmax))))
resultMA = resultMA[with(resultMA, order(resultMA$Selectionne, decreasing = TRUE)),]
View(resultMA)
MAselec = resultMA[resultMA$Selectionne == 1,]
MAselec = MAselec[with(MAselec, order(MAselec$AIC, decreasing = FALSE)),]
View(MAselec)
# On selectionne MA(13), MA(20), MA(27), MA(26)... dans l'ordre des
# AIC les plus faibles. On retient donc le MA(13) puisqu'il minimise 
# a la fois l'AIC et le BIC.

# Des modeles ARMA ? 
pmax = 20
qmax = 20
resultARMA = data.frame(t(mapply(function(p,q) modelchoice(p,q,0,0,xdd),
                         c(apply(matrix(c(0:pmax)),1, function(k) rep(k,(qmax+1)))),
                         c(0:qmax))))

resultARMA = resultARMA[with(resultARMA, order(resultARMA$Selectionne, decreasing = TRUE)),]
View(resultARMA)
ARMAselec = resultARMA[resultARMA$Selectionne == 1,]
ARMAselec = ARMAselec[with(ARMAselec, order(ARMAselec$AIC, decreasing = FALSE)),]
View(ARMAselec)
length(which(is.na(ARMAselec[,16])==0))
# On garde 75 modeles ARMA au total. 
# Les deux modeles ARMA optimaux sont ARMA(0,13) et ARMA(1,13). Ils minimisent 
# chacun un critère. 

# Des modeles SARIMA ? 
pmax = 5
qmax = 5
Pmax = 5
Qmax = 5
resultSARIMA = data.frame(t(mapply(function(p,q,P,Q) modelchoice(p,q,P,Q,xdd),
                  c(apply(matrix(c(0:pmax)),1, function(k) rep(k,(qmax+1)*(Pmax+1)*(Qmax+1)))),
                  c(apply(matrix(c(0:qmax)),1, function(k) rep(k,(Pmax+1)*(Qmax+1)))),
                  c(apply(matrix(c(0:Pmax)),1, function(k) rep(k,(Qmax+1)))),
                  c(0:Qmax))))
resultSARIMA = resultSARIMA[with(resultSARIMA, order(resultSARIMA$Selectionne, decreasing = TRUE)),]
View(resultSARIMA)
SARIMAselec = resultSARIMA[resultSARIMA$Selectionne == 1,]
SARIMAselec = SARIMAselec[with(SARIMAselec, order(SARIMAselec$AIC, decreasing = FALSE)),]
View(SARIMAselec)
length(which(is.na(SARIMAselec[,16])==0))
# On selectionne 182 modeles. Les 5 premiers sont par ordre croissant d'AIC:
# SARIMA(0,0,3,4,0,2,12), SARIMA(4,0,5,4,0,1,12), SARIMA(4,0,5,4,0,2,12)
# SARIMA(4,0,5,0,0,5,12), SARIMA(0,0,3,2,0,4,12)

# On concatene tous les modeles retenus.
resultats_retenus            = rbind(ARselec,MAselec,ARMAselec,SARIMAselec)
resultats_retenus            = resultats_retenus[with(resultats_retenus,
                                           order(resultats_retenus$AIC, decreasing = FALSE)),]
resultats_retenus["AIC+BIC"] = resultats_retenus$AIC + resultats_retenus$BIC
View(resultats_retenus)
length(which(is.na(resultats_retenus[,16])==0))
# Au total, ce sont 266 modeles qui sont retenus. Bien sur,
# il reste maintenant a selectionner le ou les meilleurs du point 
# de la minimisation des criteres.

# Pour le latex :
print(xtable(ARselec,     type = "latex", tabular.environment = "longtable"),
      file = "ARselec.tex")
print(xtable(MAselec,     type = "latex", tabular.environment = "longtable"),
      file = "MAselec.tex")
print(xtable(ARMAselec,   type = "latex", tabular.environment = "longtable"),
      file = "ARMAselec.tex")
print(xtable(SARIMAselec, type = "latex", tabular.environment = "longtable"),
      file = "SARIMAselec.tex")

###############################################
########## Selection d'un modèle ARMA #########
###############################################

# Conclusion : il nous reste a comparer les dix premiers modeles suivants : 
# celui maximisant l'AIC : SARIMA(0,0,3,4,0,2,12)
# celui maximisation le BIC : SARIMA(2,0,1,0,0,1,12)
# 8 autres modèles ayant une somme de criteres AIC + BIC parmi les plus petites :
# SARIMA(0,0,1,4,0,2,12) ; SARIMA(0,0,3,2,0,4,12)
# SARIMA(0,0,1,2,0,4,12) ; SARIMA(0,0,3,0,0,5,12)
# SARIMA(0,0,3,5,0,0,12) ; SARIMA(2,0,1,0,0,5,12)
# SARIMA(0,0,1,0,0,5,12) ; SARIMA(0,0,3,4,0,0,12)
# Ces sont valides et bien ajustes.
# Les modèles ne sont pas imbriques, on ne peut donc pas faire un LR test 
# (rapport de vraisemblance).Il nous reste un critere de prevision pour choisir 
# entre les 10 modeles. On choisit de comparer les erreurs de prediction 
# pour selectionner un unique modele.

T = length(xmlog)
trend = 1:(T-4)
xmlog_tronq = xmlog[trend]
obs = xmlog[(T-3):T]
lt = lm(xmlog[trend] ~ poly(trend,2))

# Ici on rebascule en xmlog, on doit donc prendre en compte le d = 1, D = 1 
# faisant reference a notre differenciation et desaisonnalisation.
sarima1  = arima(xmlog_tronq, order=c(0,1,3), seasonal=list(order = c(4,1,2), period=12), method="ML", optim.control=list(maxit=20000))
sarima2  = arima(xmlog_tronq, order=c(2,1,1), seasonal=list(order = c(0,1,1), period=12), method="ML", optim.control=list(maxit=20000))
sarima3  = arima(xmlog_tronq, order=c(0,1,1), seasonal=list(order = c(4,1,2), period=12), method="ML", optim.control=list(maxit=20000))
sarima4  = arima(xmlog_tronq, order=c(0,1,3), seasonal=list(order = c(2,1,4), period=12), method="ML", optim.control=list(maxit=20000))
sarima5  = arima(xmlog_tronq, order=c(0,1,0), seasonal=list(order = c(2,1,4), period=12), method="ML", optim.control=list(maxit=20000))
sarima6  = arima(xmlog_tronq, order=c(0,1,3), seasonal=list(order = c(0,1,5), period=12), method="ML", optim.control=list(maxit=20000))
sarima7  = arima(xmlog_tronq, order=c(0,1,3), seasonal=list(order = c(5,1,0), period=12), method="ML", optim.control=list(maxit=20000))
sarima8  = arima(xmlog_tronq, order=c(2,1,1), seasonal=list(order = c(0,1,5), period=12), method="ML", optim.control=list(maxit=20000))
sarima9  = arima(xmlog_tronq, order=c(0,1,0), seasonal=list(order = c(0,1,5), period=12), method="ML", optim.control=list(maxit=20000))
sarima10 = arima(xmlog_tronq, order=c(0,1,3), seasonal=list(order = c(4,1,0), period=12), method="ML", optim.control=list(maxit=20000))
modelComp = list(sarima1, sarima2, sarima3, sarima4, sarima5, sarima6, sarima7, sarima8, sarima9, sarima10)

# On introduit la fonction rmserror qui calcule la somme des carres des 
# residus sur l'ensemble de la periode tronquee. 

rmserror = function(model) sqrt(sum(model$residuals^2))

i = 0
predmodel = rep(list(list(),10))
erreur = c("SARIMA(0,0,3,4,0,2,12)minAIC"=0,
           "SARIMA(2,0,1,0,0,1,12)minBIC"=0,
           "SARIMA(0,0,1,4,0,2,12)"=0, "SARIMA(0,0,3,2,0,4,12)"=0,
           "SARIMA(0,0,1,2,0,4,12)"=0, "SARIMA(0,0,3,0,0,5,12)"=0,
           "SARIMA(0,0,3,5,0,0,12)"=0, "SARIMA(2,0,1,0,0,5,12)"=0,
           "SARIMA(0,0,1,0,0,5,12)"=0, "SARIMA(0,0,3,4,0,0,12)"=0)
RMSE = erreur
for (model in modelComp){
  i = i+1
  print(model)
  predmodel[[i]] = as.zoo(ts(predict(model,4)$pred,
                        start     = c(2017,10),
                        end       = c(2018,1),
                        frequency = 12))
  erreur[i]      = sqrt(sum((predmodel[[i]]-obs)^2)/4)
  RMSE[i]        = rmserror(model)
}
# Erreur correspond a la somme des carres des erreurs de prediction 
# sur les 4 derniers mois de la serie xmlog pour chaque modele.
# RMSE regroupe les RMSE des differents modeles.
erreur = data.frame(erreur)
colnames(erreur) = c("Erreur prediction")
RMSE = data.frame(RMSE)
colnames(RMSE) = c("RMSE")
resultsPRED = cbind(erreur,RMSE) # resultsPRED regroupe l'ensemble des erreurs
View(resultsPRED) 
print(xtable(resultsPRED, type = "latex", tabular.environment="longtable"), file = "resultsPRED.tex")

# On remarque que l'erreur de prediction est minimisee par 
# le SARIMA(0,0,1,2,0,4,12). Le RMSE est minimise
# par le SARIMA(0,0,3,4,0,2,12), i.e le SARIMA minimisant
# l'AIC. On remarque que l'erreur de prediction a horizon 4 periodes
# est plus petite pour le SARIMA minimisant le BIC. Cependant
# le RMSE de ce meme modele est plus faible.
# On a don d'un cote un modele minimisant l'AIC et 
# qui fit mieux la serie sur les donnees observees (en fait,
# parmi les modeles retenus, il fit LE mieux la serie xmlog) mais predit moins bien.
# De l'autre, un autre modele minimisant le BIC et 
# predit mieux mais qui fit moins bien la serie sur les donnees observees.
# Enfin, il existe 2 autres modeles ne minimisant aucun des 2 criteres 
# mais qui pourtant predisent mieux a horizon 4 periodes.


# On affiches les erreurs de prediction.
pred = plot(dates[(T-3):T],
             obs, 
             type = 'o',
             xlab = "Temps",
             ylim = c(1.7,2),
             ylab = 'Logarithme de la fréquentation des passagers (millions)',
             col  = 'black',
             main = 'Prévision',
             xaxt = "n")
axis(side=1, at=dates[(T-3):T], labels=format(dates[(T-3):T], '%b-%y'))
lines(predmodel[[1]],                  type = "o", col = "red",      lty = 1, lwd = 2)
lines(dates[(T-3):T], predmodel[[2]],  type = "o", col = "blue",     lty = 1, lwd = 2)
lines(dates[(T-3):T], predmodel[[3]],  type = "o", col = "green",    lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[4]],  type = "o", col = "gray",     lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[5]],  type = "o", col = "orange",   lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[6]],  type = "o", col = "darkblue", lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[7]],  type = "o", col = "cyan",     lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[8]],  type = "o", col = "brown",    lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[9]],  type = "o", col = "yellow",   lty = 1, lwd = 1)
lines(dates[(T-3):T], predmodel[[10]], type = "o", col = "pink",     lty = 1, lwd = 1)

legend(y         = 2.0,
       x         = dates[[T-1]],
       bty       = "n",
       legend    = c("Valeurs observées","SARIMA(0,0,3,4,0,2,12) minimisateur AIC",
                    "SARIMA(2,0,1,0,0,1,12) minimisateur BIC",
                    "SARIMA(0,0,1,4,0,2,12)", "SARIMA(0,0,3,2,0,4,12)",
                    "SARIMA(0,0,1,2,0,4,12)", "SARIMA(0,0,3,0,0,5,12)",
                    "SARIMA(0,0,3,5,0,0,12)" , "SARIMA(2,0,1,0,0,5,12)",
                    "SARIMA(0,0,1,0,0,5,12)", "SARIMA(0,0,3,4,0,0,12)"),
       
       col       = c("black","red", "blue", "green", "gray", "orange", "darkblue", "cyan",
                    "brown", "yellow", "pink"),
       
       text.col  = c("black","red", "blue", "green", "gray", "orange", "darkblue", "cyan",
                    "brown", "yellow", "pink"),
       lty       = c(1,1,1,1,1,1,1,1,1,1),
       pch       = c(NA,1,1,1,1,1,1,1,1,1,1),
       y.intersp = 0.5,
       ncol      = 1)

###############################################
########## PARTIE 3 : Predictions #############
###############################################

# Modele retenu

model = arima(xmlog, order=c(0,1,3), seasonal=list(order = c(4,1,2), period=12), method="ML", optim.control=list(maxit=20000))

# On realise une prevision sur 10 mois
lag        = 15
lastobs    = xmlog[T]
forecasts  = forecast(model, h=lag, level=95)
prediction = as.zoo(ts(c(lastobs,forecasts$mean), frequency = 12, start=c(2018,1)))
up         = as.zoo(ts(c(lastobs,forecasts$upper), frequency = 12, start=c(2018,1)))
down       = as.zoo(ts(c(lastobs,forecasts$lower), frequency = 12, start=c(2018,1)))
datesPRED  = seq(as.POSIXct("2018-02-01"), by = "month", length.out = lag+1)

par(mfrow=c(1,1))
plot(c(dates[(T-48):T],datesPRED), 
     c(xmlog[(T-48):T],as.zoo(ts(rep(NA,length(datesPRED)), start = c(2018,2),frequency = 12))), 
     type = 'l', 
     xlab = '', 
     ylab = 'Logarithme de la fréquentation des passagers (millions)',
     col  = "black",
     ylim = c(1.5,2.25),           
     xaxt = "n",
     main = paste0("Prévision à horizon ",lag," mois"))

axis(side=1, at=c(dates[(T-48):T],datesPRED), labels=format(c(dates[(T-48):T],datesPRED), '%b-%y'), las=2)

lines(dates[(T-48):T], fitted(model)[(T-48):T] , type = "o", col = "red",   lty = 1, lwd = 1)
lines(datesPRED, up                            , type = "o", col = "green", lty = 2, lwd = 1)
lines(datesPRED, down                          , type = "o", col = "green", lty = 2, lwd = 1)
lines(datesPRED, prediction                    , type = "o", col = "blue",  lty = 1, lwd = 1)

legend("toplef", 
       bty       = "n",
       legend    = c("Série observée", "Modèle SARIMA", "Prédiction", "Borne supérieures et inférieures à 95%"),
       col       = c("black","red", "blue", "green"),
       text.col  = c("black","red", "blue", "green"),
       lty       = c(1,1,1,2),
       pch       = c(NA,1,1,1,1),
       ncol      = 1,
       y.intersp = 0.5)

# On verifie bien que l'intervalle de confiance grandit selon que la 
# prevision se fait a horizon plus large.

plot(datesPRED, 
     up-down,
     type = 'o', 
     xlab = '', 
     ylab = "Largeur de l'intervalle de confiance a 95% (millions de passager)",
     col  = "red",
     xaxt = "n",
     main = "Largeur de l'intervalle de confiance a 95%")

axis(side=1, at=datesPRED, labels=format(datesPRED, '%b-%y'), las=2)



