#===================================//=========================================#
# Mail: <<< joseclaudio.faria@gmail.com >>>
#       <<< phgrosjean@sciviews.org >>>
#
# Plus: <<< http://zoonek2.free.fr/UNIX/48_R/all.html >>>
#===================================//=========================================#

# help.start()
# D�marre l'interface graphique de l'aide en ligne (en utilisant le navigateur
# Web par d�faut sur votre machine). Vous devriez explorer bri�vement cette aide
# pour vous familiariser avec elle si vous �tre novice en R.
# Minimisez l'aide maintenant pour pouvoir continuer la d�monstration.

x <- rnorm(50)
x
y <- rnorm(50)
y
# G�n�re deux vecteurs de nombre pseudo-al�atoires pour les coordonn�es x- et y-.

x11(w=4, h=4); bringToTop(s=T);
plot(x, y)
# Graphique en nuage de points. Une fen�tre graphique apparait automatiquement.

ls()
# Visualise les objets actuellement dans l'environnement utilisateur de R.

rm(x, y)
# Efface les objets qui ne sont plus n�cessaires.

x <- 1:20
# Cr�e x = (1, 2, . . . , 20).

w <- 1 + sqrt(x)/2
# Un vecteur 'pond�r�' d'�carts types g�n�r� de mani�re pseudo-al�atoire.

dummy <- data.frame(x=x, y=x + rnorm(x)*w)
# Dummy est un 'data frame' de deux variables, x et y.

fm <- lm(y ~ x, data=dummy)
summary(fm)
# Ajuste une r�gression lin�aire simple de y sur x, et inspecte le r�sultat.

fm1 <- lm(y ~ x, data=dummy, weight=1/w^2)
summary(fm1)
# Connaissant les �carts types, nous ajustons maintenant une r�gression pond�r�e.

lrf <- with(dummy, lowess(x, y))
# Ajuste une r�gression locale non param�trique.

with(dummy, plot(x, y))
# Nuage de points.

lines(x, lrf$y)
# Ajoute une courbe repr�sentant la r�gression locale.

abline(0, 1, lty=3)
# La droite de r�gression du mod�le (ordonn�e � l'origine 0, pente 1).

abline(coef(fm))
# La droite de r�gression non pond�r�e.

abline(coef(fm1), col="red")
# La droite de r�gression pond�r�e.

plot(fitted(fm), resid(fm),
  xlab="Fitted values",
  ylab="Residuals",
  main="Residuals vs Fitted")
# Un graphique diagnostic de base des r�sidus pour rechercher une
# h�t�rosc�dasticit�. La voyez-vous?

qqnorm(resid(fm), main="Residuals Rankit Plot")
# Un graphe quantile-quantile pour rechercher une assym�trie, un aplatissement
# ou des outliers dans la distribution des r�sidus. (rien de significatif dans
# le cas pr�sent).

rm(w, fm, fm1, lrf, x, dummy)
# Nettoye les variables � nouveau.

graphics.off()
# Ferme la ou les fen�tres graphiques ouvertes.

# q()
# Quitte R. Il vous sera demand� si vous voulez sauvegarder la session en cours.
# Dans le cadre de cette d�mo, vous ne voudrez probablement rien sauvegarder.