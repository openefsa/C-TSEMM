#===================================//=========================================#
# Mail: <<< joseclaudio.faria@gmail.com >>>
#       <<< phgrosjean@sciviews.org >>>
#
# Plus: <<< http://zoonek2.free.fr/UNIX/48_R/all.html >>>
#===================================//=========================================#

# help.start()
# Démarre l'interface graphique de l'aide en ligne (en utilisant le navigateur
# Web par défaut sur votre machine). Vous devriez explorer brièvement cette aide
# pour vous familiariser avec elle si vous être novice en R.
# Minimisez l'aide maintenant pour pouvoir continuer la démonstration.

x <- rnorm(50)
x
y <- rnorm(50)
y
# Génère deux vecteurs de nombre pseudo-aléatoires pour les coordonnées x- et y-.

x11(w=4, h=4); bringToTop(s=T);
plot(x, y)
# Graphique en nuage de points. Une fenêtre graphique apparait automatiquement.

ls()
# Visualise les objets actuellement dans l'environnement utilisateur de R.

rm(x, y)
# Efface les objets qui ne sont plus nécessaires.

x <- 1:20
# Crée x = (1, 2, . . . , 20).

w <- 1 + sqrt(x)/2
# Un vecteur 'pondéré' d'écarts types généré de manière pseudo-aléatoire.

dummy <- data.frame(x=x, y=x + rnorm(x)*w)
# Dummy est un 'data frame' de deux variables, x et y.

fm <- lm(y ~ x, data=dummy)
summary(fm)
# Ajuste une régression linéaire simple de y sur x, et inspecte le résultat.

fm1 <- lm(y ~ x, data=dummy, weight=1/w^2)
summary(fm1)
# Connaissant les écarts types, nous ajustons maintenant une régression pondérée.

lrf <- with(dummy, lowess(x, y))
# Ajuste une régression locale non paramétrique.

with(dummy, plot(x, y))
# Nuage de points.

lines(x, lrf$y)
# Ajoute une courbe représentant la régression locale.

abline(0, 1, lty=3)
# La droite de régression du modèle (ordonnée à l'origine 0, pente 1).

abline(coef(fm))
# La droite de régression non pondérée.

abline(coef(fm1), col="red")
# La droite de régression pondérée.

plot(fitted(fm), resid(fm),
  xlab="Fitted values",
  ylab="Residuals",
  main="Residuals vs Fitted")
# Un graphique diagnostic de base des résidus pour rechercher une
# hétéroscédasticité. La voyez-vous?

qqnorm(resid(fm), main="Residuals Rankit Plot")
# Un graphe quantile-quantile pour rechercher une assymétrie, un aplatissement
# ou des outliers dans la distribution des résidus. (rien de significatif dans
# le cas présent).

rm(w, fm, fm1, lrf, x, dummy)
# Nettoye les variables à nouveau.

graphics.off()
# Ferme la ou les fenêtres graphiques ouvertes.

# q()
# Quitte R. Il vous sera demandé si vous voulez sauvegarder la session en cours.
# Dans le cadre de cette démo, vous ne voudrez probablement rien sauvegarder.