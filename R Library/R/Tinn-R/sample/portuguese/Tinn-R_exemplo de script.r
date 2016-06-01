#===================================//=========================================#
# Mail: <<< joseclaudio.faria@gmail.com >>>
#       <<< phgrosjean@sciviews.org >>>
#
# Mais: <<< http://zoonek2.free.fr/UNIX/48_R/all.html >>>
#===================================//=========================================#

# help.start()
# Inicializa a interface em html do help on line do R (usando o browser
# dispon�vel no computador)

x <- rnorm(50)
x
y <- rnorm(50)
y
# Gera dois vetores de coordenadas x e y de n�meros pseudo-aleat�rios
# e inspeciona os valores gerados:

x11(w=4, h=4); bringToTop(s=T);
plot(x, y)
# Coloca os pontos em um gr�fico.
# Note que a janela gr�fica abrir� automaticamente:

ls()
# Verifica os objetos existentes na �rea de trabalho:

rm(x, y)
# Remove objetos que n�o s�o mais necess�rios:

ls()
# Verifica a remo��o

x <- 1:20
# Cria um vetor com uma seq��ncia de n�meros de 1 a 20:

w <- 1 + sqrt(x)/2
# Um vetor de pesos com os desvios padr�es de cada observa��o:

dummy <- data.frame(x=x, y=x + rnorm(x)*w)
dummy
# Monta um �data-frame� de 2 colunas, x e y, e inspecionando o objeto

fm <- lm(y ~ x, data=dummy)
summary(fm)
# Ajusta uma regress�o linear simples de y em x e examinando os resultados:

fm1 <- lm(y ~ x, data=dummy, weight=1/w^2)
summary(fm1)
# Como n�s sabemos os pesos podemos fazer uma regress�o ponderada:

lrf <- with(dummy, lowess(x, y))
# Faz uma regress�o local n�o-param�trica:

with(dummy, plot(x, y))
# Plota os pontos

lines(x, lrf$y)
# Adiciona a linha de regress�oo local:

abline(0, 1, lty=3)
# e a linha de regress�o verdadeira (intercepto 0 e inclina��o 1):

abline(coef(fm))
# a linha da regress�o sem pondera��o:

abline(coef(fm1), col="red")
# e a linha de regress�o ponderada:

plot(fitted(fm), resid(fm),
  xlab="Fitted values", ylab="Residuals",
  main="Residuals vs Fitted")
# O gr�fico diagn�stico padr�o para checar homocedasticidade

qqnorm(resid(fm), main="Residuals Rankit Plot")
# A normal scores plot to check for skewness, kurtosis and outliers. (N�o muito �til aqui)

rm(w, fm, fm1, lrf, x, dummy)
# Remove os objetos)

graphics.off()
# Fecha o gr�fico

# q()
# Sai do R. Voc� ser� se deseja salvar o "workspace", e para uma sess�o explorat�ria
# como esta, voc� provavelemente n�o ir� querer salvar.
