#===================================//=========================================#
# Mail: <<< joseclaudio.faria@gmail.com >>>
#       <<< phgrosjean@sciviews.org >>>
#
# Mais: <<< http://zoonek2.free.fr/UNIX/48_R/all.html >>>
#===================================//=========================================#

# help.start()
# Inicializa a interface em html do help on line do R (usando o browser
# disponível no computador)

x <- rnorm(50)
x
y <- rnorm(50)
y
# Gera dois vetores de coordenadas x e y de números pseudo-aleatórios
# e inspeciona os valores gerados:

x11(w=4, h=4); bringToTop(s=T);
plot(x, y)
# Coloca os pontos em um gráfico.
# Note que a janela gráfica abrirá automaticamente:

ls()
# Verifica os objetos existentes na área de trabalho:

rm(x, y)
# Remove objetos que não são mais necessários:

ls()
# Verifica a remoção

x <- 1:20
# Cria um vetor com uma seqüência de números de 1 a 20:

w <- 1 + sqrt(x)/2
# Um vetor de pesos com os desvios padrões de cada observação:

dummy <- data.frame(x=x, y=x + rnorm(x)*w)
dummy
# Monta um ‘data-frame’ de 2 colunas, x e y, e inspecionando o objeto

fm <- lm(y ~ x, data=dummy)
summary(fm)
# Ajusta uma regressão linear simples de y em x e examinando os resultados:

fm1 <- lm(y ~ x, data=dummy, weight=1/w^2)
summary(fm1)
# Como nós sabemos os pesos podemos fazer uma regressão ponderada:

lrf <- with(dummy, lowess(x, y))
# Faz uma regressão local não-paramétrica:

with(dummy, plot(x, y))
# Plota os pontos

lines(x, lrf$y)
# Adiciona a linha de regressãoo local:

abline(0, 1, lty=3)
# e a linha de regressão verdadeira (intercepto 0 e inclinação 1):

abline(coef(fm))
# a linha da regressão sem ponderação:

abline(coef(fm1), col="red")
# e a linha de regressão ponderada:

plot(fitted(fm), resid(fm),
  xlab="Fitted values", ylab="Residuals",
  main="Residuals vs Fitted")
# O gráfico diagnóstico padrão para checar homocedasticidade

qqnorm(resid(fm), main="Residuals Rankit Plot")
# A normal scores plot to check for skewness, kurtosis and outliers. (Não muito útil aqui)

rm(w, fm, fm1, lrf, x, dummy)
# Remove os objetos)

graphics.off()
# Fecha o gráfico

# q()
# Sai do R. Você será se deseja salvar o "workspace", e para uma sessão exploratória
# como esta, você provavelemente não irá querer salvar.
