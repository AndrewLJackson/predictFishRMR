set.seed(1)

y <- rnorm(100, mean = 100, sd = 2)
x <- rnorm(100, mean = 100, sd = 2)

y2 <- log10(y)
x2 <- log10(x)


m1 <- glm(log10(y)~ log10(x))
m2 <- glm(y2 ~ x2)

x_new <- 100

10^predict(m1, newdata = data.frame(x = x_new), type = "link")
10^predict(m2, newdata = data.frame(x2 = x_new), type = "link")


(10^coef(m1)[1]) * (x_new ^ coef(m1)[2]) * 10^(0.5 * 0.007832^2)


