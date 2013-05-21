source("load.data.R")

data = load.data.02()
fit = itemset.coxpath(data, trace = 2, max.steps = 100, max.time.per.step=150)
print(fit)
plot(fit, xvar = 'stepcoeffs', type='coefficients')
cv.fits = run.cv(data, training.set.ratio = 0.6, num.repeats = 10, trace = 3, min.lambda = 0.5, depth = 3, max.steps = 100)