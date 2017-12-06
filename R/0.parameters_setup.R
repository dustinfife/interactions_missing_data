		##### set up demonstration simulation parameters
demo.n = 200
demo.mean = c(0, 3, 3)
demo.sd = c(3, .4, .4)
demo.cor = .3
demo.cov = demo.cor*prod(demo.sd)
demo.slopes = c(.2, .1)
demo.iterations=1000
demo.int =-.5

demo.cov


#60% missing, c = -.5, b = .1, mu.z = 3, mu.x = 0 (standardized scale)