## simulation study 

#comparing coverage and % the causal null is outside the bounds between bounds under competing events constraint vs. without constraint

n <- 10000
coverage <- numeric(n)
exclude_0 <- numeric(n)

sim_bounds <- function(n, create_linearcausalproblem_output, bound_fun){
  set.seed(12345)
  for (i in 1:n){
    alpha <- rep(1, length(create_linearcausalproblem_output$variables))  
    
    #qs
    q_samp_i <- rdirichlet(k = length(alpha), alpha = 1)
    names(q_samp_i) <- create_linearcausalproblem_output$variables
    
    #true 
    obj_expr <- parse(text = create_linearcausalproblem_output$objective)[[1]]  
    true_i <- eval(obj_expr, envir = as.list(q_samp_i))  
    
    #ps
    p_pop_i <- as.numeric(create_linearcausalproblem_output$R[-1, ] %*% q_samp_i)
    names(p_pop_i) <- create_linearcausalproblem_output$parameters
    
    #bounds
    bounds_i <- do.call(bound_fun, as.list(p_pop_i))
    #number of time the true causal risk difference is outside the bounds  
    coverage[i] <- (true_i >= bounds_i$lower) & (true_i <= bounds_i$upper)
    exclude_0[i] <- bounds_i$upper <= 0 | bounds_i$lower >= 0
  }
  cat("% coverage:", sum(coverage)/n, "% exclude 0:", sum(exclude_0)/n)
}

CDE_0_cmpmod <- create_linearcausalproblem(cmpmod, "p{Y(M = 0, X = 1) = 1} - p{Y(M = 0, X = 0) = 1}") 
#CDE(0) under/without competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = CDE_0_cmpmod, bound_fun = CDE_0_bnds$bounds_function)


NDE_SDE_0_cmpmod <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 0), X = 1) = 1} - p{Y(M(X = 0), X = 0) = 1}")
#NDE(0)/SDE(0) under competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NDE_SDE_0_cmpmod, bound_fun = NDE_0_bnds$bounds_function)
#NDE(0)/SDE(0) without competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NDE_SDE_0_cmpmod, bound_fun = NDE_SDE_0_NoConstraint_bnds$bounds_function)

NDE_SDE_1_cmpmod <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 1), X = 0) = 1}")
#NDE(1)/SDE(1) under competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NDE_SDE_1_cmpmod, bound_fun = NDE_1_bnds$bounds_function)
#NDE(1)/SDE(1) without competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NDE_SDE_1_cmpmod, bound_fun = NDE_SDE_1_NoConstraint_bnds$bounds_function)

NIE_SIE_0_cmpmod <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 0) = 1} - p{Y(M(X = 0), X = 0) = 1}")
#NIE(0)/SIE(0) under competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NIE_SIE_0_cmpmod, bound_fun = NIE_0_bnds$bounds_function)
#NIE(0)/SIE(0) without competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NIE_SIE_0_cmpmod, bound_fun = NIE_SIE_0_NoConstraint_bnds$bounds_function)

NIE_SIE_1_cmpmod <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 0), X = 1) = 1}")
#NIE(1)/SIE(1) under competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NIE_SIE_1_cmpmod, bound_fun = NIE_1_bnds$bounds_function)
#NIE(1)/SIE(1) without competing risk constraint
sim_bounds(n = n, create_linearcausalproblem_output = NIE_SIE_1_cmpmod, bound_fun = NIE_SIE_1_NoConstraint_bnds$bounds_function)


#compare widths between bounds under competing events constraint vs. without constraint

#CDE(0)
#under competing events constraint
set.seed(12345)
simres_CDE_0 <- lapply(1:10000, \(i) {
  do.call(CDE_0_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod, simplex_sampler = function(k) {
            rdirichlet(k, alpha = 1)})))
})

simres_CDE_0 <- do.call(rbind, simres_CDE_0)
simres_CDE_0$width <- simres_CDE_0$upper - simres_CDE_0$lower
summary(simres_CDE_0$width)

#without competing events constraint
set.seed(12345)
simres_CDE_0_NoConstraint <- lapply(1:10000, \(i) {
  do.call(CDE_0_NoConstraint_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod, simplex_sampler = function(k) {
            rdirichlet(k, alpha = 1)})))
})

simres_CDE_0_NoConstraint <- do.call(rbind, simres_CDE_0_NoConstraint)
simres_CDE_0_NoConstraint$width <- simres_CDE_0_NoConstraint$upper - simres_CDE_0_NoConstraint$lower
summary(simres_CDE_0_NoConstraint$width)


#NDE(0)/SDE(0)
#under competing events constraint
set.seed(12345)
simres_NDE_SDE_0 <- lapply(1:10000, \(i) {
  do.call(NDE_0_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod)))
})
simres_NDE_SDE_0 <- do.call(rbind, simres_NDE_SDE_0)
simres_NDE_SDE_0$width <- simres_NDE_SDE_0$upper - simres_NDE_SDE_0$lower
summary(simres_NDE_SDE_0$width)

#without competing events constraint
set.seed(12345)
simres_NDE_SDE_0_NoConstraint<- lapply(1:10000, \(i) {
  do.call(NDE_SDE_0_NoConstraint_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod, simplex_sampler = function(k) {
            rdirichlet(k, alpha = 1)})))
})

simres_NDE_SDE_0_NoConstraint <- do.call(rbind, simres_NDE_SDE_0_NoConstraint)
simres_NDE_SDE_0_NoConstraint$width <- simres_NDE_SDE_0_NoConstraint$upper - simres_NDE_SDE_0_NoConstraint$lower
summary(simres_NDE_SDE_0_NoConstraint$width)


#NIE(0)/SIE(0)
#under competing events constraint
set.seed(12345)
simres_NIE_SIE_0 <- lapply(1:10000, \(i) {
  do.call(NIE_0_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod, simplex_sampler = function(k) {
            rdirichlet(k, alpha = 1)})))
})

simres_NIE_SIE_0 <- do.call(rbind, simres_NIE_SIE_0)
simres_NIE_SIE_0$width <- simres_NIE_SIE_0$upper - simres_NIE_SIE_0$lower
summary(simres_NIE_SIE_0$width)

#without competing events constraint
set.seed(12345)
simres_NIE_SIE_0_NoConstraint <- lapply(1:10000, \(i) {
  do.call(NIE_SIE_0_NoConstraint_bnds$bounds_function, 
          as.list(sample_distribution(cmpmod, simplex_sampler = function(k) {
            rdirichlet(k, alpha = 1)})))
})

simres_NIE_SIE_0_NoConstraint <- do.call(rbind, simres_NIE_SIE_0_NoConstraint)
simres_NIE_SIE_0_NoConstraint$width <- simres_NIE_SIE_0_NoConstraint$upper - simres_NIE_SIE_0_NoConstraint$lower
summary(simres_NIE_SIE_0_NoConstraint$width)


#histograms
par(mfrow=c(3, 2))
hist(simres_CDE_0_NoConstraint$width, main = "Width of CEI(0) Bounds", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5400), breaks = 20)
hist(simres_CDE_0$width, main = "Width of CEI(0) Bounds under Competing Events Constraint", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5300), breaks = 20)

hist(simres_NDE_SDE_0_NoConstraint$width, main = "Width of NDE(0)/SDE(0) Bounds", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5400), breaks = 15)
hist(simres_NDE_SDE_0$width, main = "Width of NDE(0)/SDE(0) Bounds under Competing Events Constraint", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5300), breaks = 15)

hist(simres_NIE_SIE_0_NoConstraint$width, main = "Width of NIE(0)/SIE(0) Bounds", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5400), breaks = 15)
hist(simres_NIE_SIE_0$width, main = "Width of NIE(0)/SIE(0) Bounds under Competing Events Constraint", xlab = "Width", xlim = c(0, 1.6), ylim = c(0, 5300), breaks = 15)
par(mfrow=c(1, 1))


#scatterplots
plot(simres_CDE_0$width, simres_CDE_0_NoConstraint$width, xlim = c(0,2), ylim = c(0,2), ylab = "Width of CDE(0) Bounds", xlab = "Width of CDE(0) Bounds With Constraint")
abline(a=0, b=1, col="black", lty=2)
mtext("A", side=3, line=1, adj=0, cex=1.5, font=2)
WidthComparison_CDE <- recordPlot()

plot(simres_NDE_SDE_0$width, simres_NDE_SDE_0_NoConstraint$width, xlim = c(0,2), ylim = c(0,2), ylab = "Width of NDE(0)/SDE(0) Bounds", xlab = "Width of NDE(0)/SDE(0) Bounds With Constraint")
abline(a=0, b=1, col="black", lty=2)
mtext("B", side=3, line=1, adj=0, cex=1.5, font=2)
WidthComparison_NDE_SDE <- recordPlot()

plot(simres_NIE_SIE_0$width, simres_NIE_SIE_0_NoConstraint$width, xlim = c(0,2), ylim = c(0,2), ylab = "Width of NIE(0)/SIE(0) Bounds", xlab = "Width of NIE(0)/SIE(0) Bounds With Constraint")
abline(a=0, b=1, col="black", lty=2)
mtext("C", side=3, line=1, adj=0, cex=1.5, font=2)
WidthComparison_NIE_SIE <- recordPlot()


