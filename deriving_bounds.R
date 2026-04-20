library(causaloptim)

## setting I
cmpdag <- graph_from_literal(X -+ M:Y, M -+ Y, Ur -+ Y:M) |> initialize_graph()
cmpresp <- create_response_function(cmpdag)
reporig <- create_causalmodel(cmpdag, prob.form =  list(cond = "X", out = c("M", "Y")))

# remove the invalid response patterns for Y
keepYresps <- sapply(cmpresp$Y$matrices, \(mm) {
  !any(mm$Y == 1 & mm$M == 1)
})

cmpresp$Y$values <- cmpresp$Y$values[which(keepYresps)]
cmpresp$Y$matrices <- cmpresp$Y$matrices[which(keepYresps)]
cmpresp$Y$index <- seq_along(which(keepYresps)) - 1

cmpmod <- create_causalmodel(graph = NULL, respvars = cmpresp, prob.form = list(cond = "X", out = c("M", "Y")), 
                             p.vals = subset(expand.grid(X = 0:1, M = 0:1, Y = 0:1), !(M == 1 & Y == 1)))


## setting II
cmpdag_SepEff <- graph_from_literal(X -+ Xm, X -+ Xy, Xm -+ M, Xy -+ Y, M -+ Y, Ur -+ M, Ur -+ Y) |> initialize_graph()

V(cmpdag_SepEff)$latent <- c(0, 1, 1, 0, 0, 1)

p.vals <- expand.grid(X = 0:1, M = 0:1, Y = 0:1)

cmpresp_SepEff <- create_response_function(cmpdag_SepEff)
cmpresp_SepEff$Xy$index <- cmpresp_SepEff$Xy$index[1]
cmpresp_SepEff$Xy$values <- cmpresp_SepEff$Xy$values[3]
cmpresp_SepEff$Xm$index <- cmpresp_SepEff$Xm$index[1]
cmpresp_SepEff$Xm$values <- cmpresp_SepEff$Xm$values[3]

# remove the invalid response patterns for Y
keepYresps_SepEff <- sapply(cmpresp_SepEff$Y$matrices, \(mm) {
  !any(mm$Y == 1 & mm$M == 1)
})

cmpresp_SepEff$Y$values <- cmpresp_SepEff$Y$values[which(keepYresps_SepEff)]
cmpresp_SepEff$Y$matrices <- cmpresp_SepEff$Y$matrices[which(keepYresps_SepEff)]
cmpresp_SepEff$Y$index <- seq_along(which(keepYresps_SepEff)) - 1

reporig_SepEff <- create_causalmodel(graph = NULL, respvars = cmpresp_SepEff, prob.form =  list(cond = c("X"), out = c("M", "Y")), right.vars = c("M", "Y"), 
                                     p.vals = subset(expand.grid(X = 0:1, M = 0:1, Y = 0:1), !(M == 1 & Y == 1)))



### bounds under competing events constraint
#CDE(0)
CDE_0_bnds <- create_linearcausalproblem(cmpmod, "p{Y(M = 0, X = 1) = 1} - p{Y(M = 0, X = 0) = 1}") |> optimize_effect_2()

#NDE(0)
NDE_0_bnds <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 0), X = 1) = 1} - p{Y(M(X = 0), X = 0) = 1}") |> optimize_effect_2()
#NDE(1)
NDE_1_bnds <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 1), X = 0) = 1}") |> optimize_effect_2()
#NIE(0)
NIE_0_bnds <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 0) = 1} - p{Y(M(X = 0), X = 0) = 1}") |> optimize_effect_2()
#NIE(1)
NIE_1_bnds <- create_linearcausalproblem(cmpmod, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 0), X = 1) = 1}") |> optimize_effect_2()

#SDE(0)
SDE_0_bnds <- create_linearcausalproblem(reporig_SepEff, "p{Y(Xy = 1, M(Xm = 0)) = 1} - p{Y(Xy = 0, M(Xm = 0)) = 1}") |> optimize_effect_2()
#SDE(1)
SDE_1_bnds <- create_linearcausalproblem(reporig_SepEff, "p{Y(Xy = 1, M(Xm = 1)) = 1} - p{Y(Xy = 0, M(Xm = 1)) = 1}") |> optimize_effect_2()
#SIE(0)
SIE_0_bnds <- create_linearcausalproblem(reporig_SepEff, "p{Y(Xy = 0, M(Xm = 1)) = 1} - p{Y(Xy = 0, M(Xm = 0)) = 1}") |> optimize_effect_2()
#SIE(1)
SIE_1_bnds <- create_linearcausalproblem(reporig_SepEff, "p{Y(Xy = 1, M(Xm = 1)) = 1} - p{Y(Xy = 1, M(Xm = 0)) = 1}") |> optimize_effect_2()


## bounds without competing events constraint
#CDE(0)
CDE_0_NoConstraint <- analyze_graph(cmpdag, "p{Y(M = 0, X = 1) = 1} - p{Y(M = 0, X = 0) = 1}", constraint = NULL)
CDE_0_NoConstraint_bnds <- optimize_effect_2(CDE_0_NoConstraint)

#NDE(0)/SDE(0) 
NDE_SDE_0_NoConstraint <- analyze_graph(cmpdag,effect = "p{Y(M(X = 0), X = 1) = 1} - p{Y(M(X = 0), X = 0) = 1}", constraint = NULL)
NDE_SDE_0_NoConstraint_bnds <- optimize_effect_2(NDE_SDE_0_NoConstraint)
#NDE(1)/SDE(1) 
NDE_SDE_1_NoConstraint <- analyze_graph(cmpdag, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 1), X = 0) = 1}", constraint = NULL)
NDE_SDE_1_NoConstraint_bnds <- optimize_effect_2(NDE_SDE_1_NoConstraint)
#NIE(0)/SIE(0)
NIE_SIE_0_NoConstraint <- analyze_graph(cmpdag, "p{Y(M(X = 1), X = 0) = 1} - p{Y(M(X = 0), X = 0) = 1}", constraint = NULL)
NIE_SIE_0_NoConstraint_bnds <- optimize_effect_2(NIE_SIE_0_NoConstraint)
#NIE(1)/SIE(1) 
NIE_SIE_1_NoConstraint <- analyze_graph(cmpdag, "p{Y(M(X = 1), X = 1) = 1} - p{Y(M(X = 0), X = 1) = 1}", constraint = NULL)
NIE_SIE_1_NoConstraint_bnds <- optimize_effect_2(NIE_SIE_1_NoConstraint)

