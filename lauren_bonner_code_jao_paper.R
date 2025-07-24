
# ------------------------ CODE from Lauren Bonner (email 7/22/25)

#create variance-covariance matrix
#subgroup_cpep is a subset of the dataset with the analytes we included (each column is a different analyte)
cov<-var(subgroup_cpep)

#BIC to select glasso tuning parameter
rhoarray<-exp(seq(log(0.01), log(0.5), length = 20))
BIC = rep(0, length(rhoarray))
for (rh in 1:length(rhoarray)) {
  fit.gl1 = glasso(cov, rho = rhoarray[rh])
  BIC[rh] = extendedBIC(gamma = 0, omegahat = fit.gl1$wi, 
                        S = cov, n = nrow(subgroup_cpep))
}
rho = rhoarray[which.min(BIC)]
rho

#fit glasso with 'best' tuning parameter
lasso_heu_homa_new<-glasso(cov,rho=rho, trace=TRUE)
#this contains the glasso output (list with multiple components including the covariance and inverse covariance matrix)

#Network Plots
#this next section was just getting names/colors for nodes
subgroup_cpep_names<-colnames(subgroup_cpep)
color<-labels$Color2
names(color)<-labels$omics
color_subgroup<-color[subgroup_cpep_names]

#qgraph uses the glasso output as input and generates the partial correlations - here we specified a threshold of 0.2. The output has an edgelist with nodes and weights (partial correlations)- I think I manually checked the weights with the output from the inverse covariance matrix, just to double check. 
graph<-qgraph(lasso_heu_homa_new, details=T, nodeNames=subgroup_cpep_names, color=color_subgroup, threshold=0.2,label.scale=F, layout="spring")

#this part just creates a data frame of the edgelist (I use this data frame later to get labels)
data<-data.frame(cbind(graph$Edgelist$from,graph$Edgelist$to,graph$Edgelist$weight))
#and this part puts it back in graph form (I can't remember all the details here, but I think we needed an "igraph" instead of "qgraph" object). This may all be irrelevant for you, depending on how you're plotting!
network.for.plot <- graph.data.frame(data, directed = FALSE)

# ------------------------ KAT ADJUSTING THE CODE from Lauren Bonner (email 7/22/25)

#create variance-covariance matrix
#subgroup_cpep is a subset of the dataset with the analytes we included (each column is a different analyte)
cov<-var(dat)

cov <- dat |>
  scale() |>
  as_tibble() |>
  cov(use="complete.obs") 

#BIC to select glasso tuning parameter
# https://cran.r-project.org/web/packages/iDINGO/iDINGO.pdf
rhoarray<-exp(seq(log(0.01), log(0.5), length = 20))
BIC = rep(0, length(rhoarray))
for (rh in 1:length(rhoarray)) {
  fit.gl1 = glasso(cov, rho = rhoarray[rh])
  BIC[rh] = iDINGO::extendedBIC(gamma = 0, omegahat = fit.gl1$wi, 
                        S = cov, n = nrow(dat))
}
rho = rhoarray[which.min(BIC)]
rho

ggplot(data=data.frame(x=rhoarray, y=BIC), aes(x=x, y=y)) +
  geom_point()

#fit glasso with 'best' tuning parameter
lasso_heu_homa_new<-glasso(cov,rho=rho, trace=TRUE)
#this contains the glasso output (list with multiple components including the covariance and inverse covariance matrix)

#Network Plots
#this next section was just getting names/colors for nodes
subgroup_cpep_names<-colnames(dat)
color<-c("blue", "red", "orange", "green")
names(color)<-c("one", "two", "three", "four")
color_subgroup<-color

#qgraph uses the glasso output as input and generates the partial correlations - here we specified a threshold of 0.2.
# The output has an edgelist with nodes and weights (partial correlations)- I think I manually checked the weights with the output from the inverse covariance matrix, 
# just to double check. 
graph2 <- qgraph::qgraph(lasso_heu_homa_new, details=T, nodeNames=subgroup_cpep_names
                        , color=color_subgroup, threshold=0,label.scale=F, layout="spring")
graph2
lasso_heu_homa_new$wi
#this part just creates a data frame of the edgelist (I use this data frame later to get labels)
data <- data.frame(cbind(graph$Edgelist$from,graph$Edgelist$to,graph$Edgelist$weight))

#and this part puts it back in graph form (I can't remember all the details here, but I think we needed an "igraph" instead of "qgraph" object).
# This may all be irrelevant for you, depending on how you're plotting!
network.for.plot <- graph.data.frame(data, directed = FALSE)

n_vars <- length(all_vars)
inv_sig <- lasso_heu_homa_new$wi
variance_diag <- diag(diag(inv_sig), n_vars, n_vars)

# https://online.stat.psu.edu/stat505/book/export/html/638
# https://kongres2025.stat.gov.pl/Content/Docs/Prezentacje/Sess3/01%20-%20mbogdan.pdf
# (slide 37)
pcorrs <- -sqrt(solve(variance_diag)) %*% inv_sig %*% sqrt(solve(variance_diag))
check2 <- pcorrs |>
  as_tibble() |>
  rowid_to_column() |>
  mutate(var1 = paste0("V", rowid)) |>
  select(var1, everything(), -rowid) |>
  pivot_longer(cols=-var1, names_to="var2", values_to="value") |>
  filter(var1 < var2 & value != 0)

# do get different network depending on whether standardize . . .


# --------------------- CODE FROM iDINGO PACKAGE EXAMPLE
# https://cran.r-project.org/web/packages/iDINGO/iDINGO.pdf


library(glasso)
data(gbm)
x = gbm[,1]
Y = gbm[,-1]

# Kat adding: checking if standardized. Y matrix is standardized.
colMeans(Y)
var(Y)
mean(x) # x not used in the code, so can ignore
# Estimating inverse covariance matrix using GLasso #
S = cov(Y)
rhoarray = exp(seq(log(0.001),log(1),length=100))
BIC = rep(0,length(rhoarray))
for (rh in 1:length(rhoarray)) {
  fit.gl1 = glasso(S,rho=rhoarray[rh])
  BIC[rh] = extendedBIC(gamma=0,omegahat=fit.gl1$wi,S=S,n=nrow(Y))
}
rho = rhoarray[which.min(BIC)]
fit.gl2 = glasso(S,rho=rho)
Omega = fit.gl2$wi

# --------------------- HARD CODE BIC TO CHECK
cov <- dat |>
  scale() |>
  as_tibble() |>
  cov() 

rhoarray<-exp(seq(log(0.01), log(0.5), length = 20))

fit_gl1 <- glasso(cov, rho = rhoarray[1], nobs = nrow(dat))
bic1_idingo <- iDINGO::extendedBIC(gamma = 0, omegahat = fit_gl1$wi, 
                                   S = cov, n = nrow(dat)) 
bic1_idingo

# compute via formula: -2*log likelihood PLUS
##                        |E|*log(n) PLUS
#                        4*|E|*gamma*log(p)
# where gamma is a parameter of the *extended* BIC
# and larger values penalize larger graphs. If gamma is 0
# then formula reduces to the classic BIC measure.
# and |E| is the cardinality of the edge set (i.e. the number of edges in the graph)
# https://en.wikipedia.org/wiki/Cardinality
# (this formula makes sense/is consistent with the generic BIC formula:
# https://en.wikipedia.org/wiki/Bayesian_information_criterion)
# BIC = k*ln(n) - 2*ln({\widehat {L}}), where k is the number of parameters 
#                                           estimated by the model

num_edges <- fit_gl1$wi |>
  as_tibble() |>
  rowid_to_column(var="var1") |>
  mutate(var1 = paste0("V",var1)) |>
  pivot_longer(cols=-var1, names_to="var2", values_to="value") |>
  filter(var1 < var2) |>
  as_tbl_graph() |>
  activate(edges) |>
  filter(value != 0) |>
  ecount() # count number of edges (equivalently, could use the `gsize` function)

n <- nrow(dat)
bic1_formula <- (-2*fit_gl1$loglik) + (num_edges*log(n))
bic1_formula

# hmmm got very different value!

# Gao Statistica Sinica 2012 paper
# https://en.wikipedia.org/wiki/Determinant
# https://en.wikipedia.org/wiki/Trace_(linear_algebra)
bic1_formula_b <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
                      (log(n)*num_edges)
bic1_formula_b

# same Gao paper
# sounds like EBIC is unnecessary here; for cases where p increases as n increases 
p <- ncol(cov)
ebic <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
  ((log(n)+4*log(p))*num_edges)

ebic

# well, at least bic1_formula and bic1_formula_b are super close (down to .001)

#BIC to select glasso tuning parameter
# https://cran.r-project.org/web/packages/iDINGO/iDINGO.pdf
#outc <- 1
#grp <- "Control"

# n <- nrow(dat)
# p <- ncol(cov)

rhoarray<-exp(seq(log(0.01), log(1), length = 100))

# BIC1 <- rep(0, length(rhoarray))
# BIC2 <- rep(0, length(rhoarray))
# BIC3 <- rep(0, length(rhoarray))
# EBIC <- rep(0, length(rhoarray))
# EBIC2 <- rep(0, length(rhoarray))

rho_data <- data.frame(matrix(NA, nrow = length(rhoarray), ncol = 5))
colnames(rho_data) <- c("BIC_iDINGO", "BIC_standard", "BIC_standard2"
                       , "EBIC_logn", "EBIC_shutta")
rho_data_all <- data.frame()

for (grp in c("Control", "INSTI", "PrEP")) {
  for (outc in 1:2){
    all_vars <- c(biomarkers, placenta_outcomes[outc])
    
    # (1.) compute covariance matrix
    if (is.na(grp)){
      dat <- paco |>
        select(all_of(all_vars)) |>
        # Need to do case-wise deletion (remove completely if missing one value)
        # because need to give a single sample size in glasso function
        na.omit()
    } else {
      dat <- paco |>
        filter(treatment == grp) |>
        select(all_of(all_vars)) |>
        # Need to do case-wise deletion (remove completely if missing one value)
        # because need to give a single sample size in glasso function
        na.omit()
    }  
    
    cov <- dat |>
      scale() |>
      as_tibble() |>
      cov() 
    
    n <- nrow(dat)
    p <- ncol(cov)
    
    for (rh in 1:length(rhoarray)) {
      
      fit_gl1 <- glasso(cov, rho = rhoarray[rh], nobs = nrow(dat))
      
      # when gamma is 0, the equation corresponds to the usual BIC
      rho_data$BIC_iDINGO[rh] = iDINGO::extendedBIC(gamma = log(n), omegahat = fit_gl1$wi, 
                                     S = cov, n = nrow(dat))
      
      # num_edges <- fit_gl1$wi |>
      #   as_tibble() |>
      #   rowid_to_column(var="var1") |>
      #   mutate(var1 = paste0("V",var1)) |>
      #   pivot_longer(cols=-var1, names_to="var2", values_to="value") |>
      #   filter(var1 < var2) |>
      #   as_tbl_graph() |>
      #   activate(edges) |>
      #   filter(value != 0) |>
      #   ecount() # count number of edges (equivalently, could use the `gsize` function)
      
      num_edges <-  intersect(which(fit_gl1$wi != 0)
                              , which(upper.tri(fit_gl1$wi)==TRUE)) |>
        length()
      
      # this equation is confirmed many places
      # including Shutta et al Statistics in Medicine 2022 (Equation 19)
      # "...higher values of [gamma] encourage a sparser graphical model"
      # since we don't need to encourage more sparsity, we can just use regular BIC!
      rho_data$BIC_standard[rh] <- (-2*fit_gl1$loglik) + (num_edges*log(n))
      
      rho_data$BIC_standard2[rh] <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
        (log(n)*num_edges)
      
      # compared to Eq 19 in Shutta paper, in this equation of EIC, gamma is taken
      # to equal log(n) 
      # SOMETIMES get equivalant min rho (but not values) from iDINGO function when plug log(n) as gamma . . . 
      # but not always . . . (think iDINGO wrong since also seems off when gamma = 0)
      rho_data$EBIC_logn[rh] <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
        ((log(n)+4*log(p))*num_edges)
      
      # Eq. 19 version in Shutta paper (use gamma = log(n))
      rho_data$EBIC_shutta[rh] <- (-2*fit_gl1$loglik) + (num_edges*log(n)) + (4*log(n)*num_edges*log(p))
    }
    
    rho_data_all <- rho_data_all |>
      bind_rows(rho_data |>
                  mutate(rho=rhoarray, grp=grp, outc=outc))
  }
}

# rho1 <- rhoarray[which.min(BIC1)]
# rho1
# rho2 <- rhoarray[which.min(BIC2)]
# rho2
# rho3 <- rhoarray[which.min(BIC3)]
# rho3
# rho4 <- rhoarray[which.min(EBIC)]
# rho4
# rho5 <- rhoarray[which.min(EBIC2)]
# rho5

100*3*2
# should have 600 rows - good.
count(rho_data_all, grp, outc)

rho_data_long <- rho_data_all |>
  pivot_longer(cols=-c(rho, grp, outc), names_to="measure", values_to="value")

ggplot(data=rho_data_long, aes(x=rho, y=value, color=measure, fill=measure)) +
  geom_point() +
  geom_point(data=rho_data_long |> 
               group_by(outc, grp, measure) |> 
               arrange(value) |>
               slice(1), size = 3, shape=21, color="black") +
  facet_grid(outc ~ grp)

# GOOD - BIC_standard and BIC_standard2 match exactly
ggplot(data=filter(rho_data_long, measure %in% c("BIC_standard", "BIC_standard2"))
       , aes(x=rho, y=value, color=grp, fill=grp)) +
  geom_point() +
  geom_point(data=filter(rho_data_long, measure %in% c("BIC_standard", "BIC_standard2")) |> 
               group_by(outc, grp, measure) |> 
               arrange(value) |>
               slice(1), size = 3, shape=21, color="black") +
  facet_grid(outc ~ measure)

# hmmm getting different values for the EBICs here, in the control and prep grp
# don't trust the logn version - trust the EBIC version which is a more clear extension of the BIC
ggplot(data=filter(rho_data_long, measure %in% c("EBIC_shutta", "EBIC_logn"))
       , aes(x=rho, y=value, color=measure, fill=measure)) +
  geom_point() +
  geom_point(data=filter(rho_data_long, measure %in% c("EBIC_shutta", "EBIC_logn")) |> 
               group_by(outc, grp, measure) |> 
               arrange(value) |>
               slice(1), size = 3, shape=21, color="black") +
  facet_grid(outc ~ grp)

# Even though BIC values are different, it is finding the minimum in the same place
# (same rho value minimizes EBIC)
# so this is GOOD!
ggplot(data=filter(rho_data_long, measure %in% c("BIC_iDINGO", "EBIC_shutta"))
       , aes(x=rho, y=value, color=measure, fill=measure)) +
  geom_point() +
  geom_point(data=filter(rho_data_long, measure %in% c("BIC_iDINGO", "EBIC_shutta")) |> 
               group_by(outc, grp, measure) |> 
               arrange(value) |>
               slice(1), size = 3, shape=21, color="black") +
  facet_grid(outc ~ grp)

# some of these different (just like EBIC_logn and EBIC_shutta are different)
# but I think EBIC_logn is WRONG. So EBIC_iDINGO could be correct
ggplot(data=filter(rho_data_long, measure %in% c("BIC_iDINGO", "EBIC_logn"))
       , aes(x=rho, y=value, color=measure, fill=measure)) +
  geom_point() +
  geom_point(data=filter(rho_data_long, measure %in% c("BIC_iDINGO", "EBIC_logn")) |> 
               group_by(outc, grp, measure) |> 
               arrange(value) |>
               slice(1), size = 3, shape=21, color="black") +
  facet_grid(outc ~ grp)


# --------------------- different values of gamma
# check that iDINGO and my logn equation for EBIC match across different values
# of gamma (since gamma = 0 should be equivalent to BIC and iDINGO doesn't 
# match my BIC equation, try to find out where starts to diverge . . .)


# --------------------- SEE IF PARTIAL CORR FORUMULA CAN BE 
# ----------------------  EQUIVALENTLY COMPUTED WITH COV MATRIX . . .

glasso_fit <- glasso(s = covmat
                     # set to NO regularization to start
                     # later use cross-validation to choose the reg. param.?
                     # this returns all nonzero covariance params because
                     # no regularization . . .
                     , rho = 0.01
                     , nobs = nrow(dat))

n_vars <- length(all_vars)
inv_sig <- glasso_fit$wi
precision_diag <- diag(diag(inv_sig), n_vars, n_vars)

# https://online.stat.psu.edu/stat505/book/export/html/638
# https://kongres2025.stat.gov.pl/Content/Docs/Prezentacje/Sess3/01%20-%20mbogdan.pdf
# (slide 37)
pcorrs <- -sqrt(solve(precision_diag)) %*% inv_sig %*% sqrt(solve(precision_diag))
print(pcorrs)

#NOPE NOT EQUAL! MAYBE BECAUSE the numerator is covariance i,j conditional on all
# other Xs but the denominator is wrong -- it's the variance conditional on 
# all other Xs INCLUDING the other (i or j)
var_diag <- diag(diag(glasso_fit$w), n_vars, n_vars)
pcorrs2 <- sqrt(solve(var_diag)) %*% glasso_fit$w %*% sqrt(solve(var_diag))
pcorrs
pcorrs2

