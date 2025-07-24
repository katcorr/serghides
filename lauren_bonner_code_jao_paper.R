
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
rhoarray<-exp(seq(log(0.01), log(0.5), length = 20))
BIC1 <- rep(0, length(rhoarray))
BIC2 <- rep(0, length(rhoarray))
BIC3 <- rep(0, length(rhoarray))
EBIC <- rep(0, length(rhoarray))
for (rh in 1:length(rhoarray)) {
  fit_gl1 = glasso(cov, rho = rhoarray[rh], nobs = nrow(dat))
  
  BIC1[rh] = iDINGO::extendedBIC(gamma = 0, omegahat = fit_gl1$wi, 
                                S = cov, n = nrow(dat))
  
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
  
  BIC2[rh] <- (-2*fit_gl1$loglik) + (num_edges*log(n))
  
  BIC3[rh] <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
    (log(n)*num_edges)
  
  EBIC[rh] <- (-n*log(det(fit_gl1$wi))) + (n*sum(diag(fit_gl1$wi%*%fit_gl1$w))) +
    ((log(n)+4*log(p))*num_edges)
}

rho1 <- rhoarray[which.min(BIC1)]
rho1
rho2 <- rhoarray[which.min(BIC2)]
rho2
rho3 <- rhoarray[which.min(BIC3)]
rho3
rho4 <- rhoarray[which.min(EBIC)]
rho4

rho_data <- data.frame(rho=rhoarray
                       , BIC1 = BIC1
                       , BIC2 = BIC2
                       , BIC3 = BIC3
                       , EBIC = EBIC) |>
  pivot_longer(cols=-rho, names_to = "measure", values_to = "value")

ggplot(data=rho_data, aes(x=rho, y=value, color=measure)) +
  geom_point() +
  geom_point(data=rho_data |> 
               group_by(measure) |> 
               arrange(value) |>
               slice(1), size = 2, shape="star")

ggplot() +
  geom_point(data=rho_data, aes(x=rho, y=BIC1)) +
  geom_point(data=rho_data |> 
                    arrange(BIC1) |> 
                    slice(1)
             , aes(x=rho, y=BIC1)
             , color="red")

ggplot() +
  geom_point(data=rho_data, aes(x=rho, y=BIC2)) +
  geom_point(data=rho_data |> 
               arrange(BIC2) |> 
               slice(1)
             , aes(x=rho, y=BIC2)
             , color="red")

ggplot() +
  geom_point(data=rho_data, aes(x=rho, y=BIC1), color="darkred") +
  geom_point(data=rho_data |> 
               arrange(BIC1) |> 
               slice(1)
             , aes(x=rho, y=BIC1)
             , color="red") +
  geom_point(data=rho_data, aes(x=rho, y=BIC2), color="darkblue") +
  geom_point(data=rho_data |> 
               arrange(BIC2) |> 
               slice(1)
             , aes(x=rho, y=BIC2)
             , color="blue")
