#---- PEM Survival Model --- #
#Based on the Applied Survival Models Vignette
library(purrr)
suppressMessages(library(tidyverse))
library(survival)
library(rstan)
library(assertthat)
library(corrplot)
library(cgdsr)
suppressMessages(library(dplyr))

library(ggplot2)
require(ggfortify)
theme_set(theme_bw())
###############################################
#Data obtantion
#------Obtain data by the cgdsr package from MSKCC CBioPortal ----# 

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

study_list = getCancerStudies(mycgds)

id_sutdy = getCancerStudies(mycgds)[55,1]
case_list = getCaseLists(mycgds, id_sutdy)[2,1]
clinical_data <-  tbl_df(getClinicalData(mycgds, case_list)) 
clinical_data <- clinical_data %>% tibble::rownames_to_column("sample") 

#inspect dataframe
glimpse(clinical_data)

#separate two cohorts
id_2008_sutdy = getCancerStudies(mycgds)[56,1] #cohort 2008
case_list_2008 = getCaseLists(mycgds, id_2008_sutdy)[2,1]
clinical_data_2008 <-  tbl_df(getClinicalData(mycgds, case_list_2008))
clinical_data_2008 <- clinical_data_2008 %>% tibble::rownames_to_column("sample")
# 
clinical_data <- clinical_data %>% 
  filter(sample %in% setdiff(clinical_data$sample, clinical_data_2008$sample))


####################################################################
#Data Cleaning


#convert to lower case
names(clinical_data) <- tolower(names(clinical_data)) 

#convert missig values
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
clinical_data <- clinical_data %>%
  dplyr::mutate_all(funs(convert_blank_to_na))

#inspect resulting dataframe
glimpse(clinical_data)
clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

clinical_data %>% 
  filter(is.na(os_status) | os_status == "") %>%
  select(os_months, os_status) %>%
  glimpse

clinical_data %>% 
  filter(is.na(os_status) | os_status == "" |os_months < 0 | is.na(os_months)) %>%
  select(os_months, os_status) %>%
  glimpse

#For the moment we will remove these observations from the analysis
clin_data <- tbl_df(clinical_data)
clinical_data <- 
  clin_data %>% 
  filter(!is.na(os_status) & os_status != "" )

assertthat::assert_that(nrow(clinical_data) == (nrow(clin_data) - 44))
remove(clin_data)

clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)


######################################################################
#######--------------  Data Exploration  ----------------#################
#---------   Considering overall survival   ----------------------#


########## Distribution of event times  ######################

clinical_data %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

#KM curve
mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = clinical_data %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM 2008 Cohort')

##################################################################
#------------ Prepare for fitting model--------------------------#

#Center continuos covariates

ggplot(clinical_data, aes(x = age))+
  geom_density()+
  geom_vline(xintercept = median(clinical_data$age))

centered <- function(x){
  x_centered <- x - mean (x)
  return(x_centered)
}

clinical_data <- clinical_data %>%
  dplyr::mutate( age_centered = centered(age)) 

#Impute or delete missing Covariate values

clinical_data %>%
  select(mgmt_status) %>%
  table(exclude = NULL)

clinical_data %>%
  select(idh1_status, g.cimp_methylation) %>%
  table(exclude = NULL)


#remove nas
clinical_data <- clinical_data %>%
  filter(!is.na(mgmt_status))

clinical_data %>%
  select(mgmt_status, g.cimp_methylation) %>%
  table(exclude = NULL)

#create long dataset, we are going to create a dataset for the PEM , taking into only the times where an event occured, it could also work with calendar times (picked at random), or all observed times including censored times, but it has been shown that observed times work better than random times and the dataset of event times is smaller, that is why we choose this method.

clinical_data <- clinical_data %>%
  mutate(s = seq(n())) #add sample id

#obtain the times where an event occured
times <- clinical_data %>% 
  filter(os_status == "DECEASED") %>% select(os_months)%>%
  unique %>% ungroup %>% 
  arrange(os_months) %>% 
  unlist

longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                cut = times, data = (clinical_data %>%
                                                       mutate(deceased = os_status == "DECEASED")))

#generate stan imput data
gen_stan_data <- function(data, formula = as.formula(~1)) {
  
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  
  #create time point id
  data <- data %>%
    group_by(sample) %>%
    mutate(t = seq(n())) 
  
  #covariate matrix
  x <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(x)
  
  if (M > 1){
    if("(Intercept)" %in% colnames(x))
      x <- array(x[,-1], dim = c(nrow(data), M -1))
    M <- ncol(x)
  }
  
  stan_data <- list(
    N = nrow(data),
    S = dplyr::n_distinct(data$sample),
    "T" = length(times),
    s = array(as.numeric(data$s)),
    M = M,
    t = data$t,
    event = data$deceased,
    obs_t = data$os_months,
    x = array(x, dim = c(nrow(data), M))
  )
}

##------ Starting values --------##


###------ Run Stan --------##
nChains <- 4
stanfile <- "pem_survival_model.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

pem_model <-  stan(stanfile,
                            data = gen_stan_data(longdata, '~ age_centered + mgmt_status + g.cimp_methylation'),
                            iter = 1000,
                            cores = min(1, parallel::detectCores()),
                            chains = 4
)


##--Review model convergence--#

#Fit object
print(pem_model)  #Rhat are close to 1?

#Traceplots
rstan::traceplot(pem_model, 'lp__')
rstan::traceplot(pem_model, 'beta')

if(interactive())
  shinystan::launch_shinystan(pem_model)        #Launch shiny stan. There are some divergent 


##----- Posterior predictive simulations -------#

glio_fit <- extract(pem_model, permuted = TRUE)

##--Estimated parameters---#
hist(glio_fit$baseline[,90])
hist(glio_fit$beta[,2])
hist(glio_fit$beta[,3])
hist(glio_fit$hazard[,11])

#How likely is the coefficient be greater than 0

mean(glio_fit$beta[,1] > 0)
mean(glio_fit$beta[,2] > 0)
mean(glio_fit$beta[,3] > 0)

#How well does this model fit the data

mean(glio_fit$hazard[,1])

stan_estimate <- data_frame( haz_mean = apply(glio_fit$hazard, 2, mean),
                             haz_lower_p = apply(glio_fit$hazard, 2, quantile, prob = .1),
                             haz_upper_p = apply(glio_fit$hazard, 2, quantile, prob = .9))

survdata <- longdata %>%
  select(s, deceased, os_months, age_centered, g.cimp_methylation, mgmt_status) %>%
  cbind(stan_estimate) %>%
  group_by(s) %>%
  mutate(surv_mean = exp(-cumsum(os_months*haz_mean)),
         surv_lower_p = exp(-cumsum(os_months*haz_lower_p)),
         surv_upper_p = exp(-cumsum(os_months*haz_upper_p)))


#Plot

survdata %>%
  group_by(s) %>%
  arrange(os_months) %>%
  filter(row_number() == n()) %>%
  filter(g.cimp_methylation == "non-G-CIMP", mgmt_status == "UNMETHYLATED") %>%
  ggplot(aes(x = os_months, y = surv_lower_p)) +
  geom_line()





