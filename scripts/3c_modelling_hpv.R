# Henipavirus

require(tidyverse)
require(mlr3)

dataset <- read.csv(".../data/hpv_dataset_main.csv") 

# cor(dataset[-c(1:4, 12, 32, 33, 34)]) %>% 
#   as.data.frame(.) %>%
#   rownames_to_column(var = "variable1") %>% 
#   pivot_longer(cols = 2:27, names_to = "variable2", values_to = "correlation") %>% 
#   arrange(desc(correlation)) %>% 
#   filter((correlation > 0.7 | correlation < -0.7) & correlation != 1) %>% 
#   View()

dataset <- dataset[-c(6, 30, 10, 14, 22, 17, 25, 15, 23, 8, 18, 26, 19, 27, 7)]

# Variable preparation
dataset$status <- as.factor(as.character(dataset$status))
dataset$status <- fct_relevel(dataset$status, "1")
dataset$trophic_level <- as.factor(as.character(dataset$trophic_level))
dataset$area.1 <- log10(dataset$area.1)
dataset$order <- as.factor(dataset$order)
dataset$order2 <- as.factor(dataset$order)
dataset$real_status <- as.factor(dataset$real_status)
sum(dataset$w)

dataset %>% 
  group_by(order, status) %>% 
  summarise(sum_w = sum(w))

sum(dataset$w)
require(mlr3learners)

# Defining tasks
#one with pseudo-negatives
task.pn = as_task_classif(id = "pn", reformulate(names(dataset[-c(1, 2, 4, 17)]), names(dataset[17])), data = dataset, positive = "1")
task.pn$set_col_roles("w", roles = "weight")
task.pn$set_col_roles(c("real_status", "order2"), roles = "stratum") # "real status" is used to stratify

#check
task.pn$col_roles #nice!

# Creating learners and defining hyperparameters to tune
require(mlr3tuning)
# mlr_learners$get("classif.xgboost")$param_set #check tuneable hyperparameters for the learner

learner.rf = lrn("classif.ranger",
                 mtry = to_tune(p_int(2, 15)),
                 num.trees = 1000,
                 min.node.size = to_tune(p_int(5, 90)),
                 predict_type = "prob",
                 importance = "permutation")

learner.xgb = lrn("classif.xgboost",
                  eta = to_tune(p_dbl(lower = 0.001, upper = 0.1)),
                  # max_depth = to_tune(p_int(1, 15)),
                  subsample = to_tune(p_dbl(0.5, 0.9)),
                  min_child_weight = to_tune(p_dbl(0.01, 0.1)),
                  predict_type = "prob",
                  colsample_bytree = 0.8,
                  colsample_bylevel = 0.8,                 
                  colsample_bynode = 0.6) #which is 6  #booster = "gblinear", top_k = 0

learner.nn = lrn("classif.nnet",
                 size = to_tune(p_int(lower = 1, upper = 10)),
                 maxit = to_tune(p_int(10, 150)),
                 decay = to_tune(p_dbl(1e-7, 0.1, logscale = TRUE)),
                 predict_type = "prob")

require(mlr3extralearners)
learner.gam = lrn("classif.gam",
                  formula = status ~ order.Chiroptera + order.Rodentia + order.Eulipotyphla + order.Cetartiodactyla + order.Perissodactyla+
                    trophic_level.3 + trophic_level.2 + s(adult_mass_g) + 
                    s(litters_per_year_n) + s(mean_dist) + s(area.1) + s(weaning_age_d)+
                    s(mean.wc2.1_5m_bio_9) + s(mean.wc2.1_5m_bio_18)+
                    s(stdev.wc2.1_5m_bio_19) + s(mean.wc2.1_5m_bio_19) +
                    s(stdev.wc2.1_5m_bio_18) + s(stdev.wc2.1_5m_bio_9),
                  select = TRUE,
                  method = "REML",
                  gamma = to_tune(p_dbl(1, 1)), #to_tune(p_dbl(1, 1))
                  predict_type = "prob")

# Define new performance measure (TSS)
#template code from mlr3book.mlr-org.com Chapter 10
MeasureClassifTSS = R6::R6Class("MeasureClassifTSS",
                                inherit = mlr3::MeasureClassif,
                                public = list(
                                  initialize = function() { # initialize class
                                    super$initialize(
                                      id = "classif.tss", # unique ID
                                      packages = character(), # no package dependencies
                                      properties = character(), # no special properties
                                      predict_type = "response", # measures response prediction
                                      range = c(-1, 1), # results in values between (-1, 1)
                                      minimize = FALSE # larger values are better
                                    )
                                  }
                                ),
                                
                                private = list(
                                  # define score as private method
                                  .score = function(prediction, ...) {
                                    # define tss
                                    tss.fun = function(truth, response) {
                                      conf_mat = table(response, truth)
                                      conf_mat[1, 1]/(conf_mat[1, 1]+conf_mat[2, 1])+conf_mat[2, 2]/(conf_mat[2, 2]+conf_mat[1, 2])-1
                                    }
                                    # call tss function
                                    tss.fun(prediction$truth, prediction$response)
                                  }
                                )
)

# add tss to library
mlr_measures$add("classif.tss", MeasureClassifTSS)


# Variable pre-processing:
require(mlr3pipelines)
require(igraph)

non_pred_selector <- selector_name(c("status", "real_status", "w"))
pred_selector <- selector_invert(non_pred_selector)

# I'm putting all the pre-processing steps into one graph
graph = gunion(list(po("select", selector = non_pred_selector, id = "select_non_predictors") %>>% po("nop"),
                    po("select", selector = pred_selector, id = "select_predictors") %>>% po("encode", method = "one-hot", affect_columns = selector_type("factor"), id = "factor_encoding") %>>%  # 1) factor encoding
                    po("scale", affect_columns = selector_type("numeric"), center = TRUE, id = "numeric_scaling"))) %>>%  # 2) scaling continuous predictors
  po("featureunion")

graph$plot(horizontal = "TRUE")

# Tuning learners
require(mlr3mbo)

# defining optimisation algorithm: bayesian optimisation
tuner = tnr("mbo")

# defining resampling strategy for tuning
inner_resampling = rsmp("cv", folds = 3, id = "inner_cv")

# defining resampling strategy for validation
outer_resampling = rsmp("repeated_cv", folds = 5, repeats = 20, id = "outer_cv") # nested res (change repeats to 20)
# outer_resampling = rsmp("cv", folds = 5, id = "outer_cv")

# preprocessing each task
prepro_task.pn <- graph$train(task.pn)

tuned_xgb = auto_tuner(tuner = tuner,
                       learner = learner.xgb,
                       resampling = inner_resampling,
                       measure = msr("classif.tss"),
                       term_evals = 50, #50
                       store_models = F)

tuned_rf = auto_tuner(tuner = tuner,
                      learner = learner.rf,
                      resampling = inner_resampling,
                      measure = msr("classif.tss"),
                      term_evals = 30, #30
                      store_models = F)

tuned_gam = auto_tuner(tuner = tuner,
                       learner = learner.gam,
                       resampling = inner_resampling,
                       measure = msr("classif.tss"),
                       term_evals = 3, #5
                       store_models = F)

tuned_nn = auto_tuner(tuner = tuner,
                      learner = learner.nn,
                      resampling = inner_resampling,
                      measure = msr("classif.tss"),
                      term_evals = 30, #30
                      store_models = F)

# does work
set.seed(222)
tuned_xgb$train(prepro_task.pn[[1]])
tuned_rf$train(prepro_task.pn[[1]])
tuned_gam$train(prepro_task.pn[[1]])
tuned_nn$train(prepro_task.pn[[1]])

gr_stacking <- gunion(list(tuned_rf,
                           tuned_xgb,
                           tuned_gam,
                           tuned_nn)) %>>%
  po("classifavg", id = "w_average", innum = 4, param_vals = list(weights = c(tuned_rf$tuning_result %>% summarise(mean(classif.tss)) %>% .[[1]],
                                                                              tuned_xgb$tuning_result %>% summarise(mean(classif.tss)) %>% .[[1]],
                                                                              tuned_gam$tuning_result %>% summarise(mean(classif.tss)) %>% .[[1]],
                                                                              tuned_nn$tuning_result %>% summarise(mean(classif.tss)) %>% .[[1]])))

# gr_stacking$plot()
tuned_stack <- as_learner(gr_stacking)
tuned_stack$id <- "stacked.learner"
tuned_stack$predict_sets=c("test")

# Parallelization over resamples
library(future)
future::plan("multisession", workers = 20)
set.seed(1234567)
# Benchmarking experiment
bench_stack_results <- benchmark(benchmark_grid(tasks = prepro_task.pn,
                                                learners = tuned_stack,
                                                resamplings = outer_resampling),
                                 store_models = F)

saveRDS(bench_stack_results, file = ".../data/hpv_stack_model.rds")

#### Predicting out of sample (species with unknown status) ####
# Loading benchmark result (contains performance metrics and the tuned ensemble)
ensemble_model <- as.data.table(bench_stack_results)

# Loading dataset to predict
dataset_for_prediction <- read.csv(".../data/all_predictors.csv")
phylo <- read_csv(".../data/hpvmean_phylo_dist_to_known_hosts.csv")
dataset_for_prediction <- dataset_for_prediction %>% merge.data.frame(., phylo, by.x = 1, by.y = 1)
dataset_for_prediction$area.1 <- log10(dataset_for_prediction$area.1)

# selecting orders with at least 1 positive species
predictable_orders <- unique(dataset$order)

dataset_for_prediction2 <- dataset_for_prediction %>%
  filter(order %in% predictable_orders) %>% 
  dplyr::select(c(1, 2, 3, 4, 7, 11, 16, 18, 19, 
                  20, 22, 45, 61:79)) %>% 
  filter(complete.cases(.))

pred_tibble <- data.frame(matrix(ncol = 100, nrow = length(dataset_for_prediction2$phylacine_binomial)))

library(future)
future::plan("multisession", workers = 20)

set.seed(123)

for (i in c(12:100)){
  
  graph_learner <- as_learner(graph %>>% ensemble_model$learner[[i]])
  pred_tibble[i] <- as.data.table(graph_learner$train(task.pn)$predict_newdata(dataset_for_prediction2)) %>% dplyr::select(prob.1) %>% .[[1]]
  
}

predictions <- cbind(dataset_for_prediction2[c(1:4)], pred_tibble)

write.csv(predictions, ".../data/hpv_predictions_outofsample.csv", row.names = F)

#### Predictions for in sample species (high- and low- evidence hosts, pseudo-negatives) ####
dataset <- read.csv(".../data/hpv_dataset_main.csv")
bench_stack_results <- readRDS(".../data/hpv_stack_model.rds")
dataset_predictions_withpseudo <- as.data.table(bench_stack_results) %>% .[1:100]

dataset_predictions_withpseudo$prediction %>%
  lapply(as.data.table) %>%
  data.table::rbindlist(idcol = "fold") -> dataset_predictions_withpseudo

dataset_predictions_withpseudo$iteration <- rep(1:20, each = length(unique(dataset_predictions_withpseudo$row_ids)))

tab_rep <- NULL

for (i in c(1:20)){
  tab <- dataset_predictions_withpseudo %>%
    group_by(iteration) %>%
    arrange(row_ids, .by_group = T) %>%
    filter(iteration == i) %>%
    ungroup() %>%
    cbind(., dataset[c(1, 2, 3, 4, 32, 33)]) %>%
    filter(iucn2020_binomial %in% dataset$iucn2020_binomial)
  
  tab_rep <- rbind(tab_rep, tab)
}

# Putting them together
'%!in%' <- function(x, y)!('%in%'(x, y))
predictions_outofsample <- read.csv(".../data/hpv_predictions_outofsample.csv")
predoutsample <- predictions_outofsample %>% filter(iucn2020_binomial %!in% dataset$iucn2020_binomial) %>% rowwise() %>%
  mutate(mean_pred = mean(c_across(5:104), na.rm = T),
         sd_pred = sd(c_across(5:104), na.rm = T),
         hard_class = ifelse(mean_pred >= 0.5, 1, 0)) %>% select(1, 2, 3, 4, 105, 106, 107) %>% mutate(real_status = "unknown")
predinsample <- tab_rep %>% group_by(phylacine_binomial, iucn2020_binomial, family, order, status, real_status) %>% summarise(mean_pred = mean(prob.1), sd_pred = sd(prob.1), hard_class = ifelse(mean_pred > 0.5, 1, 0))

strict_pred_fam <- unique(predinsample %>% filter(real_status == "high-evidence host") %>% pull(family))
pred_fam <- unique(predinsample %>% filter(status == 1) %>% pull(family))

predictions_both <- rbind(predoutsample %>% ungroup(), predinsample %>% ungroup() %>% select(-5)) %>%
  mutate(predictable_fam = ifelse(family %in% pred_fam, 1, 0),
         strict_predictable_family = ifelse(family %in% strict_pred_fam, 1, 0))

write.csv(predictions_both, ".../data/hpv_predictions_main.csv", row.names = F)
