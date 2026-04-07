# Orthonairovirus
require(tidyverse)
require(mlr3)

dataset <- read.csv(".../data/ohv_dataset_main.csv")

# cor(dataset[c(5:9, 11:24)]) %>% 
#   as.data.frame(.) %>%
#   rownames_to_column(var = "variable1") %>% 
#   pivot_longer(cols = 2:20, names_to = "variable2", values_to = "correlation") %>% 
#   arrange(desc(correlation)) %>% 
#   filter(correlation > 0.7 | correlation < -0.7) %>% View(.)
# Remove variables with a pairwise |correlation| >= 0.7

dataset <- dataset[-c(10, 7, 22, 14, 25, 23, 15, 18, 19, 20, 17, 27, 28, 9, 6, 11, 30)]

# Variable preparation
dataset$status <- as.factor(as.character(dataset$status))
dataset$status <- fct_relevel(dataset$status, "1")
dataset$trophic_level <- as.factor(as.character(dataset$trophic_level))
dataset$area.1 <- log10(dataset$area.1)
dataset$order <- as.factor(dataset$order)
dataset$order2 <- as.factor(dataset$order)
dataset$real_status <- as.factor(dataset$real_status)
sum(dataset$w)

sum(dataset$w)

require(mlr3learners)

# Defining tasks
#one with pseudo-negatives
task.pn = as_task_classif(id = "pn", reformulate(names(dataset[-c(1, 2, 4, 15)]), names(dataset[15])), data = dataset, positive = "1")
task.pn$set_col_roles("w", roles = "weights_learner")
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
                  eta = to_tune(p_dbl(lower = 0.1, upper = 0.3)),
                  # max_depth = to_tune(p_int(1, 15)),
                  subsample = to_tune(p_dbl(0.5, 0.9)),
                  min_child_weight = to_tune(p_dbl(0.01, 0.1)),
                  predict_type = "prob")

learner.nn = lrn("classif.nnet",
                 size = to_tune(p_int(lower = 1, upper = 10)),
                 maxit = to_tune(p_int(10, 150)),
                 decay = to_tune(p_dbl(1e-7, 0.1, logscale = TRUE)),
                 predict_type = "prob")

require(mlr3extralearners)
learner.gam = lrn("classif.gam",
                  formula = status ~ + order.Cetartiodactyla + order.Chiroptera + order.Rodentia + order.Eulipotyphla+ order.Primates+ order.Perissodactyla+ order.Proboscidea+ order.Lagomorpha+ order.Carnivora+
                    trophic_level.3 + trophic_level.2 + s(adult_mass_g) +
                    s(litter_size_n) +
                    s(mean.wc2.1_5m_bio_9) + s(mean.wc2.1_5m_bio_19)+
                    s(stdev.wc2.1_5m_bio_9)+
                    s(stdev.wc2.1_5m_bio_16)+
                    s(stdev.wc2.1_5m_bio_19)+s(area.1),
                  select = TRUE,
                  method = "REML",
                  gamma = to_tune(p_dbl(1, 5)),
                  predict_type = "prob")

# Define new performance measure (TSS)
#template code from mlr3book.mlr-org.com Chapter 10
MeasureClassifTSS = R6::R6Class("MeasureClassifTSS",
                                inherit = mlr3::MeasureClassif,
                                public = list(
                                  initialize = function() { # initialize class
                                    super$initialize(
                                      id = "classif.tss", 
                                      packages = character(), 
                                      properties = character(), 
                                      predict_type = "response", 
                                      range = c(-1, 1), 
                                      minimize = FALSE 
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

# stacking_graph <- ppl("stacking",
#                       base_learners = c(tuned_gam,
#                                         tuned_nn,
#                                         tuned_rf,
#                                         tuned_xgb),
#                       super_learner = lrn("classif.log_reg"),
#                       id = "stack")

# Stacking models
# cool but does not work
# gr_stacking <- gunion(list(tuned_rf %>>% mlr_pipeops$get("copy", outnum = 2, id = "copy_rf"),
#                    tuned_xgb %>>% mlr_pipeops$get("copy", outnum = 2, id = "copy_xgb"))) %>>%
#   gunion(list(po("nop", id = "nope_xgb"),
#               po("classifavg", id = "w_average", innum = 2, param_vals = list(weights = 1)),
#               po("nop", id = "nope_rf"))) %>>%
#   po("unbranch")
#
# tuned_stack <- as_learner(gr_stacking)

# does work
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


# gr_stacking = gunion(list(tuned_rf %>>% gunion(list(mlr_pipeops$get("copy", outnum = 1, id = "copy_rf"),
#                                                     po("classifavg", id = "stack", innum = 2, param_vals = list(weights = 1)))),
#                           tuned_xgb %>>% mlr_pipeops$get("copy", outnum = 1, id = "copy_xgb") %>>% po("select", id = "select_xgb"))) %>>%
#   gunion(list(po("classifavg", id = "stack", innum = 2, param_vals = list(weights = 1))))

# res_results <- resample(prepro_task.pn[[1]], tuned_stack, outer_resampling, store_models = TRUE)
#
# res_results$score(msrs(c("classif.auc", "classif.prauc", "classif.tnr", "classif.tpr", "classif.npv", "classif.ppv", "classif.mcc", "classif.tss"))) %>%
#   summarise(mean_TSS = mean(classif.tss),
#             sd_TSS = sd(classif.tss),
#             mean_AUC = mean(classif.auc, na.rm = T),
#             sd_AUC = sd(classif.auc),
#             mean_PRAUC = mean(classif.prauc,  na.rm = T),
#             sd_PRAUC = sd(classif.prauc),
#             mean_TPR = mean(classif.tpr),
#             sd_TPR = sd(classif.tpr),
#             mean_TNR = mean(classif.tnr),
#             sd_TNR = sd(classif.tnr),
#             mean_NPV = mean(classif.npv),
#             sd_NPV = sd(classif.npv),
#             mean_PPV = mean(classif.ppv),
#             sd_PPV = sd(classif.ppv),
#             mean_MCC = mean(classif.mcc),
#             sd_MCC = sd(classif.mcc))
#
# extract_inner_tuning_results(res_stack_results)[,list(learner_id, classif.tss)]
#
# # extract stuff (i.e., metrics and predictions)
# res_stack_results$model$learner$graph_model$pipeops$ranger.ranger$predict(list(task))
#
# autoplot(res_results, type = "boxplot", measure = msr("classif.auc"))


# mlr3misc::map(as.data.table(res_stack_results)$learner, "model") -> data
#
# pred <- res_stack_results$prediction(predict_sets = "test")
# data <- as.data.table(res_stack_results, reassemble_learners = TRUE, convert_predictions = TRUE, predict_sets = "test")
#
# map(as.data.table(res_stack_results)$learner, "model") -> g_mod
# g_mod[[1]]$classif.ranger.tuned$model$learner$model$
#
# $graph_model$pipeops$classif.ranger.tuned$output$
# graph$pipeops$variance$.result$output$data()
# map(data$learner, "model") -> data_models
# map(data_models, "classif.ranger.tuned") -> ranger_mod
#
# map(ranger_mod, "model") -> ranger_mod2
#
# map(ranger_mod2, "learner") -> ranger_lrn
#
# res_data <- as.data.table(res_stack_results)
#
# res_stack_results$learners[[1]]$model$classif.ranger.tuned$model$learner$
#
# res_stack_results$learners
#
# data[[1]]$classif.ranger.tuned
# data[[1]]$classif.xgboost.tuned$model$tuning_instance$archive$benchmark_result$score()
# data[[1]]$classif.xgboost.tuned$model$learner$model
#
# as.data.table(as_benchmark_result(res_stack_results)) -> bench_table
# map(bench_table$learner, "model") -> model_tab
#
# map(model_tab, "classif.nnet.tuned") -> xgboost_tab

# res_stack_results$learner$graph_model$pipeops$classif.ranger.tuned$learner$learner$model

# data[[1]]$classif.ranger.tuned$model$tuning_instance$archive
# data[[1]]$classif.ranger.tuned$model$learner$model$predictions

# # Check predictions
# learner <- map(data$learner, "learner")
# tuning <- map(data$learner, "tuning_result")
# archive <- map(data$learner, "archive")

# Parallelization over resamples
library(future)
future::plan("multisession", workers = 20)
set.seed(2025)
# Benchmarking experiment
bench_stack_results <- benchmark(benchmark_grid(tasks = prepro_task.pn,
                                                learners = tuned_stack,
                                                resamplings = outer_resampling),
                                 store_models = F)
# Resampling experiment (only stack bu can't extract performance of individual classifiers)
# res_stack_results <- resample(task = prepro_task.pn[[1]], learner = tuned_stack, resampling = outer_resampling, store_models = TRUE)

saveRDS(bench_stack_results, file = ".../data/ohv_stack_model.rds")

#### Predicting out of sample (species with unknown status) ####
# Loading benchmark result (contains performance metrics and the tuned ensemble)
ensemble_model <- as.data.table(bench_stack_results)

# Loading dataset to predict
dataset_for_prediction <- read.csv(".../data/all_predictors.csv")
phylo <- read_csv(".../data/ohvmean_phylo_dist_to_known_hosts.csv")
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

set.seed(42)

for (i in c(1:100)){
  
  graph_learner <- as_learner(graph %>>% ensemble_model$learner[[i]])
  pred_tibble[i] <- as.data.table(graph_learner$train(task.pn)$predict_newdata(dataset_for_prediction2)) %>% dplyr::select(prob.1) %>% .[[1]]
  
}

predictions <- cbind(dataset_for_prediction2[c(1:4)], pred_tibble)

write.csv(predictions, ".../data/ohv_predictions_outofsample.csv", row.names = F)
# bench_stack_results <- readRDS(file = "ohv_stack_mod.rds")

#### Predictions for in sample species (high- and low- evidence hosts, pseudo-negatives) ####
dataset <- read.csv(".../data/ohv_dataset_main.csv")
bench_stack_results <- readRDS(".../data/ohv_stack_model.rds")
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
predictions_outofsample <- read.csv(".../data/ohv_predictions_outofsample.csv")
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

write.csv(predictions_both, ".../data/ohv_predictions_main.csv", row.names = F)
