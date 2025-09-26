#BMI_siginicance in saliva

t_env <- trans_env$new(dataset = dataset_saliva, 
                       add_data = dataset_saliva$sample_table["BMI"])
t_env$cal_diff(group = "Group", method = "t.test")
t_env$res_diff

write.csv(t_env$res_diff, file ="BMI_significance_saliva",row.names = TRUE)
#====
#BMI_siginicance in gut

t_env <- trans_env$new(dataset = dataset_gut, 
                       add_data = dataset_gut$sample_table["BMI"])
t_env$data_env[] <- lapply(t_env$data_env, function(x) as.numeric(as.character(x)))
t_env$data_env$BMI
t_env$cal_diff(group = "Group", method = "t.test")
t_env$res_diff

write.csv(t_env$res_diff, file ="BMI_significance_gut",row.names = TRUE)
