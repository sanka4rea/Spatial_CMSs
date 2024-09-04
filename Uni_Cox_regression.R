# An example for univariate cox regression analyses of multi-features
dfs_time <- Clinical_info$rfs.delay
dfs_event <- Clinical_info$rfs.event
res <- Surv(dfs_time,dfs_event)
univ_formulas <- sapply(factor_name,function(x) as.formula(paste('res ~', x)))
univ_models <- lapply(univ_formulas, function(x){summary(coxph(x, data=Clinical_info) ) })
univ_results <- sapply(1:length(univ_models), function(i){
    x <- univ_models[[i]]
    p <- x$logtest[3]
    HR <- x$conf.int[1,c(1,3,4)]
    p <- as.numeric(p)
    aa <- c(HR,p)
    return(aa)
})
colnames(univ_results) <- factor_name
rownames(univ_results) <- c("HR","L95","H95","logtest_P")