#!/usr/bin/env R

# for generic functions used often
# ported from gen_sapt_efp.r
# and scs_rad.r etc

# IN PARTICULAR FOR SCS PROJECTS
# functions copied from sapt_efp
# will not be removed from original
# scripts, instead don't create 
# separate copies in other SCS files

# ONLY FUNCTIONS THAT WILL NOT BE ALTERED 
# FUNCTIONS LIKE coef.scs.fn WHICH WILL BE
# TWEAKED FOR PARTICULAR PROJECTS AND PURPOSES
# E.G. BOOTSTRAPPING, SHOULD BE LEFT IN THEIR
# ORIGINAL SCRIPTS
# THIS WAY, CHANGES TO THOSE FUNCTIONS WILL
# NOT IMPACT OTHER PROJCETS

# STILL UP IN THE AIR ABOUT THE BOOTSTRAP
# FUNCTION, AS IT DOESN'T CHANGE MUCH 

classic.stats <- function(x) {
    c(mean = mean(x, na.rm = TRUE)
      , med = median(x, na.rm = TRUE)
      , sd = sd(x, na.rm = TRUE)
      , min = min(x, na.rm = TRUE)
      , max = max(x, na.rm = TRUE))
}

# mean average error
mae <- function(x) {
    mean(abs(x), na.rm = TRUE)
}

# gives mae, med, sd, min, max (min.name, max.name)
basic.stats <- function(x, iden = NULL) {
    
    if (missing(iden)) {
        a = c(MAE = mae(x)
                 , Med = median(x, na.rm = TRUE)
                 , SD = sd(x, na.rm = TRUE),
                 Min = min(x, na.rm = TRUE)
                 , Max = max(x, na.rm = TRUE))
    } else {
        iden = as.list(iden)
        a = c(MAE = mae(x)
                 , Med = median(x, na.rm = TRUE)
                 , SD = sd(x, na.rm = TRUE)
                 , Min = min(x, na.rm = TRUE)
                 , Max = max(x, na.rm = TRUE)
                 , Min.Name = iden[which.min(x)]
                 , Max.Name = iden[which.max(x)])
    }
    #return(as.list(a))
    return(a)
}

mue.stats <- function(x, iden = NULL) {
    if (missing(iden)) {
        return(as.list(basic.stats(x)))
    } else {
        return(as.list(basic.stats(x = x, iden = iden)))
    }
}



# gives mae, sd, max (max.name)
bare.stats <- function(x, iden = NULL) {
    # somehow need double square brackets to pick up number only
    # otherwise rowname will show up in max column name
    # not sure why this doesn't affect basic.stats
    
    if (missing(iden)) {
        mymax = which.max(abs(x))
        a = c(MAE = mae(x)
                 , SD = sd(x, na.rm = TRUE)
                 , Max = x[[mymax]])
    } else {
        iden = as.list(iden)
        a = c(MAE = mae(x), SD = sd(x, na.rm = TRUE)
                 , Max = x[[which.max(abs(x))]]
                 , Max.Name = iden[which.max(abs(x))])
    }
    return(as.list(a))
}

# no median
nomed.stats <- function(x) {
    a = c(MAE = mae(x)
          , SD = sd(x, na.rm = TRUE),
          Min = min(x, na.rm = TRUE)
          , Max = max(x, na.rm = TRUE))
    return(a)
}


# from gen_sapt_efp.r
coef.simple.fn <- function(dt, x, y, iden.name = NULL) {
    # expand dependent x's & construct formula
    dep = paste(y, collapse = " + ")
    my.formula = paste(x, " ~ ", dep, sep = "")
    
    m = lm(formula = my.formula, data = dt, na.action = na.exclude)
    model = summary(m)
    
    coefficients = coef(m)
    setNames(coefficients, names(coef(m)))
    # return vector of interested values
    #return(c(coef(model), adj_r2, basic.stats(resid(model))))
    # model$adj.r.squared 
    if (missing(iden.name)) {
        a = c(coefficients, "R.sq" = model$r.squared, basic.stats(resid(model)), n = nrow(dt))
    } else {
        a = c(coefficients, "R.sq" = model$r.squared, basic.stats(resid(model), iden.name), n = nrow(dt))
    }
    return(as.list(a))
}


coef.stats.fn <- function(dt, x, y, iden.name = NULL) {
    # expand dependent x's & construct formula
    dep = paste(y, collapse = " + ")
    my.formula = paste(x, " ~ ", dep, sep = "")
    
    model = summary(lm(formula = my.formula, data = dt, na.action = na.exclude))
    
    coefficients = coef(model)
    # create vector of names for named list
    coef.names = dimnames(coefficients)[1]
    other.stats = c("coef", "stderr", "tval", "pval")
    l.n = c()
    #coef.names = c("OS.ncp", "SS.ncp")
    for (i in other.stats) {
        l.n = c(l.n, mapply(paste, coef.names, i, sep = "_"))
    }
    # flatten table
    coefficients = c(coefficients)
    # name the coefficients
    names(coefficients) = l.n
    # return vector of interested values
    #return(c(coef(model), adj_r2, basic.stats(resid(model))))
    # model$adj.r.squared 
    if (missing(iden.name)) {
        a = c(coefficients, "R.sq" = model$r.squared, basic.stats(resid(model)), n = nrow(dt))
    } else {
        a = c(coefficients, "R.sq" = model$r.squared, basic.stats(resid(model), iden.name), n = nrow(dt))
    }
    return(as.list(a))
}

fast.stats <- c("OS.nonCP", "SS.nonCP", "MAE", "SD", "Min", "Max")

fast.coef.fn <- function(dt, indices, x, y) {
    # indep and dep swapped, old mistake
    dt = dt[indices, ]
    tindep = data.matrix(dt[, x, with = F])
    tdep = data.matrix(dt[, y, with = F])
    
    mod = lm.fit(x = tdep, y = tindep)
    cc = mod$coefficients
    
    numbers = c(cc, nomed.stats(resid(mod)))
    
    #return(as.list(numbers))
    return(numbers)
}

predict.fn <- function(dt, x, y) {
    # copied from above
    dep = paste(y, collapse = " + ")
    my.formula = paste(x, " ~ ", dep, sep = "")
    
    model = lm(formula = my.formula, data = dt, na.action = na.exclude)
    
    return(data.table(predict(model, newdata = dt)))
}

scaled.fn <- function(dt, x, y, i) {
    # assume no intercept, formula of form a + b + 0
    # i is the name (specified in y) 
    # of the column to scale in dt
    # it is a string, e.g. "OS.ncp"
    dep = paste(y, collapse = " + ")
    my.formula = paste(x, " ~ ", dep, sep = "")
    
    m = lm(formula = my.formula, data = dt, na.action = na.exclude)
    
    # column names in dt correspond to y
    return(dt[, get(i) * coef(m)[i]])
}

# function to factorize integers
factor_integer <- function(x) {
    x = as.integer(x)
    div = seq_len(abs(x))
    factors = div[x %% div == 0L]
    factors = list(neg = -factors, pos = factors)
    return(factors)
}


### ============= BOOTSTRAP FUNCTION =============== ###
# copied from bootstrap.r from the sapt_efp.Rproj 
# function to generate  (99%) confidence intervals
# NOTE: boot_comp takes in list of interested statistics and their positions returned by the statistic fn
# paired tuple of names of statistic to bootstrap & it's position
ci.extract.fn <- function(dt, indep, dep, boot_comp, stat.fn) {
    boot.out <- boot(data = dt, R = n_iter, statistic = stat.fn,
                     x = indep, y = dep,
                     parallel = "multicore", ncpus = 12)
    n_stats = length(boot_comp)                     # number of statistics, length of first row
    boot.stat.dt <- data.table(Component = as.factor(boot_comp), 
                               stat = vector(mode = "numeric", length = n_stats),
                               lower = vector(mode = "numeric", length = n_stats),
                               upper = vector(mode = "numeric", length = n_stats))
    # get value, and upper and lower bounds of 95% confidence interval for interested stats
    # output of boot.ci looks like
    # t* mean sd min max
    for (i in 1:n_stats) {
        ### REMEMBER TO CHANGE CI.OUT$BCA/BASIC WHEN CHANGING TYPE!!! ###
        ci.out <- boot.ci(boot.out, index = as.integer(i), type = "basic", conf = 0.99)
        int_val <- c(ci.out$t0, ci.out$basic[4:5])  #can't use tail for some reason
        boot.stat.dt[Component == boot_comp[i], c("stat", "lower", "upper") := as.list(int_val)]
    }
    return(boot.stat.dt)
}

scs.bs.fn <- function(dt, indices, x, y) {
    # SAPT (x) on vertical, EFP (y) taking on coefficient
    # expand dependent x's & construct formula
    dep = paste(y, collapse = " + ")
    my.formula = paste(x, " ~ ", dep, sep = "")
    
    d = dt[indices, ]
    
    # using summary since resid() can extract residuals from both 
    # lm() and summary(lm())
    # and other two groups need summary()
    model <- summary(lm(formula = my.formula, data = d))
    
    return(c(coef(model), model$r.squared, basic.stats(resid(model))))
}

wide.ci.printfn <- function(dt) {
    d = gather(data = dt, key = Stats, value = values, ... = stat:upper) 
    d$Component = factor(d$Component, levels = l.stats)
    d$Stats = factor(d$Stats, levels = c("stat", "lower", "upper"))
    print(d)
    d = unite(data = d, col = comp.stat, ... = c(Component, Stats), sep = ".") %>%
        spread(key = comp.stat, value = values) %>%
        data.table()
    return(d)
}

collate <- function(a,b,c) {
    # stat, lower and upper for CI data.table's
    sprintf("% .3g, (% .3g, % .3g)", a, b, c) 
}

perc_err.fn <- function(dt, indep, dep) {
    # calculates the percentage error
    d = copy(dt)
    d[, predEn := predict.fn(dt = .SD, x = indep, y = dep)]
    a = d[, bare.stats(-100 * abs(corrEn - predEn) / corrEn)]
    return(a)
}

# for plotting different conditions
# zeroN and zeroA now that corrEn < -4 no longer used
# function to match cond to appropriate title according to list given
cond_to_title <- function(cond, inList) {
    return(inList[[cond]])
}

# for plotting points side by side
#dodge <- position_dodge(width = 0.02)

# first written in "perc_err.r", appropriated here 
# since it's first mention is used here (Mon 02 May 2016)
interested <- c("os.coef", "ss.coef", "MAE", "SD", "Min", "Max")


# for plotting ci plots of the pointrange type
# moved here from srs.r (Mon 02 May 2016)
srs.ci.pltfn <- function(dt, titleStr, ggsave.prefix) {
    plt <- ggplot(data = dt[Component %in% interested]) +
        geom_pointrange(aes(x = cutoff, y = stat, ymin = lower, ymax = upper, colour = r, shape = method)
                        , position = dodge, size = 0.5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_grid(Component ~ basis, scales = "free_y") +
        ggtitle(titleStr)
    print(plt)
    
    ggsave(plot = plt, path = "~/GoogleDrive/scs-it/images/cutoffs/"
           , filename = paste(ggsave.prefix, ".pdf", sep = "")
           , width = 12 ,height = 7)
}