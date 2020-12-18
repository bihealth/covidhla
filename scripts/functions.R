comb_fisher <- function(p) {

  k <- length(p)
  chi2 <- -2 * sum(log(p))
  return(pchisq(chi2, df=2*k, lower.tail=FALSE))
}


hla2wide <- function(hla, id="ID", a1="allele1", a2="allele2", minN=0) { 

  tmp <- hla %>% mutate(a1=paste0(locus, "*", .[[a1]]), a2=paste0(locus, "*", .[[a2]])) %>% select(ID=all_of(id), a1, a2) %>% 
    pivot_longer(cols=a1:a2, names_to="a1a2", values_to="allele") %>% 
    filter(!duplicated(paste0(ID, allele))) %>% select(ID, allele) %>% 
    mutate(value=1) %>% 
    pivot_wider(id_cols="ID", names_from="allele", values_from=value) 
  tmp[ is.na(tmp) ] <- 0

  hla_wide <- tmp %>% select(-ID) %>% select_if(~ sum(.) > minN)

  hla_wide <- cbind(ID=tmp$ID, hla_wide)
  return(hla_wide)
}


adjust_by_freq <- function(pval, freq, N, parameter, min.freq=0.05, min.patient=5) {

  sel <- freq >= min.freq & (1-freq) >= min.freq & N >= min.patient

  ret <- rep(NA, length(pval))
  pars <- unique(parameter)
  for(p in pars) {
    sel2 <- p == parameter
    ret[sel & sel2] <- p.adjust(pval[sel & sel2])
  }
  return(ret)
}


freq_allele <- function(locus, a1, a2) {

  hla <- data.frame(locus=locus, a1=a1, a2=a2)
  hla %>% pivot_longer(cols=2:3, names_to="pos", values_to="allele") %>% 
    drop_na %>%
    select(locus, allele) %>% group_by(locus) %>%
    mutate(locn=n()) %>% group_by(locus, allele) %>% mutate(n=n(), freq=n()/locn) %>% 
    filter(!duplicated(paste0(locus, allele))) %>% ungroup() %>%
    mutate(allele_n=paste(locus, allele))
}

run_assoc_test_new <- function(hla, traits, loci=NULL, patients=NULL, id="ID", allele1="allele1_trunc", allele2="allele2_trunc", model="d", log=TRUE, ctest="t.test") {
  cat("run_assoc_test:")
  if(is.null(loci)) {
    loci <- unique(hla$locus)
  }
  if(!is.null(patients)) {
    hla <- merge(hla, patients, by=id)
  }

  if(!all(traits %in% colnames(hla))) {
    stop("All traits must be in hla/patients data frames")
  }
  res <- map(set_names(traits), ~ {
    trait <- .
    message(trait)
    cat(".")
    res <- map(set_names(loci), hla_assoc_test_new, hla, ., log=log, id=id, allele1=allele1, allele2=allele2, ctest=ctest)
    res <- res[!map_lgl(res, is.null)]
  })
  cat("\n")
  res <- res[!map_lgl(res, is.null)]
  return(res)
}

run_assoc_test <- function(hla, traits, loci=NULL, patients=NULL, id="ID", allele1="allele1_trunc", allele2="allele2_trunc", model="d") {
  cat("run_assoc_test:")
  if(is.null(loci)) {
    loci <- unique(hla$locus)
  }
  if(!is.null(patients)) {
    hla <- merge(hla, patients, by=id)
  }

  if(!all(traits %in% colnames(hla))) {
    stop("All traits must be in hla/patients data frames")
  }
  res <- map(set_names(traits), ~ {
    trait <- .
    cat(".")
    res <- map(set_names(loci), hla_assoc_test, hla, ., log=TRUE, id=id, allele1=allele1, allele2=allele2)
    res <- res[!map_lgl(res, is.null)]
  })
  cat("\n")
  res <- res[!map_lgl(res, is.null)]
  return(res)
}

get_assoc_summary <- function(res, min_freq, min_patient) {
  colsel <- c("allele", "parameter", "freq", "test", "[-/-]", "[-/h,h/h]", "T.[-/-]", "T.[-/h,h/h]", "E", "h.pval", "p.val", "p.adj")

  res <- res %>% 
    imap( ~ { 
      param <- .y 
      .x %>% imap(~ { .x[["locus"]] <- .y ; .x }) %>% 
        { Reduce(rbind, .) } %>% 
        mutate(parameter=param) 
    }) %>% { Reduce(vmerge, .) } %>%
    mutate("T.[-/-]"=ifelse(is.na(.data[["%.[-/-]"]]), .data[["avg.[-/-]"]], .data[["%.[-/-]"]])) %>%
    mutate("T.[-/h,h/h]"=ifelse(is.na(.data[["%.[-/h,h/h]"]]), .data[["avg.[-/h,h/h]"]], .data[["%.[-/h,h/h]"]])) %>%
    mutate("test"=ifelse(is.na(fisher.p), "ttest", "chisq")) %>%
    mutate(p.val=ifelse(is.na(fisher.p), ttest.p, fisher.p)) %>%
    mutate(p.adj=adjust_by_freq(p.val, freq, .data[["[-/h,h/h]"]], parameter, min.freq=min_freq, min.patient=min_patient)) %>%
    select(all_of(colsel))

  return(res)
}

get_assoc_summary_new <- function(res, min_freq, min_patient) {
  colsel <- c("locus", "allele", "parameter", "freq", "type", "[-/-]", "[-/h,h/h]", "T.[-/-]", "T.[-/h,h/h]", "E", "p.val", "p.adj")

  res <- res %>% 
    imap( ~ { 
      param <- .y 
      .x %>% imap(~ { .x[["locus"]] <- .y ; .x }) %>% 
        { Reduce(rbind, .) } %>% 
        mutate(parameter=param) 
    }) %>% { Reduce(vmerge, .) }

  if(!"%.[-/-]" %in% colnames(res)) {
    res[["%.[-/-]"]] <- NA
    res[["%.[-/h,h/h]"]] <- NA
  }
  if(!"avg.[-/-]" %in% colnames(res)) {
    res[["avg.[-/-]"]] <- NA
    res[["avg.[-/h,h/h]"]] <- NA
  }
     
  res <- res %>%
    mutate("T.[-/-]"=ifelse(is.na(.data[["%.[-/-]"]]), .data[["avg.[-/-]"]], .data[["%.[-/-]"]])) %>%
    mutate("T.[-/h,h/h]"=ifelse(is.na(.data[["%.[-/h,h/h]"]]), .data[["avg.[-/h,h/h]"]], .data[["%.[-/h,h/h]"]])) %>%
    mutate(p.adj=adjust_by_freq(p.val, freq, .data[["[-/h,h/h]"]], parameter, min.freq=min_freq, min.patient=min_patient)) %>%
    select(all_of(colsel))

  return(res)
}




#' Association test for a single trait
#' @param locus.sel selected locus, e.g. HLA_A
#' @param hla data frame with alleles and traits
#' @param trait name of the column for the association test
#' @param allele1, allele2 names of the columns with the alleles
#' @param model which type of model (dominant, additive, recessive, genotype)
#' @return data frame with results
hla_assoc_test_new <- function(locus.sel, hla, trait, id="Patient_ID", log=FALSE, allele1="allele1_trunc", allele2="allele2_trunc", model="d", ctest="t.test") {
  hla <- hla %>% filter(locus == locus.sel) 
  #message(sprintf("HLA assoc test for variable %s and locus %s", trait, locus.sel))

  if(is.null(hla[[trait]])) stop(paste0("no such trait: ", trait))

  type <- "continuous"


  if(!is.numeric(hla[[trait]])) {
    #message("Converting trait ", trait, " to 0/1")
    type <- "categorical"
    hla[[trait]] <- factor(hla[[trait]])
    if(length(levels(hla[[trait]])) > 2) {
      stop(paste0("Too many levels for ", trait))
    }
    hla[[trait]] <- as.numeric(hla[[trait]]) - 1
  } else if(log) {
    x <- hla[[trait]]
    sel <- x == 0
    x[sel] <- min(x[!sel], na.rm=TRUE)
    hla[[trait]] <- log(x)
  }

  hla <- hla[ !is.na(hla[[trait]]) & !is.na(hla[[allele2]]) & !is.na(hla[[allele1]]), ]
  if(nrow(hla) < 5) {
    warning("Too few samples (<5) with this trait and complete alleles")
    return(NULL)
  }

  hla[["test"]]      <- hla[[trait]]
  hla[["sample.id"]] <- hla[[id]]
  hla[["locus"]] <- locus.sel

# allele_freq <- freq_allele(hla[["locus"]], hla[[allele1]], hla[[allele2]]) %>%
#   mutate(allele=paste(locus, allele)) %>% select(allele, freq, n, locn)


# HL <- hlaAllele(hla[["sample.id"]], H1=hla[[allele1]], H2=hla[[allele2]], locus=gsub("HLA_", "", locus.sel), assembly="hg38")
#
# res <- hlaAssocTest(HL, test ~ h, data=hla, model=model, verbose=FALSE) %>% rownames_to_column %>% colorDF %>%
#   mutate(locus=locus.sel, allele=paste(locus, rowname))
#
  res2 <- hla_assoc_test_allele(hla[[allele1]], hla[[allele2]], hla[[trait]], model, type, ctest=ctest)
  if(is.null(res2)) {
    return(NULL)
  }
  res2 <- cbind(locus=locus.sel, type=type, res2) %>% mutate(allele=paste(locus, allele))

# if(!(any(res[[2]] > 5 & res[[3]] > 5))) return(NULL)
# if("h.pval" %in% colnames(res)) res <- res %>% arrange(h.pval)
# if("fisher.p" %in% colnames(res)) res <- res %>% arrange(fisher.p)
# if("ttest.p" %in% colnames(res)) res <- res %>% arrange(ttest.p)
# if("anova.p" %in% colnames(res)) res <- res %>% arrange(anova.p)
#
# res <- merge(res, allele_freq, by="allele", all.x=TRUE)
#
  res2
}

hla_assoc_test_fisher <- function(ret, levs, trait, var, all) {
  levs <- sort(set_names(levs, paste0("%.", levs)))
  ret <- cbind(ret, 
    map_dfr(all, ~ { .all <- .x ; 
           ret <- map_dbl(levs, ~ sum(var[[.all]] == .x & trait == "1", na.rm=TRUE)/sum(var[[.all]] == .x, na.rm=TRUE)) })) 
  #ret$E <- ret[[names(levs)[2]]] / ret[[names(levs)[1]]]
  ret_2 <- map(all, ~ {
         ct <- table(var[[.x]], trait)
         if(!(nrow(ct) > 1 && ncol(ct) > 1)) {
           return(data.frame(E=NA, p.val=NA))
         }
         chsq <- chisq.test(ct)
         return(data.frame(
          E=sqrt(chsq$statistic/sum(chsq$observed)),
          p.val=fisher.test(ct)$p.value
          ))
  }) %>% reduce(rbind)
  ret <- cbind(ret, ret_2)
  return(ret)
}

hla_assoc_test_wilcox <- function(ret, levs, trait, var, all) {
  levs <- set_names(levs, paste0("avg.", levs))
  ret <- cbind(ret, 
    map_dfr(all, ~ { .all <- .x ; 
           ret <- map_dbl(levs, ~ mean(trait[ var[[.all]] == .x], na.rm=TRUE)) })) 
  ret_2 <- map(all, ~ {
    vv <- var[[.x]]
    ss <- summary(factor(vv))
    if(length(ss) != length(levs) || any(summary(factor(vv)) < 3)) {
      return(data.frame(E=NA, p.val=NA))
    }

    dd <- data.frame(trait=trait, var=var[[.x]])

    return(data.frame(
      E    =rstatix::wilcox_effsize(dd, trait ~ var)$effsize,
      p.val=rstatix::wilcox_test   (dd, trait ~ var)$p))
  }) %>% reduce(rbind)

  ret <- cbind(ret, ret_2)
  return(ret)
}


hla_assoc_test_ttest <- function(ret, levs, trait, var, all) {
  levs <- set_names(levs, paste0("avg.", levs))
  ret <- cbind(ret, 
    map_dfr(all, ~ { .all <- .x ; 
           ret <- map_dbl(levs, ~ mean(trait[ var[[.all]] == .x], na.rm=TRUE)) })) 
  ret_2 <- map(all, ~ {
    vv <- var[[.x]]
    ss <- summary(factor(vv))
    if(length(ss) != length(levs) || any(summary(factor(vv)) < 3)) {
      return(data.frame(E=NA, p.val=NA))
    }
    return(data.frame(E=effsize::cohen.d(trait ~ var[[.x]], na.rm=TRUE)$estimate,
           p.val= t.test(trait ~ var[[.x]])$p.value))
  }) %>% reduce(rbind)

  ret <- cbind(ret, ret_2)
  return(ret)
}

hla_assoc_test_allele <- function(allele1, allele2, trait, model, type, ctest="t.test") {

  sel <- !is.na(trait)
  allele1 <- allele1[sel]
  allele2 <- allele2[sel]
  trait <- trait[sel]
  all <- unique(c(allele1, allele2))

  zygo <- map(set_names(all), ~ ifelse(allele1 == .x & allele2 == .x, "[h/h]", ifelse(allele1 == .x | allele2 == .x, "[-/h]", "[-/-]")))

  var <- map(zygo, ~ switch(model,
    d=ifelse(.x %in% c("[-/h]", "[h/h]"), "[-/h,h/h]", "[-/-]"),
    r=ifelse(.x == "[h/h]", .x, "[-/-,-/h]"),
    g=.x))

  levs <- set_names(unique(unlist(var)))
  if(length(levs) < 2) {
    return(NULL)
  }
  ret <- cbind(allele=all, map_dfr(all, ~ { .all <- .x ; map_int(levs, ~ sum(var[[.all]] == .x, na.rm=TRUE)) })) %>%
    mutate(freq=map_dbl(zygo, ~ (sum(.x == "[-/h]") + sum(.x == "[h/h]") * 2) / (2 * length(.x)))) %>%
    mutate(carrier_freq=map_dbl(zygo, ~ sum(.x %in% c("[-/h]", "[h/h]")) / length(.x)))

  if(type == "categorical") {
    ret <- hla_assoc_test_fisher(ret, levs, trait, var, all)
  } else {
    if(ctest == "t.test") {
      ret <- hla_assoc_test_ttest(ret, levs, trait, var, all)
    } else {
      ret <- hla_assoc_test_wilcox(ret, levs, trait, var, all)
    }
  }

  return(ret)

}



#' Association test for a single trait
#' @param locus.sel selected locus, e.g. HLA_A
#' @param hla data frame with alleles and traits
#' @param trait name of the column for the association test
#' @param allele1, allele2 names of the columns with the alleles
#' @param model which type of model (dominant, additive, recessive, genotype)
#' @return data frame with results
hla_assoc_test_old <- function(locus.sel, hla, trait, id="Patient_ID", log=FALSE, allele1="allele1_trunc", allele2="allele2_trunc", model="d") {
  hla <- hla %>% filter(locus == locus.sel) 
  #message(sprintf("HLA assoc test for variable %s and locus %s", trait, locus.sel))

  if(is.null(hla[[trait]])) stop(paste0("no such trait: ", trait))

  if(!is.numeric(hla[[trait]])) {
    #message("Converting trait ", trait, " to 0/1")
    hla[[trait]] <- factor(hla[[trait]])
    if(length(levels(hla[[trait]])) > 2) {
      stop(paste0("Too many levels for ", trait))
    }
    hla[[trait]] <- as.numeric(hla[[trait]]) - 1
  } else if(log) {
    x <- hla[[trait]]
    sel <- x == 0
    x[sel] <- min(x[!sel], na.rm=TRUE)
    hla[[trait]] <- log(x)
  }

  hla <- hla[ !is.na(hla[[trait]]) & !is.na(hla[[allele2]]) & !is.na(hla[[allele1]]), ]
  if(nrow(hla) < 5) {
    warning("Too few samples (<5) with this trait and complete alleles")
    return(NULL)
  }

  hla[["test"]]      <- hla[[trait]]
  hla[["sample.id"]] <- hla[[id]]
  hla[["locus"]] <- locus.sel

  allele_freq <- freq_allele(hla[["locus"]], hla[[allele1]], hla[[allele2]]) %>%
    mutate(allele=paste(locus, allele)) %>% select(allele, freq, n, locn)


  HL <- hlaAllele(hla[["sample.id"]], H1=hla[[allele1]], H2=hla[[allele2]], locus=gsub("HLA_", "", locus.sel), assembly="hg38")

  res <- hlaAssocTest(HL, test ~ h, data=hla, model=model, verbose=FALSE) %>% rownames_to_column %>% colorDF %>%
    mutate(locus=locus.sel, allele=paste(locus, rowname)) %>%
    mutate(freq_me=allele_freq$freq[ match(allele, allele_freq$allele) ])
  if(!(any(res[[2]] > 5 & res[[3]] > 5))) return(NULL)
  if("h.pval" %in% colnames(res)) res <- res %>% arrange(h.pval)
  if("fisher.p" %in% colnames(res)) res <- res %>% arrange(fisher.p)
  if("ttest.p" %in% colnames(res)) res <- res %>% arrange(ttest.p)
  if("anova.p" %in% colnames(res)) res <- res %>% arrange(anova.p)

  res <- merge(res, allele_freq, by="allele", all.x=TRUE)
  res
}

hla_assoc_test <- hla_assoc_test_old

vmerge <- function(x, y) {
  xcol <- colnames(x)
  ycol <- colnames(y)

  newx <- setdiff(ycol, xcol)
  for(c in newx) {
    x[[c]] <- NA
  }

  newy <- setdiff(xcol, ycol)
  for(c in newy) {
    y[[c]] <- NA
  }

  allc <- union(xcol, ycol)
  x <- x[ , allc, drop=FALSE]
  y <- y[ , allc, drop=FALSE]

  return(rbind(x, y))
}

