suppressMessages(library('dplyr'))
suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("mvnfast"))

LOG_INFO <- "INFO"
LOG_WARNING <- "WARNING"
LOG_ERROR <- "ERROR"

LOG <- function(x, level=LOG_INFO) {
    str <- paste0("[%D %H:%M:%S - ", level, "]")
    cat(format(Sys.time(), str), x, "\n")
}

log_sum_exp <- function(a, b) {
    x <- c(a, b)
    # Computes log(sum(exp(x))
    # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
    if ( max(abs(x)) > max(x) ) {
        offset <- min(x)
    } else {
        offset <- max(x)
    }
    log(sum(exp(x - offset))) + offset
}

annotate_cred_set <- function(df, prb=0.90) {
    # the call to abs  function may seem weird, but there is numerical issues with 1 - cumsum for strange reason.
    # w/o it sometimes small negative numbers appear and throw off CS computation
    return (df %>% group_by(BLOCK) %>% mutate(NPIP=PIP / sum(PIP)) %>% group_by(BLOCK) %>% arrange(BLOCK, NPIP) %>%
            mutate(IN.CRED.SET=abs((1 - cumsum(NPIP))) <= prb) %>% select(-NPIP))
}

get_independent <- function(CHR, ID, P0, P1, regions) {
    # Sometimes labels are the same bc of bugs...
    P1 = P0 + 1

    # compute overlaps
    g_ranges <- GRanges(seqnames = CHR, ranges= IRanges::IRanges(start = P0, end = P1, names = ID))
    r_ranges <- GRanges(seqnames = regions$CHR, ranges= IRanges::IRanges(start = regions$START, end = regions$STOP))
    overlaps <- findOverlaps(g_ranges, r_ranges, select="arbitrary", maxgap=1e4)

    # some annotations don't overlap exactly
    # just get nearest for now
    if (sum(is.na(overlaps)) > 0) {
        subset <- g_ranges[is.na(overlaps)]
        miss <- nearest(subset, r_ranges, select="arbitrary")
        overlaps[is.na(overlaps)] <- miss
    }

    # prettify regions and output them as 'blocks'
    pranges <- paste0(regions$CHR, ":", regions$START, "..", regions$STOP)
    pranges[overlaps]
}

get_local_params <- function(cur.FILE, cur.ID, cur.MODEL, genos) {
    # load weights into matrix after QCing...
    Wlist <- lapply(1:length(cur.FILE), function(i) {
        load( cur.FILE[i] )
        wgt.matrix[is.na(wgt.matrix)] = 0

        # Match up the SNPs and weights
        m = match( snps[,2] , genos$bim[,2] )
        m.keep = !is.na(m)
        snps = snps[m.keep,]
        wgt.matrix = wgt.matrix[m.keep,]

        cur.genos = genos$bed[,m[m.keep]]
        cur.bim = genos$bim[m[m.keep],]

        # Flip WEIGHTS for mismatching alleles
        qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
        wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

        # Predict into reference
        mod = cur.MODEL[i]

        wgt <- tibble(SNP=rownames(wgt.matrix), WGT=as.double(wgt.matrix[, mod]))
        colnames(wgt)[2] <- cur.ID[i]
        wgt
    })

    W <- purrr::reduce(Wlist, full_join, by="SNP")
    snps <- W$SNP
    W$SNP <- NULL
    W[is.na(W)] <- 0
    W <- as.matrix(W)

    # check if we have single gene
    if (is.null(dim(W)) && length(W) > 0) {
        W <- t(t(W))
    }
    rownames(W) <- snps

    # scale weights and compute LD
    m <- match(rownames(W), genos$bim[, 2])
    X <- genos$bed[,m]
    S <- apply(X %*% W, 2, sd)
    if (length(S) == 1) {
        flags <- S != 0
        if (S != 0) {
            S <- t(t(1 / S))
        }
    } else {
        # drop genes with 0 genetic covariance
        flags <- S != 0
        S <- S[flags]
        W <- W[, flags]
        S <- diag(1 / S)
    }
    SW <- W %*% S
    LD <- cor(X)

    return (list(SW=SW, LD=LD, FLAGS=flags))
}

fine_map <- function(cur.Z, cur.ID, cur.W, cur.LD, prb, prior_chisq, intercept=F, compute.grad=F, posterior_check=0, tol=2.220446e-14, verbose=F) {

    # check length
    m <- length(cur.Z)
    if (m > 1) {
        S <- t(cur.W) %*% cur.LD %*% cur.W
    } else {
        S <- t(1)
    }
    rownames(S) <- cur.ID
    colnames(S) <- cur.ID

    if (!compute.grad) {
        cur.grad <- NA
        cur.hess <- NA
    } 

    if (m > 1) {
        # this isn't the optimal solution, but _is_ unbiased.
        # optimal would require computing GLS during REML solution
        # and potentially induce large overhead
        if (intercept) {
            tmp_v <- rowSums(t(cur.W) %*% cur.LD)
            inter.est <- (t(tmp_v) %*% cur.Z) / sum(tmp_v ** 2)
            inter.s2 <- sum((cur.Z - tmp_v * inter.est)^2) / (m - 1)
            inter.se <- sqrt(inter.s2 / sum(tmp_v ** 2))
            inter.z <- inter.est / inter.se
            inter.p <- 2 * pnorm(abs(inter.z), lower.tail=F)
            mu <- tmp_v * inter.est
            Zp <- (cur.Z - mu)
        } else {
            mu <- rep(0, m)
            inter.z <- 0
            inter.p <- 1
            Zp <- cur.Z
        }

        pip <- rep(0, m) 
        null.pip <- m * log(1 - prb)
        log.marginal <- null.pip # initialize with NULL prior-weighted BF

        if (compute.grad) {
            cur.grad <- 0
            cur.hess <- 0
        }

        k <- min(5, m) # this needs to be a parameter at some point
        pset <- unlist(sapply(1:k, function(x) combn(m, x, simplify=F)), recursive=F)
        for (idx_set in pset) {
            # only need genes in the causal configuration using FINEMAP BF trick
            cur.S <- S[idx_set, idx_set]
            cur.Zp <- Zp[idx_set]

            # compute SVD for robust estimation
            # if rank deficient, drop corresponding eigenvectors/values
            res <- svd(cur.S)
            keep <- res$d > tol
            nc <- sum(keep)
            cur.eig <- res$d[keep]
            cur.U <- res$u[,keep]
            cur.chi2 <- prior_chisq / nc
            cur.scaled.Zp <- (t(cur.Zp) %*% cur.U)^2

            # log BF + log prior
            cur.log.BF <- 0.5 * -sum(log(1 + cur.chi2 * cur.eig)) +
                          0.5 *  sum((cur.chi2 / (1 + cur.chi2 * cur.eig)) * cur.scaled.Zp) +
                          nc * log(prb) + (m - nc) * log(1 - prb)

            # keep track for marginal likelihood
            log.marginal <- log_sum_exp(log.marginal, cur.log.BF)

            if (compute.grad) {
                cur.grad <- cur.grad + exp(cur.log.BF) * 0.5 * (-sum(cur.eig / (nc + prior_chisq * cur.eig)) + 
                                                                sum(nc * cur.scaled.Zp / (nc + prior_chisq * cur.eig)^2))
                cur.hess <- cur.hess + exp(cur.log.BF) * (sum(cur.eig^2 / (nc + prior_chisq * cur.eig)^2) -
                                                          sum((nc * cur.eig * cur.scaled.Zp) / (nc + prior_chisq * cur.eig)^3))
            } 

            # marginalize the posterior for marginal-posterior-inclusion probability (pip) on causals
            for (idx in idx_set) {
                if (pip[idx] == 0) {
                    pip[idx] <- cur.log.BF
                } else {
                    pip[idx] <- log_sum_exp(pip[idx], cur.log.BF)
                }
            }
        }
        
        # convert logpips to pips
        pip <- exp(pip - log.marginal)
        null.pip <- exp(null.pip - log.marginal)

        if (compute.grad) {
            cur.grad <- cur.grad / exp(log.marginal)
            cur.hess <- cur.hess / exp(log.marginal)
        }
        if (posterior_check > 0) {
            res <- svd(S)
            eigs <- res$d
            U <- res$u
            sim.mu <- mu
            simulation <- bind_rows(lapply(1:posterior_check,
                                           function(x) {
                                               sim.caus <- rbinom(n=m, size=1, prob=pip)
                                               sim.chi2 <- ifelse(sum(sim.caus) > 0, prior_chisq / sum(sim.caus), 0)
                                               sim.V <- S %*% diag(prior_chisq * sim.caus) %*% S + S
                                               sim.Zs <- rmvn(n=1, mu=sim.mu, sigma=sim.V)
                                               #sim.Zs <- t(rnorm(n=m, mean=sim.mu, sd=sqrt(eigs * prior_chisq * sim.caus + 1))) %*% U
                                               data.frame(ID=cur.ID,
                                                          DATA.Z=cur.Z,
                                                          CAUS=t(t(sim.caus)),
                                                          SIM.Z=t(sim.Zs))
                                           }))
        } else {
            simulation <- NA
        }

    } else {
        marginal <- 0
        inter.z <- NA
        inter.p <- NA
        BF <- dnorm(cur.Z, mean=0, sd=sqrt(prior_chisq + 1)) / dnorm(cur.Z, mean=0, sd=1) * prb
        marginal <-  (1 - prb) + BF
        pip <- BF / marginal
        null.pip <- (1 - prb) / marginal
        log.marginal <- log(marginal)
        if (compute.grad) {
            cur.grad <- pip * (0.5 * (-1/(prior_chisq + 1)) + (cur.Z^2 / (prior_chisq + 1)^2)) + null.pip * 0.5 * (-1 / (prior_chisq + 1)) 
            cur.hess <- pip * ((1 / (prior_chisq + 1)^2) - (cur.Z^2 / (prior_chisq + 1)^3)) + null.pip * (1 / (prior_chisq + 1)^2)
        } 
        Zp <- cur.Z
        if (posterior_check > 0) {
            sim.caus <- rbinom(posterior_check, size=1, prob=pip)
            sim.Z <- rnorm(posterior_check, mean=0, sd=sqrt(prior_chisq * sim.caus + 1))
            simulation <- data.frame(ID=rep(cur.ID, posterior_check),
                                     CAUS=t(t(sim.caus)),
                                     DATA.Z=rep(cur.Z, posterior_check), 
                                     SIM.Z=t(t(sim.Z)))
        } else {
            simulation <- NA
        }
    }

    return (list(RESID.Z=Zp, INTER.Z=inter.z, INTER.P=inter.p, PIP=pip, NULL.PIP=null.pip,
                 MARG.LOG.LIKE=log.marginal, GRAD=cur.grad, HESS=cur.hess, SIM=simulation))
}

allele.qc = function(a1,a2,ref1,ref2) {
	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

option_list = list(
  make_option("--input", action="store", default=NA, type='character',
              help="Path to TWAS test output [required]"),
  make_option("--regions", action="store", default=NA, type='character',
              help="Path to gene blocks [required]"),
  make_option("--sig_regions", action="store", default=NA, type='character',
              help="Path to gene blocks [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--minp_input", action="store", default=1.0, type='double',
              help="Minimum p-value to include TWAS Z-score in analysis [default: %default]"),
  make_option("--use_intercept", action="store_true", default=F, type='logical',
              help="Flag to include intercept for posterior calculation"),
  make_option("--prior_chi2", action="store", default=40, type='double',
              help="Prior chi-sq parameter for causal genes [default: %default]"),
  make_option("--learn_prior_chi2", action="store_true", default=F, type='logical',
              help="Flag to learn the prior chi-square parameter using Empirical Bayes"),
  make_option("--prior_prob", action="store", default=1e-3, type='double',
              help="Prior probability for a gene to be causal [default: %default]"),
  make_option("--learn_prior_prob", action="store_true", default=F, type='logical',
              help="Flag to learn the prior causal probability using Empirical Bayes"),
  make_option("--cs_prob", action="store", default=0.9, type='double',
              help="Probability for credible set computation. Annotates genes with in/not-in flag [default: %default]"),
  make_option("--posterior_check", action="store", default=0, type='integer',
              help="Number of posterior simulations for model checking [default: %default]"),
  make_option("--tol", action="store", default=2.220446e-14, type='double',
              help="Error tolerance to determine for non-zero singular values"),
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
options( digits = 5 )

max_iter <- 10

chr = opt$chr
prior_chisq <- opt$prior_chi2
prb <- opt$prior_prob
tol <- opt$tol

# read all weights in list and compute union of SNPs to analyze
LOG("Loading TWAS summary statistics")
wgtlist = read.table(opt$input, as.is=T, head=T)
wgtlist = wgtlist[ wgtlist$CHR == chr & !is.na(wgtlist$TWAS.P) & wgtlist$TWAS.P <= opt$minp_input, ]

blocks <- read.table(opt$regions, h=T)
blocks <- blocks %>% filter(CHR == chr)
wgtlist$BLOCK <- get_independent(wgtlist$CHR, wgtlist$ID, wgtlist$P0, wgtlist$P1, blocks) 
wgtlist <- wgtlist[order(wgtlist$BLOCK, wgtlist$P0),]

sig_regions <- read.table(opt$sig_regions, h=F)

# load in genotype files by chromosome, restrict to matching SNPs and combine
LOG("Loading reference LD panel")
genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")
genos$bed = scale(genos$bed)


uids <- unique(wgtlist$BLOCK)
B <- length(uids)
N = length(unique(wgtlist$ID))

ord <- order(unlist(lapply(strsplit(uids, ":"), function(x) as.numeric(strsplit(x[2], "\\.\\.")[[1]][1]))))
ids <- uids[ord]

if (opt$learn_prior_prob || opt$learn_prior_chi2) {
    iter <- 1
    like_diff <- 1
    log_like <- 0
    old_log_like <- 0

    # Perform Empirical Bayes to fit the prior probability of causality
    # using the observed Z scores and correlation structure.
    start <- Sys.time()
    LOG("Performing Empirical Bayes estimation of prior probability for causality")
    while (like_diff > 1e-3 && iter < max_iter) {

        old_log_like <- log_like
        log_like <- 0
        sum_prob <- 0
        m_gene <- 0
        cur.GRAD <- 0
        cur.HESS <- 0

        for ( i in ids ) {

            ge.keep <- wgtlist$BLOCK == i

            cur.FILE <- wgtlist$FILE[ge.keep]
            cur.ID <- wgtlist$ID[ge.keep]
            cur.MODEL <- wgtlist$MODEL[ge.keep]
            cur.CHR <- wgtlist$CHR[ge.keep]
            cur.P0 <- wgtlist$P0[ge.keep]
            cur.P1 <- wgtlist$P0[ge.keep]
            cur.Z <- wgtlist$TWAS.Z[ge.keep]
            cur.P <- wgtlist$TWAS.P[ge.keep]

            # region can be skipped
            # theoretically we can still use 'null' regions to help fit prior
            # question is whether we should since data are not similar to 'hit' regions
            if (!i %in% sig_regions$V1) {
                next
            }

            params <- get_local_params(cur.FILE, cur.ID, cur.MODEL, genos)
            cur.LD <- params$LD
            cur.SW <- params$SW
            cur.keep <- params$FLAGS

            res <- fine_map(cur.Z, cur.ID, cur.SW, cur.LD, prb, prior_chisq, intercept=opt$use_intercept, compute.grad=T, tol=opt$tol)

            cur.PIP <- res$PIP
            cur.MARG.LOG.LIKE <- res$MARG.LOG.LIKE
            cur.HESS <- cur.HESS + res$HESS
            cur.GRAD <- cur.GRAD + res$GRAD

            m_gene <- m_gene + length(cur.PIP)
            sum_prob <- sum_prob + sum(cur.PIP)
            log_like <- log_like + cur.MARG.LOG.LIKE
        }
        if (opt$learn_prior_prob) {
            prb <- sum_prob / m_gene
            LOG(paste("Prior prob for iteration", iter, "=", prb))
        }
        if (opt$learn_prior_chi2) {
            prior_chisq <- prior_chisq + cur.GRAD / cur.HESS # take Newton-ascent step => maximize log-like
            LOG(paste("Prior causal-effect var GRAD for iteration", iter, "=", cur.GRAD))
            LOG(paste("Prior causal-effect var HESS for iteration", iter, "=", cur.HESS))
            LOG(paste("Prior causal-effect var estimate for iteration", iter, "=", prior_chisq))
        }

        like_diff <- abs(log_like - old_log_like)
        LOG(paste("Log-likelihood for iteration", iter, "=", log_like))
        if (log_like < old_log_like & iter > 1) {
            LOG("Log-likelihood decreased. Something is screwy", LOG_WARNING)
        }
        iter <- iter + 1
    }
    done <- difftime(Sys.time(), start, units="secs")
    LOG(paste("Completed Empirical Bayes in", done[[1]], "seconds"))
} 

start <- Sys.time()
LOG(paste("Performing fine mapping with prior probability for causality =", prb, "and prior chi2 =", prior_chisq))
out.tbl = NULL
out.sim.tbl = NULL
# Perform fine mapping over regions with at least one TWS hit using the prior determined
# either from EmpBayes or whatever the user supplied
for ( i in ids ) {

	ge.keep <- wgtlist$BLOCK == i

    cur.FILE <- wgtlist$FILE[ge.keep]
    cur.MODEL <- wgtlist$MODEL[ge.keep]
    cur.ID <- wgtlist$ID[ge.keep]
    cur.CHR <- wgtlist$CHR[ge.keep]
    cur.P0 <- wgtlist$P0[ge.keep]
    cur.P1 <- wgtlist$P0[ge.keep]
    cur.Z <- wgtlist$TWAS.Z[ge.keep]
    cur.P <- wgtlist$TWAS.P[ge.keep]

    if (!i %in% sig_regions$V1) {
        LOG(paste("Skipping region", i, "with", length(unique(cur.ID)), "genes"))
        next
    }
    LOG(paste("Finemapping region", i, "with", length(unique(cur.ID)), "genes"))

    params <- get_local_params(cur.FILE, cur.ID, cur.MODEL, genos)
    cur.LD <- params$LD
    cur.SW <- params$SW
    cur.keep <- params$FLAGS

    # ffs just make a data.frame already...
    cur.FILE <- cur.FILE[cur.keep]
    cur.MODEL <- cur.MODEL[cur.keep]
    cur.ID <- cur.ID[cur.keep]
    cur.CHR <- cur.CHR[cur.keep]
    cur.P0 <- cur.P0[cur.keep]
    cur.P1 <- cur.P1[cur.keep]
    cur.Z <- cur.Z[cur.keep]
    cur.P <- cur.P[cur.keep]

    res <- fine_map(cur.Z,
                    cur.ID,
                    cur.SW,
                    cur.LD,
                    prb,
                    prior_chisq,
                    intercept=opt$use_intercept,
                    posterior_check=opt$posterior_check,
                    verbose=T,
                    tol=opt$tol)
    
    cur.INTERCEPT.Z <- res$INTER.Z
    cur.INTERCEPT.P <- res$INTER.P
    cur.RESID.Z <- res$RESID.Z
    cur.PIP <- res$PIP
    cur.NULL.PIP <- res$NULL.PIP

    tbl <- data.frame(FILE=c("NULL", "INTERCEPT", cur.FILE),
                      ID=c(paste0("NULL.", i), paste0("INTERCEPT.", i), as.character(cur.ID)),
                      CHR=chr,
                      P0=c(0, 0, cur.P0),
                      P1=c(0, 0, cur.P1),
                      BLOCK=i,
                      Z=c(0, cur.INTERCEPT.Z, cur.Z),
                      P=c(1, cur.INTERCEPT.P, cur.P),
                      RESID.Z=c(0, 0, cur.RESID.Z),
                      PIP=c(cur.NULL.PIP, 0, cur.PIP))

    if (opt$posterior_check > 0) {
        sim.tbl <- res$SIM %>% mutate(BLOCK = i)
    }
    if (is.null(out.tbl)) {
        out.tbl <- tbl
        out.sim.tbl <- sim.tbl
    } else {
        out.tbl <- rbind(out.tbl, tbl)
        out.sim.tbl <- rbind(out.sim.tbl, sim.tbl)
    }
}
done <- difftime(Sys.time(), start, units="secs")
LOG(paste("Completed fine mapping in", done[[1]], "seconds"))

# remove regions w/o finemapping and add credible set flags
if (!is.null(out.tbl)) {
    out.tbl <- annotate_cred_set(out.tbl, opt$cs_prob)
    write.table(out.tbl, quote=F , row.names=F , sep='\t' , file=opt$out)
}
if (!is.null(out.sim.tbl)) {
    write.table(out.sim.tbl, quote=F, row.names=F, sep='\t', file=paste0(opt$out, ".sim.txt"))
}
