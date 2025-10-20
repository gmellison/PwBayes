library(ape)

# this is the base simulation. Just simulate the 15 taxa tree, 
# and estimate the node ages


# n_site = c(1000, 10000, 100000)
# al = 0.1, 1, 0 
# j = 1, ..., 400

options(scipen=999)
source("functions.R")
sites <- c(1000, 10000, 50000) 
#sites <- c(100000) 

als <- c(0, 0.5) 
n_constraints <- 2

j <- commandArgs(trailingOnly=TRUE)[1]

tre <- read.tree("sim15.tree")
taxa_in_order <- sort(tre$tip.label)

tre_nobl <- tre
tre_nobl$edge.length <- NULL 
tree_str <- write.tree(tre_nobl)
tree_str_unr <- write.tree(unroot(tre_nobl))

nsts <- c(1,6)

constraints <- c("constraint root = 1-15;",
                 "constraint outgroup = 1 3-15;",
                 "constraint fossil1 = 4 5;",        # 23,20
                 "constraint fossil2 = 1 6 7 10;",   # 3,26,29,68
                 "constraint fossil3 = 7 10;",       # 29,68
                 "constraint fossil4 = 14 11;",      # 116,86 
                 "constraint fossil5 = 15 13 3 8;",  # 107,16,58,117
                 "constraint fossil6 = 9 12;"        # 59,94
                 )

constraints <- c("constraint root = 1-15;",
                 sprintf("constraint outgroup = %s;", paste0(taxa_in_order[c(1,3:15)], sep="", collapse=" ")),
                 sprintf("constraint fossil1 = %s;", paste0(taxa_in_order[c(4,5)], sep="", collapse=" ")),        # 23,20
                 sprintf("constraint fossil2 = %s;", paste0(taxa_in_order[c(1,6,7,10)], sep="", collapse=" ")),   # 3,26,29,68
                 sprintf("constraint fossil3 = %s;", paste0(taxa_in_order[c(7,10)], sep="", collapse=" ")),       # 29,68
                 sprintf("constraint fossil4 = %s;", paste0(taxa_in_order[c(11,14)], sep="", collapse=" ")),      # 116,86 
                 sprintf("constraint fossil5 = %s;", paste0(taxa_in_order[c(3,8,13,15)], sep="", collapse=" ")),  # 107,16,58,117
                 sprintf("constraint fossil6 = %s;", paste0(taxa_in_order[c(9,12)], sep="", collapse=" "))         # 59,94
                 )


calibrations <- c(
                  "calibrate root      = uniform(267,287);",
                  "calibrate outgroup  = uniform(83,93);  ",
                  "calibrate fossil1   = offsetexponential(28.3 , 100);",
                  "calibrate fossil2   = offsetexponential(51.58, 100);",
                  "calibrate fossil3   = offsetexponential(48.4,  100);",
                  "calibrate fossil4   = offsetexponential(51.58, 100);",    
                  "calibrate fossil5   = offsetexponential(60.5 , 100);",    
                  "calibrate fossil6   = offsetexponential(51.58, 100);"    
)
simdir <- "sim1a"
if (!dir.exists(sprintf("%s/%s", simdir, j))) dir.create(sprintf("%s/%s", simdir, j))

for (site in sites) {
    for (al in 0.5) {
        for (nst in 6) {

            aln_file <- sprintf("%s/%s/aln.%s.%s", simdir, j, al, site)

            if (!file.exists(aln_file)) {
                print(sprintf("simulating alignment %s", aln_file))
                if (nst == 1) {
                    aln_str <- sim_dna("iln.tree", site, 15, rep(1,6), al=al, aln_file=aln_file,sort=TRUE)
                } else if (nst == 6) {
                    aln_str <- sim_dna("iln.tree", site, 15, 1:6/6, al=al, aln_file=aln_file,sort=TRUE)
                }
j
                writeLines(aln_str, aln_file)
            } else  {
                print(sprintf("reading alingment %s", aln_file))
                aln_str <- readLines(aln_file)
            }

            for (mod in c("pww1", "pww0", "full")) {
                    for (clock in c(TRUE, FALSE)) {
                fname <- sprintf("%s/%s.%s.%s%s",j,al,mod,site,ifelse(clock, ".clock",""))
                if (nst==6) fname <- sprintf("%s.gtr",fname)

                if (file.exists(sprintf("%s/%s.nex.vstat", simdir, fname))) {
                        print(sprintf("skipping %s", fname))
                        next
                }
                print(sprintf("running %s", fname))

                if (clock == TRUE) {
                mcmc <- c("propset ParsSPRClock(Tau,V)$prob=0 ExtSPRClock(Tau,V)$prob=0 NNIClock(Tau,V)$prob=0;",
                          "startvals tau=simtree;",
                          "mcmc ngen=100000;")  

                prior <- c("prset brlenspr=clock:birthdeath;",
                           "prset clockvarpr=iln;",
                           "prset clockratepr=normal(0.001, 0.01);",
                           "begin trees;",
                           sprintf("tree simtree= [&R] %s", tree_str),
                           "end;",
                           constraints[1:(2+n_constraints)],
                           calibrations[1:(2+n_constraints)],
                           sprintf("prset topologypr=constraints(root, outgroup, %s);", paste0("fossil", 1:n_constraints, collapse=", ")),
                           ifelse(mod != "full", "mcmcp initsubmod=yes initrunngen=5000 initsamplefreq=25 initrunburnin=1000;", ""),
                           "prset nodeagepr=calibrated;"
                )


                } else {
                mcmc <- c("propset ExtTbr(Tau,V)$prob=0 ExtSPR(Tau,V)$prob=0;" ,
                          "propset ParsTbr(Tau,V)$prob=0 ParsSPR1(Tau,V)$prob=0 NNI(Tau,V)$prob=0;",
                          "startvals tau=simtree;",
                          "mcmc ngen=100000;")  
 
                prior <- c(
                           "begin trees;",
                           sprintf("tree simtree= [&U] %s", tree_str_unr),
                           "end;",
                           ifelse(mod != "full", "mcmcp initsubmod=yes initrunngen=5000 initsamplefreq=25 initrunburnin=1000;", "")
                )

               
                }

                model <- c("lset nst=6;")

                if (mod %in% c("pww1", "pww0")) model <- c(model, "lset usepairwise=yes;")
                if (mod == "pww1") model <- c(model, "lset pwweights=1 nsplits=5;")

                if (al > 0) model <- c(model, "lset rates=gamma;")

                setup_mrb(aln_str, 15, site, prior, model, mcmc, out_dir=simdir, fname=sprintf("%s.nex",fname))
                system(sprintf("../MrBayes/src/mb %s/%s.nex > %s/%s.log", simdir, fname, simdir, fname))
            }
        }
    }
} 
}

