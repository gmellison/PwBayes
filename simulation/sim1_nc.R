library(ape)
source("functions.R")
source("binaries.R")

# Simulation 1: branch length estimation

# this file takes one argument, simulation number, 
# from the command line. it simulates the data using seq-gen 
# and then writes and runs the mrbayes file

options(scipen=999)
sites <- c(1000, 10000, 50000) 

als <- c(0, 0.5) 
n_constraints <- 2

j <- commandArgs(trailingOnly=TRUE)[1]

simdir <- "sim1_nc"

if (!dir.exists(simdir)) dir.create(simdir)
if (!dir.exists(sprintf("%s/%s", simdir, j))) dir.create(sprintf("%s/%s", simdir, j))

tree_iln_file <- "iln.tree"
tre <- read.tree("iln.tree")
taxa_in_order <- sort(tre$tip.label)

tre_nobl <- unroot(tre)
tre_nobl$edge.length <- NULL 
tree_str <- write.tree(tre_nobl)

nsts <- c(1,6) # 1: JC  6: GTR

for (site in sites) {
    for (al in als) {
        for (nst in nsts) {

            aln_file <- sprintf("%s/%s/aln.%s.%s", simdir, j, al, site)
            if (nst == 6) {
                aln_file <- sprintf("%s/%s/aln.%s.%s.gtr", simdir, j, al, site)
            }

            if (!file.exists(aln_file)) {
                print(sprintf("simulating alignment %s", aln_file))
                if (nst == 1) {
                    aln_str <- sim_dna(tree_iln_file, site, 15, rep(1,6), al=al, aln_file=aln_file,sort=TRUE)
                } else if (nst == 6) {
                    aln_str <- sim_dna(tree_iln_file, site, 15, 1:6/6, al=al, aln_file=aln_file,sort=TRUE)
                }

                writeLines(aln_str, aln_file)
            } else  {
                print(sprintf("reading alignment %s", aln_file))
                aln_str <- readLines(aln_file)
            }

            for (mod in c("pww2", "pww1", "pww0", "full")) {
                   
                fname <- sprintf("%s/%s.%s.%s",j,al,mod,site)
                if (nst==6) fname <- sprintf("%s.gtr",fname)

                if (file.exists(sprintf("%s/%s.nex.vstat", simdir, fname))) {
                        print(sprintf("skipping %s", fname))
                        next
                }
                print(sprintf("running %s", fname))

                mcmc <- c("propset ExtTBR(Tau,V)$prob=0 ExtSPR(Tau,V)$prob=0;",
                          "propset ParsTBR(Tau,V)$prob=0 ParsSPR1(Tau,V)$prob=0;",
                          "propset NNI(Tau,V)$prob=0;",
                          "startvals tau=simtree;",
                          "mcmc ngen=100000;") 
                          
                if (nst == 1) model <- c("lset nst=1;")
                else if (nst == 6) model <- c("lset nst=6;")

                if (mod %in% c("pww1","pww2","pww3","pww0")) model <- c(model, "lset usepairwise=yes;")
                if (mod == "pww1") model <- c(model, "lset pwweights=1 nsplits=20;")
                if (mod == "pww2") model <- c(model, "lset pwweights=2 nsplits=20;")

                if (al > 0) model <- c(model, "lset rates=gamma;")

                prior <- c("begin trees;",
                           sprintf("tree simtree= [&U] %s", tree_str),
                           "end;"
                )

                if (al > 0) prior <- c(prior, sprintf("prset shapepr=fixed(%s);", al))
                if (nst == 6)  prior <- c(prior, sprintf("prset revmatpr=fixed(%s);", paste0(1:6/6,collapse=",")))

                setup_mrb(aln_str, 15, site, prior, model, mcmc, out_dir=simdir, fname=sprintf("%s.nex",fname))
                system(sprintf("%s %s/%s.nex > %s/%s.log", pwb_path, simdir, fname, simdir, fname))
            }
        }
    }
}

