library(ape)
source("functions.R")
source("binaries.R")

# Simulation 2: Node Age Estimation 
options(scipen=999)
sites <- c(1000, 10000, 50000) 
als <- c(0) 
n_constraints <- 2

sim_dir <- "sim2_altfossil"

# read the simulation number from the command line
j <- commandArgs(trailingOnly=TRUE)[1]

tree <- read.tree("iln.tree")
taxa_in_order <- sort(tree$tip.label)

tree_nobl <- tree
tree_nobl$edge.length <- NULL 
tree_str <- write.tree(tree_nobl)

constraints <- c("constraint root = 1-15;",
                 "constraint outgroup = 1 3-15;",
                 "constraint fossil1 = 4 5;",        
                 "constraint fossil2 = 1 6 7 10;"   
                 )

constraints <- c("constraint root = 1-15;",
                 sprintf("constraint outgroup = %s;", paste0(taxa_in_order[c(1,3:15)], sep="", collapse=" ")),
                 sprintf("constraint fossil1 = %s;", paste0(taxa_in_order[c(4,5)], sep="", collapse=" ")),        
                 sprintf("constraint fossil2 = %s;", paste0(taxa_in_order[c(1,6,7,10)], sep="", collapse=" "))   
                 )

# Scenario 1: narrow
cals1 <- c(
          "calibrate root      = uniform(267, 287);",
          "calibrate outgroup  = uniform(83, 93);  ",
          "calibrate fossil1   = offsetexponential(28.3, 100);",
          "calibrate fossil2   = offsetexponential(51.58, 100);"
)

# Scenario 2: wide
cals2 <- c(
          "calibrate root      = uniform(255.9, 299.8);",
          "calibrate outgroup  = uniform(66.7, 94.3);  ",
          "calibrate fossil1   = offsetexponential(28.3, 100);",
          "calibrate fossil2   = offsetexponential(51.58, 100);"
)

# Scenario 3: wide , wrong 
cals3 <- c(
          "calibrate root      = uniform(255.9, 299.8);",
          "calibrate outgroup  = uniform(66.7, 94.3);  ",
          "calibrate fossil1   = offsetexponential(28.3 , 100);",
          "calibrate fossil2   = offsetexponential(61.58, 100);"
)

# Scenario 4: narrow, wrong 
cals4 <- c(
          "calibrate root      = uniform(267, 287);",
          "calibrate outgroup  = uniform(83, 93);  ",
          "calibrate fossil1   = offsetexponential(28.3 , 100);",
          "calibrate fossil2   = offsetexponential(61.58, 100);"
)


# Scenario 1: narrow
cals5 <- c(
          "calibrate root      = uniform(267, 287);",
          "calibrate outgroup  = offsetlognormal(83, 88.9, 1.78);  ",
          "calibrate fossil1   = offsetexponential(28.3, 100);",
          "calibrate fossil2   = offsetexponential(51.58, 100);"
)

# Scenario 2: wide
cals6 <- c(
          "calibrate root      = uniform(255.9, 299.8);",
          "calibrate outgroup  = offsetlognormal(66.7, 81.3, 5.3);  ",
          "calibrate fossil1   = offsetexponential(28.3, 100);",
          "calibrate fossil2   = offsetexponential(51.58, 100);"
)

# Scenario 3: wide , wrong 
cals7 <- c(
          "calibrate root      = uniform(255.9, 299.8);",
          "calibrate outgroup  = offsetlognormal(66.7, 81.3, 5.3);  ",
          "calibrate fossil1   = offsetexponential(28.3 , 100);",
          "calibrate fossil2   = offsetexponential(61.58, 100);"
)

# Scenario 4: narrow, wrong 
cals8 <- c(
          "calibrate root      = uniform(267, 287);",
          "calibrate outgroup  = offsetlognormal(83, 88.9, 1.78);  ",
          "calibrate fossil1   = offsetexponential(28.3 , 100);",
          "calibrate fossil2   = offsetexponential(61.58, 100);"
)

cals <- list(cals1, cals2, cals3, cals4, cals5, cals6, cals7, cals8)

if (!dir.exists(sprintf("%s/%s", sim_dir, j))) dir.create(sprintf("%s/%s", sim_dir, j), recursive=TRUE)

for (site in sites) {
    aln_file <- sprintf("%s/%s/aln.%s", sim_dir, j, site)
    if (!file.exists(aln_file)) {
        print(sprintf("simulating alignment %s", aln_file))
        aln_str <- sim_dna("iln.tree", site, 15, rep(1,6), al=0, aln_file=aln_file, sort=TRUE)
        writeLines(aln_str, aln_file)
    } else  {
        print(sprintf("reading alingment %s", aln_file))
        aln_str <- readLines(aln_file)
    }

    for (mod in c("pww2", "pww1", "pww0", "full")) {

        for (cal in 5:8) {

            calibrations <- cals[[cal]]

            fname <- sprintf("%s/%s.%s.%s",j,mod,site,cal)

            if (file.exists(sprintf("%s/%s.nex.vstat", sim_dir, fname))) {
                    print(sprintf("skipping %s", fname))
                    next
            }
            print(sprintf("running %s", fname))

            mcmc <- c("propset ParsSPRClock(Tau,V)$prob=0 ExtSPRClock(Tau,V)$prob=0 NNIClock(Tau,V)$prob=0;",
                      "startvals tau=simtree;",
                      "mcmc ngen=100000;") 
                      
            model <- c("lset nst=1;")
            if (mod %in% c("pww2","pww1","pww0")) model <- c(model, "lset usepairwise=yes;")
            if (mod == "pww1") model <- c(model, "lset pwweights=1 nsplits=5;")
            if (mod == "pww2") model <- c(model, "lset pwweights=2 nsplits=5;")

            prior <- c("prset brlenspr=clock:birthdeath;",
                       "prset clockvarpr=iln;",
                       "begin trees;",
                       sprintf("tree simtree= [&R] %s", tree_str),
                       "end;",
                       constraints[1:(2+n_constraints)],
                       calibrations[1:(2+n_constraints)],
                       sprintf("prset topologypr=constraints(root, outgroup, %s);", paste0("fossil", 1:n_constraints, collapse=", ")),
                       "prset clockvarpr=strict;",
                       "prset clockratepr=normal(0.001, 0.01) ;",
                       "prset nodeagepr=calibrated;")
                

            setup_mrb(aln_str, 15, site, prior, model, mcmc, out_dir=sim_dir, fname=sprintf("%s.nex",fname))
            system(sprintf("%s %s/%s.nex > %s/%s.log", pwb_path, sim_dir, fname, sim_dir, fname))
        }
    }
}

