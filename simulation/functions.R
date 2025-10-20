
## simulate dna sequence alignment
# tree: tree string in newick format
# n_site: integer number of sites
# n_tax: integer number of taxa in tree
# aln_file: file name for output
sim_dna <- function(tree,n_site,n_tax,rate,
                     aln_file="aln.txt", 
                     al=0, freq=c(1,1,1,1)/4, 
                     missing_rate = 0, sort=TRUE) {

    alignment_str <- NULL 

    if (n_site < 100000) {
       # simulate some data using seqgen  
       sg_call <- ifelse(al > 0, 
                    paste(seqgen_path, " -mGTR ", 
                       " -f ", paste0(freq,collapse=" "),  
                       " -r ", paste0(rate,collapse=" "), 
                       " -l ", n_site, 
                       " -a ", al, 
                       " -g4 ",
                       " -op ", 
                       " -z ", floor(runif(1)*3928109+2817615), 
                       " < ", tree,
                       " > ", sprintf("%s.txt",aln_file), sep=""),
                    paste(seqgen_path, " -mGTR ", 
                       " -f ", paste0(freq,collapse=" "),  
                       " -r ", paste0(rate,collapse=" "), 
                       " -l ", n_site, 
                       " -op", 
                       " -z ", floor(runif(1)*3928109+2817615), 
                       " < ", tree,
                       " > ", sprintf("%s.txt", aln_file), sep=""))
       system(sg_call)

       aln_from_file <- phybase::read.dna.seq(sprintf("%s.txt",aln_file),"phylip")
       taxa_names  <- aln_from_file$name
       alignment_char <- aln_from_file$seq

       if (missing_rate > 0) {
            tot_chars <- nrow(alignment) * ncol(alignment)
            del_sites <- sample(1:tot_chars, round(tot_chars * missing_rate,0))
            alignment[del_sites] <- "-"
            alignment_char[del_sites] <- "-"
        }

        alignment_str <- apply(alignment_char, 1, 
                           function(x) paste0(x,collapse=""))

        alignment_str <- paste(taxa_names, alignment_str, sep="  ")
        if (sort) alignment_str <- sort(alignment_str)

        return(alignment_str)
    }

    n_sim_so_far <- 0
    sim_rounds <- 0
    if (n_site %% 10000 != 0) stop("simulate data in multiples of 10000 please")
    while (n_sim_so_far < n_site) {
       # simulate some data using seqgen  
       sg_call <- ifelse(al > 0, 
                    paste(seqgen_path, " -mGTR ", 
                       " -f ", paste0(freq,collapse=" "),  
                       " -r ", paste0(rate,collapse=" "), 
                       " -l ", 10000, 
                       " -a ", al,
                       " -g4 ",
                       " -z ", floor(runif(1)*3928109+2817615), 
                       " < ", tree,
                       " > ", sprintf("%s.%s.txt",aln_file,sim_rounds), sep=""),
                    paste(seqgen_path, " -mGTR ", 
                       " -f ", paste0(freq,collapse=" "),  
                       " -r ", paste0(rate,collapse=" "), 
                       " -l ", 10000, 
                       " -z ", floor(runif(1)*3928109+2817615), 
                       " < ", tree,
                       " > ", sprintf("%s.%s.txt", aln_file, sim_rounds), sep=""))
       system(sg_call)

       aln_from_file <- phybase::read.dna.seq(sprintf("%s.%s.txt",aln_file,sim_rounds),"phylip")
       taxa_names  <- aln_from_file$name
       alignment_char <- aln_from_file$seq

    if (missing_rate > 0) {
        tot_chars <- nrow(alignment) * ncol(alignment)
        del_sites <- sample(1:tot_chars, round(tot_chars * missing_rate,0))
        alignment[del_sites] <- "-"
        alignment_char[del_sites] <- "-"
    }

    alignment_str_chunk <- apply(alignment_char, 1, 
                            function(x) paste0(x,collapse=""))
    alignment_str_chunk <- paste(taxa_names, alignment_str_chunk, sep="  ")
    if (sort) alignment_str_chunk <- sort(alignment_str_chunk)

    alignment_str <- c(alignment_str, alignment_str_chunk)

    sim_rounds <- sim_rounds + 1
    n_sim_so_far <- n_sim_so_far + 10000

    }

    return(alignment_str)
}


# helper to write mrbayes file to disk
setup_mrb <- function(aln_str,
                      n_taxa, n_site,
                      prior,
                      model,
                      mcmc,
                      out_dir="sim",
                      fname="",
                      datatype="dna",
                      interleave="yes",
                      burninfrac=0.5,
                      write_file=TRUE
                      ) {

    ## now write the alignment:
    mb_lines <- c("#NEXUS\n",

    # write the data block
                "begin data;", 
                sprintf("dimensions ntax=%s nchar=%s;",n_taxa,n_site),
                sprintf("format datatype=%s interleave=%s gap=- missing=?;",datatype, interleave),
                "matrix",
                aln_str,
                "\t;\nend;\n",

    # write the mrbayes block
                 "begin mrbayes;",
                "log start filename=log.txt;",
                 prior,
                 model,
                 mcmc,
                 sprintf("sump outputname=%s/%s burninfrac=%s;",out_dir,fname,burninfrac),
                 sprintf("sumt outputname=%s/%s burninfrac=%s conformat=figtree;",out_dir,fname,burninfrac),
                 "log stop;",
                 "end;"

                 )

    if (write_file) {
        writeLines(mb_lines, sprintf("%s/%s",out_dir,fname))
    }

}

# wrapper
run_mrb <- function(mrb_path, nex_file, out_file) {
        system(sprintf("%s %s &> %s", mrb_path, nex_file,  out_file), intern=FALSE)
}

# functions for getting mrbayes results 
get_parm_ests <- function(sump_fname) {
    parm_ests <- read.csv(sump_fname, skip=1, head=TRUE, sep="\t")
    return(parm_ests)
}

# compute time
get_compute_time <- function(out_fname) {
    # get compute time 
    l <- readLines(out_fname)
    time_line <- l[grepl("Analysis used", l)]
    if (length(time_line) == 0) {
            print(out_fname)
            return(NA)
    }
    #min_sec <- try(as.numeric(stringr::str_extract_all(time_line, "[1-9\\.]+")[[1]]))
    min_sec <- try(regmatches(time_line, gregexpr("[[:digit:]]+(\\.[[:digit:]]+)?", time_line)))
    if (inherits(min_sec, "try_error")) return(NA) 
    else min_sec <- unlist(min_sec)
    mins <- ifelse(length(min_sec) > 1, min_sec[1], 0)
    secs <- ifelse(length(min_sec) > 1, min_sec[2], min_sec)
    sec <- try(as.numeric(mins)*60 + as.numeric(secs))
    return(sec)
}

