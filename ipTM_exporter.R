##### Install packages if missing
using<-function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        install.packages(need)
        lapply(need,require,character.only=TRUE)
    }
}
packages <- c('dplyr', 'tidyr', 'ggplot2', 'jsonlite')
using(packages)

arg <- commandArgs(TRUE)

#if(length(arg) < 3) {
 #   arg[3] <- 1 } #Default to monomer if no n-mer specified

map <- read.csv(arg[2]) 
azip <- arg[1]

filenames <- unzip(azip, list = T)[,1]
macosx_stuff <- grep("__MACOSX", filenames)
folder_level <- 1
if (length(macosx_stuff) > 0){
    filenames <- filenames[-macosx_stuff]
    folder_level <- 2
}
suppressWarnings(
    folders <- do.call("rbind", strsplit(filenames,"/"))[,folder_level]
    )
unique_folders <- unique(folders)
unique_folders <- unique_folders[grep('Pool', ignore.case = T, unique_folders)]

table_of_iptm <- data.frame()
  #for each folder
for (i in seq_along(unique_folders)){
    pool <- gsub('Pool_', '', unique_folders[i], ignore.case = T)
      #get the summary condifence files of top 5 models for each pool
    files_we_want <- filenames[intersect(which(folders == unique_folders[i]),grep("summary_confidences", filenames))]
    #read the iptms in
    list_of_runs <- list()
    list_of_mins <- list()
    list_of_scores <- list()
    for (j in seq_along(files_we_want)){ #read each model
        jsonfile <- jsonlite::fromJSON(unz(azip,files_we_want[j]))
        list_of_runs[[j]] <- jsonfile$chain_pair_iptm
        list_of_mins[[j]] <- jsonfile$chain_pair_pae_min
        chain_ids <- jsonfile$chain_ids #chain ids should be the same for every file, so it's ok to store just the last one
        list_of_scores[[j]] <- jsonfile$ranking_score
    }
      #average of the top 5 models
    iptms <- apply(simplify2array(list_of_runs), c(1, 2), mean)
    paes <- apply(simplify2array(list_of_mins), c(1, 2), mean)
    sds <- apply(simplify2array(list_of_runs), c(1, 2), sd)
    n <- length(chain_ids)
    idx <- expand.grid(k = seq_len(n), l = seq_len(n))
# Matrix indexing with a two-column matrix pulls out [k, l] for each row
    long_df <- data.frame(
        chain_1       = chain_ids[idx$k],
        chain_2       = chain_ids[idx$l],
        pair_iptm     = iptms[cbind(idx$k, idx$l)],
        pair_pae_min  = paes[cbind(idx$k, idx$l)],
        iptm_sd       = sds[cbind(idx$k, idx$l)],
 # chain_1_ptm   = chain_ptm[idx$i],
  #chain_1_iptm  = chain_iptm[idx$i],
#  chain_2_ptm   = chain_ptm[idx$j],
 # chain_2_iptm  = chain_iptm[idx$j],
        stringsAsFactors = FALSE,
        row.names = NULL
    )
#Add prediction score and pool number
    long_df$Pool <- as.numeric(pool)
    long_df$ranking_score <- mean(unlist(list_of_scores))
    table_of_iptm <- rbind(long_df, table_of_iptm)
}
#
table_of_iptm <- table_of_iptm %>%
    left_join(map %>%
              select(Pool, Chain_id, Protein_id, Length) %>%
              rename(chain_1 = Chain_id,
                     Protein_1 = Protein_id,
                     Length_1 = Length),
              by = join_by(Pool, chain_1)) %>%
    left_join(map %>%
              select(Pool, Chain_id, Protein_id, Length) %>%
              rename(chain_2 = Chain_id,
                     Protein_2 = Protein_id,
                     Length_2 = Length),
              by = join_by(Pool, chain_2)) %>%
    mutate(pair_length = Length_1 + Length_2) %>%
    mutate(exp_iptm = -0.036255571 + 0.004470512 * sqrt(pair_length)) %>%
    mutate(corr_iptm = pair_iptm - exp_iptm) %>%
    select(chain_1, chain_2, corr_iptm, Protein_1, Length_1, Protein_2, Length_2, pair_length, Pool, pair_iptm, exp_iptm, iptm_sd, pair_pae_min, ranking_score)

#Plot iptm correction 
c <- -0.036255571
m <- 0.004470512
expected_line <- data.frame(pair_length = min(table_of_iptm$pair_length):max(table_of_iptm$pair_length))
expected_line$exp_iptm <- c + m * sqrt(expected_line$pair_length)
table_of_iptm %>%
    ggplot(aes(x = pair_length, y = pair_iptm)) + geom_point() +
    geom_line(aes(y = exp_iptm), data = expected_line, color = 'red') +
    labs(x = 'Pair length',
         y = 'Pair ipTM') -> plot

#Get name of bait protein for folder name

bait <- table_of_iptm[table_of_iptm$chain_1 == "A", "Protein_1"]
bait <- unique(bait)
results_folder <- paste0(bait, '_ipTMs')
dir.create(results_folder, showWarnings = F)
#Export files
table_of_iptm %>%
    write.csv(file.path(results_folder, 'All_ipTMs.csv'), row.names = F)
table_of_iptm %>%
    filter(chain_1 == "A") %>%
    select(chain_2, corr_iptm, Protein_2, Length_2, pair_length, Pool, pair_iptm, exp_iptm, iptm_sd, pair_pae_min, ranking_score) %>%
    arrange(desc(corr_iptm)) %>%
    write.csv(file.path(results_folder, paste0(bait, '_only.csv')), row.names = F)
ggsave(file.path(results_folder, 'ipTM_correction.png'), height = 7, width = 7)

