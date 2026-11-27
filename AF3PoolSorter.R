#### Input syntax
## Rscript PooledAF3_dimer.R proteome_file.faa 'Protein_ID' number_of_protein_monomers[numeric]
#### Load libraries
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
packages <- c('dplyr', 'tidyr', 'stringr', 'seqinr', 'jsonlite')
using(packages)

#### Load functions

#Function to distribute a proteome in equal-sized pools
distr_seqs <- function(table, #Data frame with all proteins to sort. Must contain 'Length' column, and must not include query protein
                       n #Max size of pool. Should be 5000 - length of query protein
                       ){
    stopifnot(all(c("Length") %in% names(seqs))) #Check that Length column is present
    seqs <- table %>%
        filter(Length < n)
    if(nrow(seqs) < nrow(table)) {
        warning('Sequences larger than max pool size were removed')}
    
    seqs %>%
        arrange(desc(Length)) %>% head
    no_pools <- ceiling(sum(seqs$Length) / n) #Ideal number of pools
    seqs$Pool <- rep(1:no_pools, nrow(seqs))[1:nrow(seqs)]
    
    #which pools exceed the max pool size
    seqs %>%
        group_by(Pool) %>%
        summarise(PoolSize = sum(Length)) %>%
        arrange(desc(PoolSize)) %>%
        filter(PoolSize > n) -> ex_pools
    #which pools don't exceed max pool size
    seqs %>%
        group_by(Pool) %>%
        summarise(PoolSize = sum(Length)) %>%
        arrange(PoolSize) %>%
        filter(PoolSize <= n)-> non_ex_pools

    message('Sorting ', nrow(seqs), ' sequences into ', no_pools, ' pools')
    new_seqs <- seqs
    l <- 1
                                        #
    while(l < 10) { #Give up after 10 attempts
        if(nrow(ex_pools) == 0) {
            message('Finished reassigning all pools')
            break}
                                        #    
        for(k in 1:nrow(ex_pools)) {
            i <- ex_pools[k,]$Pool #Exceeding pool 
            j <- 1#protein index within pool
                                        #
            while(j < 10) {
                #Check that j protein to add would not exceed the pool size
                if(new_seqs %>%
                   group_by(Pool) %>%
                   mutate(Index = row_number()) %>%
                   filter(Pool == i & Index == j) %>%
                   pull(Length) > (n - non_ex_pools[k,]$PoolSize) ) {
                 #   print(paste('No room for new protein', l))
                    break}
                #Re-assign j protein from k exceeding pool to k non-exceeding pool
                new_seqs <- new_seqs %>%
                    group_by(Pool) %>%
                    mutate(Index = row_number()) %>%
                    mutate(Pool = ifelse(Pool == i & Index == j,
                                         non_ex_pools[k,]$Pool, Pool))
                #Update index
                new_seqs <- new_seqs %>%
                    group_by(Pool) %>%
                    mutate(Index = row_number()) 
                #Update non-exceeding and exceeding pools
                new_seqs %>%
                    group_by(Pool) %>%
                    summarise(PoolSize = sum(Length)) %>%
                    arrange(PoolSize) %>%
                    filter(PoolSize <= n)-> non_ex_pools
                #Stop if current exceeding pool is no longer exceeding
                if(i %in% non_ex_pools$Pool) {
                 #   print(paste('Finished reassigning pool', i, ', Iteration', l))
                    break}
            }
        }
        #Update exceeding pools
        new_seqs %>%
            group_by(Pool) %>%
            summarise(PoolSize = sum(Length)) %>%
            arrange(desc(PoolSize)) %>%
            filter(PoolSize > n) -> ex_pools
        #Update sequence order
        new_seqs <- new_seqs %>%
            arrange(desc(Length))
        l <- l + 1
    }
    
    new_seqs %>%
        filter(Index == 1 & Pool %in% ex_pools$Pool) %>%
        ungroup-> bigseqs
    new_seqs %>%
        ungroup() %>%
        anti_join(bigseqs, by = join_by(Protein, Sequence, Length, Pool, Index)) -> pools
    result <- list('pools'= pools, 'bigseqs' = bigseqs)
    return(result)
}

#Function to check that the resulting pools do not exceed max pool length
check_pools <- function(table, n) {
    table %>%
        group_by(Pool) %>%
        summarise(PoolSize = sum(Length)) %>%
        filter(PoolSize > n)
}

#Clean sequence: remove spaces and change to upper case
clean_seq <- function(s) toupper(gsub("[[:space:]*]", "", s))

#Functions to export alphafold jobs 
make_job <- function(pool) {
    job_name    <- paste("Pool", pool)
    model_seeds <- 1
    df <- final_table %>%
        filter(Pool == pool) 

    all_ids  <- c(query_ids, df$Protein)
    all_seqs <- c(query_seq, df$Sequence)
    no_seqs        <- length(all_seqs)

    sequences <- lapply(seq_len(no_seqs), function(i) {
        list(proteinChain = list(sequence = all_seqs[i], count = 1))
    })

    
    list(
        name       = job_name,
        modelSeeds = list(model_seeds),
        sequences  = sequences,
        dialect    = "alphafoldserver",
        version    = 1
    )
}
#And maps of the proteins in the jobs
make_map <- function(pool) {
     df <- final_table %>%
        filter(Pool == pool) 

     all_ids  <- c(query_ids, df$Protein)
     all_seqs <- c(query_seq, df$Sequence)
     no_seqs       <- length(all_seqs)

     #Chain IDs: A..Z
     chain_ids <- LETTERS
     if (no_seqs > length(chain_ids)) stop("Too many chains for the available chain IDs.")
     chain_ids <- chain_ids[seq_len(no_seqs)]

     total_res <- sum(nchar(all_seqs))
     
     #Save the mapping so you can match chains back to proteins
     chain_map <- data.frame(Chain_id = chain_ids,
                            Protein_id = all_ids,
                            Length = nchar(all_seqs),
                            Pool = pool,
                            stringsAsFactors = FALSE)

     if(total_res > 5000 ) {
         message('Pool', pool, 'too big, ', total_res, ' total residues')
     }
     chain_map
}





#### Script

## separate arguments
args <- commandArgs(TRUE)

### Get integers for copy numbers
### Taken from varhandle::check.numeric
regexp_pattern <- "(^(-|\\+)?\\d+$)|(^(-|\\+)?(\\d*)e(-|\\+)?(\\d+)$)"
isnumeric <- grepl(pattern = regexp_pattern, x = args)
args_numeric <- as.numeric(args[isnumeric])
args <- args[!isnumeric]
file <- args[grep('*.txt|*.fasta|*.faa', args)]
query <- args[! args %in% file]

if (length(file) > 1) stop('Too many sequence files provided')
if(length(file) == 0) stop('Proteome file not found')
#If no copy numbers provided, assume it's 1 for every protein
if (length(args_numeric) == 0) {
    message('Copy number not provided. Assumed 1')
    args_numeric <- rep(1, length(query))
}
if (length(args_numeric) != length(query)) stop('Mismatched number of proteins/copy number')

#Make table with copy numbers
copy_numbers <- data.frame(Protein = query,
                           CopyNumber = args_numeric)

#Read proteome file
x <- seqinr::read.fasta(file, forceDNAtolower = FALSE, as.string = TRUE)
seqs <- data.frame(Protein = unlist(x) %>% names,
                   Sequence = unlist(x),
                   Length = sapply(x, nchar)) %>%
    as_tibble
#Check if all query proteins are found in proteome
if (length(query) != nrow( seqs[seqs$Protein %in% query,]) ) {
    message("Query protein not found in proteome")
    print(query[!query %in% seqs$Protein])
    stop()
}
query <- seqs[seqs$Protein %in% query,] #Find query protein
##Find size of complex
query <- query %>%
    left_join(copy_numbers, by = join_by(Protein)) 
n <- sum(query$CopyNumber * query$Length)
message("Query complex")
print(copy_numbers)
message("Total length: ", n)
n <- 5000 - n  #Define max pool length, minus query complex
seqs <- seqs[! seqs$Protein %in% query$Protein,] #Remove query protein from proteome
if (nrow(filter(seqs, Length > n)) > 0 ) { #Check if any protein doesn't fit in the pools at all
    seqs_too_big <- seqs %>%
        filter(Length > n)
    write.csv(seqs_too_big, 'Seqs_too_big.csv', row.names = F)
    message(nrow(seqs_too_big), ' sequences too big for current pool size')
}
#Remove sequences that are too big to fit in pools
seqs <- seqs %>%
    filter(Length <= n)



#### Pending
### Fix name of seqs_too_big.csv




#Run initial sorting algorithm
results <- distr_seqs(seqs, n)
#Initialize final table
final_table <- results$pools
#Sometimes it fails to keep the pool size because the biggest protein in some pools doesn't fit in any other pool
#Transfer such proteins from final_table to bigseqs to assign to new pools
misfits <- final_table %>%
    filter(Pool %in% check_pools(final_table, n)$Pool) %>%
    group_by(Pool) %>%
    filter(Length == max(Length)) %>%
    ungroup
final_table <- final_table %>%
    anti_join(misfits, by = join_by(Protein, Sequence, Length, Pool, Index))
results$bigseqs <- results$bigseqs %>%
    bind_rows(misfits)
message("Reassigning big proteins to new pools")
#Iterate algorithm on too big sequences until they are all in a pool 
while(nrow(results$bigseqs) != 0) {
    results <- distr_seqs(results$bigseqs, n)
    results$pools <- results$pools %>%
        mutate(Pool = Pool + max(final_table$Pool)) 
    final_table <- final_table %>%
        bind_rows(results$pools)
}

if(nrow(seqs) == nrow(final_table)) {
    message('All sequences processed succesfully')
}
if(nrow(check_pools(final_table, n)) == 0) {
    message('All pools within size')
} else {message('Pools exceed size limit')
    print(check_pools(final_table, n))}
                                       


##### Export to json

## Check sequences

valid_aa  <- "^[ACDEFGHIKLMNPQRSTVWY]+$"

final_table$Sequence <- vapply(final_table$Sequence, clean_seq, character(1))
final_table$valid    <- grepl(valid_aa, final_table$Sequence)
if (any(!final_table$valid)) {
  warning("Skipping invalid sequences: ",
          paste(final_table$Protein[!final_table$valid], collapse = ", "))
}
final_table <- final_table[final_table$valid, ]


## Query_id for file name, query_ids for map and job functions

#Define query sequence
query_id <- paste(query$Protein, collapse = '_')
query_ids <- character()
query_seq <- character()
#
for( line in 1:nrow(query)) {
    query_ids <- append(query_ids, rep(query$Protein[line], query$CopyNumber[line]))
    query_seq <- append(query_seq, rep(query$Sequence[line], query$CopyNumber[line]))
}


##Check query sequence
query_seq <- clean_seq(query_seq)
if (!all(grepl(valid_aa, query_seq))) stop("Query sequence contains non-standard residues.")

#Define file names and folders for json files
mer <- paste(query$CopyNumber, collapse = '_')
no_pools <- length(unique(final_table$Pool)) #numebr of pools
folder_json <- paste0(query_id, "_", mer,"_pools_json")
folder_map <- paste0(query_id, "_", mer,"_pools_map")

dir.create(folder_json, showWarnings = F)
dir.create(folder_map, showWarnings = F)

#Define pools to process
pools <- unique(final_table$Pool)
pools <- pools[order(pools)]




### apply functions to sets of 30 jobs

no_files <- ceiling(length(pools) / 30) #Number of files to make

for(file in seq_len(no_files) ) {
    start <- ((file-1) * 30) + 1
    end <- min( file * 30, length(pools))
    pools_to_process <- pools[start:end]

    out_json    <- paste0("af3_File", stringr::str_pad(file, 2, 'left', '0'),
                          "_Pools_", start, "-", end, ".json")
    map_file    <- paste0("af3_Map", stringr::str_pad(file, 2, 'left', '0'), ".csv")

    #Export json
    jobs <- lapply(pools_to_process, make_job)
    write_json(jobs,
               file.path(folder_json, out_json), auto_unbox = TRUE, pretty = TRUE)
    message("Wrote ", out_json)
    #
    #Export maps
    maps <- lapply(pools_to_process, make_map)
    maps <- bind_rows(maps)
    maps$File <- file
    write.csv(maps, file.path(folder_map, map_file), row.names = FALSE)   
}

message("Wrote ", no_files, ' files with ', length(pools), ' pools from ',
        nrow(final_table), ' proteins in total.')

#### Put together pool maps


files <- list.files(path = folder_map, pattern = '*.csv')

maps <- lapply(file.path(folder_map, files), read.csv)
maps <- bind_rows(maps) %>% arrange(Pool)
write.csv(maps, paste0(folder_map, '.csv'), row.names = F)
unlink(folder_map, recursive = T)
