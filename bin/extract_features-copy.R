#!/usr/bin/env Rscript

library(optparse)
library(ape) # parses gff
library(string) # working with strings
library(Biostrings) # Parse fasta

#### Get arguments ##
option_list <- list(
  make_option(c("-d", "--work_dir"),      type = "character",   default = ".",
              metavar = "character",      help = "Working directory containing input files. Additionally, output will be generated here"),
  make_option(c("-s", "--regard_strands"),  type = "logical",     default = FALSE, action="store_true",
              metavar = "logical",        help = "If true, only siRNAs aligneing to the oppsing strand will be reported"),
  make_option(c("-c", "--get_coding_sequences"),  type = "logical",     default = FALSE, action="store_true",
              metavar = "logical",        help = "If true, coding sequence from a cds file will be extracted"),
  make_option(c("-o", "--group_name"),  type = "character",   default = ".",
              metavar = "character",      help = "Name of group of current organisms (all organisms currently under investigation)"),
  make_option(c("-n", "--matching_nucleotides"),  type = "numeric",   default = "16",
              metavar = "numeric",      help = "Least amount of nucleotides matching to report a hit")
)
opt_parser  <- OptionParser(option_list = option_list)
opt         <- parse_args(opt_parser)

dir_base      <- opt$work_dir
regard_strandiness <- opt$regard_strands
extract_coding_sequences <- opt$get_coding_sequences
organism_name <- opt$organism_name
minimum_matching_nucleotides <- opt$matching_nucleotides


dir_output <- paste(dir_base,"output/",sep='')



 # "halyomorpha" # #  #  # "mamestra" 



if (!dir.exists(dir_output)) {
  dir.create(dir_output)
}

# functions
extract_locus_tag <- function(string,pattern) {
  matches <- str_match(string, pattern)
  if (!is.na(matches[2])) {
    return(matches[2])
  } else {
    return(NA)
  }
}

raw_table <- read.table(file = paste("/home/patrick/PhD/projects/download_data/GFP_vs_", organism_name, ".tsv", sep = ""),
                        header = FALSE,
                        col.names = c("query_id","ref_id","matches","mismatches","gaps",
                                      "strand","alignment_length","query_start","query_end","query_seq",
                                      "ref_start","ref_end","ref_seq","evalue"),
                        sep='\t')

annotation_info <- read.table(file = paste("/home/patrick/PhD/projects/download_data/info_annotation_", organism_name, ".tsv", sep=''),
                                          header = TRUE,
                                          sep='\t')

siRNA_sequences <- readDNAStringSet("/home/patrick/PhD/projects/download_data/GFP_long_siRNAs.fa")
siRNA_names <- names(siRNA_sequences)
sequences <- paste(siRNA_sequences)
siRNA_sequences <- data.frame(seq = sequences)
row.names(siRNA_sequences) <- siRNA_names

annotation_list <- list()
for (row in 1:nrow(annotation_info)) {
  if (annotation_info$annotation[row]) {
    gff_tmp <- read.gff(annotation_info$path[row], na.strings = c("."),GFF3 = TRUE)
    annotation_list[[ annotation_info$accession[row] ]] = gff_tmp
  }
}

# no gaps in siRNAs
edit_table <- raw_table[raw_table$gaps == 0,] 
# minimum of 16 matches
edit_table <- edit_table[edit_table$matches >= minimum_matching_nucleotides,]
# query needs to start at pos 2 since it's the start of the seed sequence
edit_table <- edit_table[edit_table$query_start <= 2,] 
# Check if seed sequence is 100% correct
edit_table <- edit_table[ifelse(edit_table$query_start == 1, substr(edit_table$query_seq,2,8) == substr(edit_table$ref_seq,2,8), substr(edit_table$query_seq,1,7) == substr(edit_table$ref_seq,1,7) ),]
# Turn strand 'plus' to '+' and 'minus' to '-'
edit_table$strand <- ifelse(edit_table$strand == 'plus', '+', '-')

# Start empty data frame to put out results
df <- data.frame(siRNA_name=character(),
                 siRNA_seq=character(),
                 number_hit=numeric(),
                 annotation_available=logical(),
                 target_organism=character(),
                 target_chromosome=character(),
                 posistion_start=numeric(),
                 position_end=numeric(),
                 strand=factor(levels = c("+","-")),
                 target_id=character(),
                 target_rna_type=character(), # maybe convert to factor later on (allow only unique(gff$type))
                 target_biotype=character(),
                 target_parent=character(),
                 product=character(),
                 feature_start=numeric(),
                 feature_end=numeric(),
                 feature_strand=factor(levels = c("+","-")),
                 alignment_length=numeric(),
                 matches=numeric(),
                 mismatches=numeric(),
                 mismatch_locations=numeric())


# Count number of entries 
entry_count <- 1

for (row in 1:nrow(edit_table)) {
  tmp_split <- unlist(str_split(edit_table[row,"ref_id"],"_",n=2))
  edit_table[row,"ref_id"] <- tmp_split[1]
  edit_table[row,"accession"] <- tmp_split[2]
  siRNA_seq <- siRNA_sequences[edit_table[row,"query_id"],"seq"]
  rm(tmp_split)
  if(edit_table[row,"accession"] %in% names(annotation_list)){
    gff <- annotation_list[[edit_table$accession[row]]]
    if(edit_table[row,"ref_id"] %in% unique(gff[,"seqid"]) ){
      # Get all entries of the gff that overlap with the siRNA target site
      tmp_gff <- gff[(gff$start <= edit_table[row,"ref_end"]) & 
                 (gff$end >= edit_table[row,"ref_start"]) &
                 (gff$seqid == edit_table[row,"ref_id"]),]
      # remove region as it describes the whole chromosome
      tmp_gff <- tmp_gff[tmp_gff$type != "region",]
      if(nrow(tmp_gff) > 0){
      
        
        for(entry in 1:nrow(tmp_gff)){
          # Check if strandiness shall be considered, if yes only entries on opposing strands are used (as siRNA need to be complementary)
          if(regard_strandiness & edit_table[row,"strand"] == tmp_gff[entry,"strand"]){
            next
          }
          # Parse necessary attributes
          attributes <- unlist(str_split(tmp_gff$attributes[entry], pattern = ';'))
          feature_id <- attributes[which(startsWith(attributes,'ID='))]
          feature_id <- unlist(str_split(feature_id,pattern='='))[2]
          
          # Check if biotype is available. If not NA is given
          feature_biotype <- attributes[which(startsWith(attributes,'gene_biotype='))]
          if(length(feature_biotype) == 1) {
            feature_biotype <- unlist(str_split(feature_biotype,pattern='='))[2]
          } else {
            feature_biotype <- NA
          }
          
          # Check if product is available. If not NA is given
          feature_product <- attributes[which(startsWith(attributes,'product='))]
          if(length(feature_product) == 1) {
            feature_product <- unlist(str_split(feature_product,pattern='='))[2]
          } else {
            feature_product <- NA
          }
          
          # Check if parent is available. If not NA is given
          feature_parent <- attributes[which(startsWith(attributes,'Parent='))]
          if(length(feature_parent) == 1) {
            feature_parent <- unlist(str_split(feature_parent,pattern='='))[2]
          } else {
            feature_parent <- NA
          }
          df <- rbind(df,data.frame(siRNA_name   = edit_table$query_id[row],
                              siRNA_seq          = siRNA_seq,
                              number_hit         = entry_count,
                              annotation_available=TRUE,
                              target_organism    = edit_table$accession[row],
                              target_chromosome  = edit_table$ref_id[row],
                              posistion_start    = edit_table$ref_start[row],
                              position_end       = edit_table$ref_end[row],
                              strand             = edit_table$strand[row],
                              target_id          = feature_id,
                              target_rna_type    = tmp_gff[entry,"type"],
                              target_biotype     = feature_biotype,
                              target_parent      = feature_parent,
                              product            = feature_product,
                              feature_start      = tmp_gff[entry,"start"],
                              feature_end        = tmp_gff[entry,'end'],
                              feature_strand     = tmp_gff[entry,'strand'],
                              alignment_length   = edit_table$alignment_length[row],
                              matches            = edit_table$matches[row],
                              mismatches         = edit_table$mismatches[row],
                              mismatch_locations = NA)) # TODO: get from both strings
        }
      } else {
        df <- rbind(df,data.frame(siRNA_name   = edit_table$query_id[row],
                                  siRNA_seq          = siRNA_seq,
                                  number_hit         = entry_count,
                                  annotation_available=TRUE,
                                  target_organism    = edit_table$accession[row],
                                  target_chromosome  = edit_table$ref_id[row],
                                  posistion_start    = edit_table$ref_start[row],
                                  position_end       = edit_table$ref_end[row],
                                  strand             = edit_table$strand[row],
                                  target_id          = NA,
                                  target_rna_type    = NA,
                                  target_biotype     = NA,
                                  target_parent      = NA,
                                  product            = NA,
                                  feature_start      = NA,
                                  feature_end        = NA,
                                  feature_strand     = NA,
                                  alignment_length   = edit_table$alignment_length[row],
                                  matches            = edit_table$matches[row],
                                  mismatches         = edit_table$mismatches[row],
                                  mismatch_locations = NA)) # TODO: get from both strings
      }
    }
    rm(gff)
  } else {
    df <- rbind(df,data.frame(siRNA_name   = edit_table$query_id[row],
                              siRNA_seq          = siRNA_seq,
                              number_hit         = entry_count,
                              annotation_available=FALSE,
                              target_organism    = edit_table$accession[row],
                              target_chromosome  = edit_table$ref_id[row],
                              posistion_start    = edit_table$ref_start[row],
                              position_end       = edit_table$ref_end[row],
                              strand             = edit_table$strand[row],
                              target_id          = NA,
                              target_rna_type    = NA,
                              target_biotype     = NA,
                              target_parent      = NA,
                              product            = NA,
                              feature_start      = NA,
                              feature_end        = NA,
                              feature_strand     = NA,
                              alignment_length   = edit_table$alignment_length[row],
                              matches            = edit_table$matches[row],
                              mismatches         = edit_table$mismatches[row],
                              mismatch_locations = NA)) # TODO: get from both strings
  }
  entry_count <- entry_count + 1
}
df_gene <- df[df$target_rna_type=="mRNA" & !is.na(df$target_rna_type=="mRNA"),]

# remove annotations to open up RAM 
rm(annotation_info, annotation_list,siRNA_sequences)

# Sort df by organism names to only load one genome at a time in order to preserve RAM
occurences <- data.frame()
if(extract_coding_sequences){
  for (organism in unique(df$target_organism[df$annotation_available==TRUE])) {
    file_path_cds <- paste("/home/patrick/PhD/projects/download_data/", organism, '/cds_from_genomic.fna', sep='')
    tmp_cds_sequences <- readDNAStringSet(file_path_cds)
    pattern <- "\\[locus_tag=([^\\]]+)\\]"
    
    sequence_names <- names(tmp_cds_sequences)
    locus_tags <- unlist(lapply(sequence_names, function(x) extract_locus_tag(x,pattern)))
    
    # Extract all entries from the current organism which also have an rna_type
    tmp_df <- df[df$target_organism==organism & !is.na(df$target_rna_type),]
    # Extract all entries with gene as rna type
    tmp_df_gene <- tmp_df[tmp_df$target_rna_type == 'gene',]
    collect_indices <- list()
    
    for (row in 1:nrow(tmp_df_gene)) {
      tmp_locus_tag <- gsub("gene-","",tmp_df_gene[row,"target_id"])
      index_cds <- which(locus_tags == tmp_locus_tag)
      collect_indices <- c(collect_indices,index_cds)
    }
    
    tmp_occurences <- cbind(as.data.frame(table(locus_tags[unlist(collect_indices)])), organism=organism)
    #tmp_occurences <- cbind(tmp_occurences, cds_name = unique(sequence_names[unlist(collect_indices)]))
    occurences <- rbind(occurences,tmp_occurences)
    
    collect_indices <- unlist(unique(collect_indices))
    collect_cds <- tmp_cds_sequences[collect_indices]
    writeXStringSet(collect_cds, paste(dir_output, organism, '_hit_cds.fna', sep=''))
    rm(file_path_cds,tmp_cds_sequences,locus_tags)
  }
}


occurences <- occurences[order(occurences$organism,occurences$Freq,decreasing = TRUE),]
colnames(occurences) <- c('locus_tag','number_hits','organism')
occurences <- occurences[,c('organism','locus_tag','number_hits')]
df_gene <- df_gene[order(df_gene$target_organism,df_gene$target_chromosome,df_gene$posistion_start,decreasing = TRUE),]

write.table(df, file=paste(dir_output, 'results_GFP_vs_',organism_name, '.tsv',sep = ""), quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
write.table(df_gene, file=paste(dir_output,'genes_GFP_vs_',organism_name, '.tsv',sep = ""), quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
write.table(occurences, file=paste(dir_output,'amount_hits_GFP_vs_',organism_name, '.tsv',sep = ""), quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)




