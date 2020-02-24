library(data.table)
library(seqinr)
library(stringr)
library(stats)

p = "IZ"

#path_mydata <- "/home/jane/PhD_Skoltech/Allergy/data_grouped/"
clone_path <- sprintf("/home/Evgeniia.Alekseeva/allergy/clones/%s_ChangeO_vj", p)
path_mydata <- "/home/Evgeniia.Alekseeva/allergy/data/"
path_ig_libraries <- "/home/Evgeniia.Alekseeva/igblast/database/"
path_ig_optional <- "/home/Evgeniia.Alekseeva/igblast/optional_file/"

d_segments <- read.table("/home/Evgeniia.Alekseeva/allergy/data/D_database.fas")
#d_segments <- read.table("/home/jane/PhD_Skoltech/Allergy/D_database.fas")
nd <- nrow(d_segments)
d_df <- cbind(as.character(d_segments[(2*(1:(nd/2))-1),]), as.character(d_segments[(2*(1:(nd/2))),]))
d_df <- data.frame(d_df)
colnames(d_df) <- c("D","seq")

setwd(path_mydata)
patient_files <- sprintf("%s_grouped.txt", p)
timepoints <- fread("timepoints.csv", data.table = FALSE)

#isotype indexing
isotypes_order <- c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2")
## M + D = 1, G1234 = 2, A12 = 3, E = 4 
isotypes_corr <- c(1, 1, 2, 2, 3, 2, 2, 4, 3)
cell_type <- c("Bmem", "PBL",  "PL")
cell_type_code <- c("B", "P", "L")    
borders_table <- c()


df <- fread(patient_files, data.table = FALSE)
df <- data.frame(df)
## remove singletons
singleton_freq <- min(df$cloneFraction)
df <- df[df$cloneFraction > singleton_freq,]
## remove errors of sequencing 
df <- df[df$errorClone == FALSE,]
  
##### by ChangeO
table_of_repeats <- c()
clones <- unique(df$clonalGroupId)
print(sprintf("Number of %s clones:", p))
print(length(clones))

p_b <- length(df$subpop[df$subpop == "Bmem"])
p_p <- length(df$subpop[df$subpop == "PBL"])
p_l <- length(df$subpop[df$subpop == "PL"])

prop_b <- 0
prop_p <- 0
prop_l <- 0

borders_table <- c()
for (cl in clones){
  tryCatch({
    ##################################
    ##################################	
    clone <- df[df$clonalGroupId == cl,]
    size <- length(unique(clone$nSeqVDJRegion))
    print("Clone number:")
    print(which(cl == clones))
    print("Size:")
    print(size)

    if (size > 20){
    
      productive <- c()
      clone_table <- c()
    
      for (seq in unique(clone$nSeqVDJRegion)){
        ########## parameters of the sequence
        ########## if seq repeats
        if (length(clone$nSeqVDJRegion[clone$nSeqVDJRegion == seq]) == 1){
          id <- clone$cloneId[clone$nSeqVDJRegion == seq]
          isot <- clone$bestCHit[clone$nSeqVDJRegion == seq]
          isot <- substring(isot, 1, (nchar(isot) - 3))
          isot <- isotypes_corr[which(isotypes_order == isot)]
          type <- cell_type_code[which(cell_type == clone$subpop[clone$nSeqVDJRegion == seq])]
          date <- timepoints[timepoints$Patient == p,]
          date <- date$Date[date$Timepoint == clone$timepoint[clone$nSeqVDJRegion == seq]]
          replica <- clone$repnum[clone$nSeqVDJRegion == seq]
          seq_name <- paste(id, replica, "isot", isot, type, date, sep = "_")  
        } else {
          rows_of_seq <- which(clone$nSeqVDJRegion == seq)
          names_of_rep <- c()
          for (gg in rows_of_seq){
            id <- as.character(clone$cloneId[gg])
            isot <- as.character(clone$bestCHit[gg])
            isot <- as.character(substring(isot, 1, (nchar(isot) - 3)))
            isot <- as.character(isotypes_corr[which(isotypes_order == isot)])
            type <- as.character(cell_type_code[which(cell_type == clone$subpop[gg])])
            date <- timepoints[timepoints$Patient == p,]
            date <- date$Date[date$Timepoint == clone$timepoint[gg]]
            replica <- as.character(clone$repnum[gg])
            seq_name_ind <- paste(id, replica, "isot", isot, type, date, sep = "_")
            names_of_rep <- c(names_of_rep, seq_name_ind)
          }
          seq_name <- paste(id, "5", "isot", isot, type, date, sep = "_")
          table_of_repeats <- rbind(table_of_repeats,
                                  c(p, cl, seq_name, paste(names_of_rep, collapse = "R"), paste(p, "clone", cl)))
        }
      
      
        write.fasta(as.list(seq), as.list(seq_name), file.out = sprintf("%sfor_blast_%s.fas", clone_path, p))
        command <- sprintf("igblastn -germline_db_V %sIGHV_edit.fas -germline_db_J %sIGHJ_edit.fas -germline_db_D %sIGHD_edit.fas -organism human -query %sfor_blast_%s.fas -auxiliary_data %shuman_gl.aux -outfmt 19 > %sfrom_blast_%s.csv",
                           path_ig_libraries, path_ig_libraries, path_ig_libraries, clone_path, p, path_ig_optional, clone_path, p)
        system(command)
        blast_dt <- fread(sprintf("%sfrom_blast_%s.csv", clone_path, p), data.table = FALSE)


        if (!is.na(blast_dt$stop_codon) & !is.na(blast_dt$productive)){
          if(blast_dt$stop_codon == "FALSE" & blast_dt$productive == "TRUE") {
            productive <- rbind(productive, c(seq_name, seq))
          } 
        }
      }

      #### constructing germline 
      borders <- unlist(str_split(clone$refPoints[1], ":"))
      ## trimming v
      p11 <- as.numeric(borders[11])
      v <- clone$germlineVSeq[1]
      if (p11 < 0) {v <- substring(v, 1, (nchar(v) + p11))}
      
      ## trimming d
      if (nchar(clone$bestDHit[1]) > 0){	
      	d <- str_extract(clone$bestDHit[1], "IGHD[0-9]*-[0-9]*")
      	d <- paste0(">", d, "*01")
      	d <- toupper(as.character(d_df$seq[d_df$D == d]))
      	p14 <- as.numeric(borders[14])
      	p15 <- as.numeric(borders[15])
      	if (nchar(d) > 0 & (nchar(d) + p14) > 0 & (nchar(d) + p15) > 0 & (nchar(d) + p14 + p15) > 0){
        	if (p14 < 0) {
          		d <- substring(d, (1 - p14), nchar(d))
        	} else if (p15 < 0){
          		d <- substring(d, 1, (nchar(d) + p15))
        	}
      	} else {d <- ""}
      } else {d <- ""}

      ## trimming j
      p20 <- as.numeric(borders[18])
      j <- clone$germlineVSeq[1]
      if (p20 < 0) {j <- substring(j, (1 - p20), nchar(j))}
      
      
      
      germ <- paste0(v, d, j)
      germ_name <- paste("germline", clone$bestVHit[1], clone$bestDHit[1], clone$bestJHit[1], sep = "_")  
      
      productive <- rbind(productive, c(germ_name, germ))
      productive <- data.frame(productive)
      colnames(productive) <- c("id","seq")
     
      setwd(clone_path)
      if (nrow(productive) > 20) {

         prop_b <- prop_b + length(clone$subpop[clone$subpop == "Bmem"])/p_b
	       prop_p <- prop_p + length(clone$subpop[clone$subpop == "PBL"])/p_p
         prop_l <- prop_l + length(clone$subpop[clone$subpop == "PL"])/p_l


        clone_name <- paste(p, "clone", cl, "size", nrow(productive), "vj", sep = "_")
        write.fasta(as.list(as.character(productive$seq)), as.list(as.character(productive$id)), file.out = paste(clone_name, "fas", sep = "."))
        print(clone_name)
        borders_table <- rbind(borders_table, c(clone_name, unlist(clone$refPoints[1]))) 
      }
    }
  }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})  
} # for (cl in clones)
	

table_of_repeats <- data.frame(table_of_repeats)
colnames(table_of_repeats) <- c("pat", "cl", "sum_name", "un_names")
saveRDS(table_of_repeats, file = sprintf("table_of_repeats_%s.rds", p))

borders_table <- data.frame(borders_table)
colnames(borders_table) <- c("clone", "ref_points")
saveRDS(borders_table, file = sprintf("borders_table_%s.rds", p))
