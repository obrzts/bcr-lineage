suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(igraph))
suppressMessages(library(Biostrings))
suppressMessages(library(stringdist))

run_mir <- function(dt, filename, output_path, mir_path, U = 3, u = 8, d = 3, k = 10){
  write.table(dt %>% select(count, freq, cdr3nt, cdr3aa, v, d, j, VEnd, DStart, DEnd, JStart),
              file=paste0(output_path, "mir/input/", filename), quote=F, sep="\t", row.names = F, na = "")
  command = paste0("java -Xmx20G ",
                   "-cp ", mir_path, 
                   " com.milaboratory.mir.scripts.Examples cdr3nt-graph ",
                   "-F VDJtools -S Human -U ", U, 
                   " -u ", u, 
                   " -d ", d, 
                   " -k ", k, 
                   " -G IGH ",
                   "-I ", output_path, "mir/input/", filename,
                   " -O ", output_path, "mir/output/", filename)
  system(command)
  fread(paste0(output_path, "mir/output/", filename))
}

process_mir_out <- function(dt.mo, dt.cdr3, dt.clns){
  dt.mo <- dt.mo %>%
    mutate(from.id = as.character(from.id), to.id = as.character(to.id))
  
  g <- graph_from_data_frame(dt.mo, directed = FALSE)
  cmp <- decompose(g)
  
  clusters.cdr3 <- lapply(1:length(cmp), 
                          function(x) data.table(cdr3.id = V(cmp[[x]])$name, 
                                                 cluster = rep(x, length(V(cmp[[x]]))))) %>%
    rbindlist
  
  clusters.cdr3 <- clusters.cdr3 %>%
    mutate(cdr3.id = as.integer(cdr3.id)) %>%
    merge(dt.cdr3 %>% select(cdr3.id, cdr3nt), all.y = T)
  
  # add cdr3 singletons if any
  clusters.cdr3.2 <- clusters.cdr3[is.na(clusters.cdr3$cluster),]
  if (nrow(clusters.cdr3.2) > 0){
    clusters.cdr3.2 <- clusters.cdr3.2 %>% 
      mutate(cluster = (length(cmp)+1):(length(cmp)+sum(is.na(clusters.cdr3$cluster))))
    clusters.cdr3 <- rbind(clusters.cdr3[!is.na(clusters.cdr3$cluster),],
                           clusters.cdr3.2)
  }
  
  clusters.clns <- dt.clns %>%
    merge(clusters.cdr3 %>% dplyr::rename(nSeqCDR3=cdr3nt))
  
  clusters.clns
}

align_nseq <- function(seq1, seq2, root1, root2){
  data.table(root1 = root1, root2 = root2,
             score = Biostrings::score(pairwiseAlignment(seq1, seq2, type="global")))
}

root_distance <- function(dt.roots, ncores, root_dissim_threshold = 0.28){
  dt.roots$nseq <- as.character(dt.roots$nseq)
  dt.roots$cloneId <- as.character(dt.roots$cloneId)
  
  # find pairs of sequences w/ small Levenshtein distance
  good.pairs <- stringdistmatrix(dt.roots$nseq, dt.roots$nseq,
                                 method = "lv", useName=T,
                                 nthread = ncores) %>%
    melt %>%
    mutate(Seq1 = as.character(Var1),
           Seq2 = as.character(Var2),
           dist = value,
           dissim = dist / nchar(Seq1)) %>%
    filter(dissim < root_dissim_threshold, Seq1 != Seq2) %>%
    merge(dt.roots %>% select(cloneId1 = cloneId, Seq1 = nseq)) %>%
    merge(dt.roots %>% select(cloneId2 = cloneId, Seq2 = nseq))
  
  # align selected pairs
  nseqs <- dt.roots$nseq
  nseqs <- DNAStringSet(nseqs)
  names(nseqs) <- dt.roots$cloneId
  
  # print(paste0(nrow(good.pairs), " small LVD pairs"))
  
  idx <- mapply(c, good.pairs$cloneId1, good.pairs$cloneId2, SIMPLIFY=FALSE)
  idx %>%
    mclapply(function(x) align_nseq(nseqs[[x[1]]], nseqs[[x[2]]], x[1], x[2]), mc.cores = ncores) %>%
    rbindlist -> scores
  scores
}

merge_roots <- function(dt.clusters, ncores, glob_aln_threshold = -40){
  dt.clusters <- dt.clusters %>%
    mutate(nseq = str_sub(targetSequences, vEnd+1, jStart-1),
           cloneId = as.character(cloneId)) # for igraph functions
  
  dt.roots <- dt.clusters %>%
    filter(!is.na(nseq)) %>%
    mutate(mut.total = str_count(mutationsDetailedInVRegion, ",") + str_count(mutationsDetailedInJRegion, ",")) %>%
    group_by(cluster) %>%
    arrange(mut.total) %>%
    summarise(min.mut = cloneId[1], size=n()) %>%
    ungroup()
  
  dt.roots2 <- dt.roots %>%
    mutate(cloneId = as.character(min.mut)) %>%
    merge(dt.clusters %>% select(cloneId, nseq))
  
  print(paste0(nrow(dt.roots), " roots"))
  root_dist = root_distance(dt.roots2, ncores = ncores)
  
  root_dist2 <- root_dist %>%
    filter(score > glob_aln_threshold)
  
  if (nrow(root_dist2) > 0){
    g_root <- graph_from_data_frame(root_dist2, directed = FALSE)
    cmp_root <- decompose(g_root)
    clusters_root <- lapply(1:length(cmp_root), function(x) data.table(root_id = V(cmp_root[[x]])$name, cmp = x)) %>%
      rbindlist
    clusters_root2 <- clusters_root %>%
      merge(dt.roots2 %>% select(root_id = cloneId, size, cluster)) %>%
      group_by(cmp) %>%
      arrange(-size) %>%
      mutate(max_size_cl = cluster[1]) %>%
      ungroup
    
    clusters_merged <- dt.clusters
    for (i in 1:nrow(clusters_root2)){
      cl_old = clusters_root2$cluster[i]
      cl_new = clusters_root2$max_size_cl[i]
      if (cl_old %in% clusters_merged$cluster){
        clusters_merged[clusters_merged$cluster==cl_old,]$cluster = cl_new
      }
    }
    clusters_merged %>% select(-nseq)
  } else {
    dt.clusters %>% select(-nseq)
  }
}

define_clusters_donor <- function(donor_id, dt.clones, dt.cdr3, output_path, mir_path,
                                  merge_by_root, ncores,
                                  dissim_threshold = 0.2,
                                  substs_threshold = 3,
                                  indel_threshold = 3){
  vj <- dt.clones %>%
    filter(donor == donor_id) %>%
    group_by(v, j) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    arrange(-n)
  
  clusters <- data.table()
  
  for (i in 1:nrow(vj)){
    print(paste0("Sample: ", donor_id, ", processing vj pair ", i, " out of ", nrow(vj)))
    
    vt = vj$v[i]
    jt = vj$j[i]
    print(vt)
    print(jt)
    
    .dt.cdr3 <- dt.cdr3[donor==donor_id & v==vt & j==jt]
    .dt.cdr3$cdr3.id = 1:nrow(.dt.cdr3)
    
    if (nrow(.dt.cdr3) > 1){
      filename = paste0(str_remove(donor_id, ".txt"), "_", vt, "_", jt, ".txt")
      filename = str_remove(filename, "/") # for some clumsy gene names
      
      if (file.exists(paste0(path, "mir/output/", filename))){
        print("mir output file exists")
        dt.mo = fread(paste0(path, "mir/output/", filename)) %>%
          mutate(dissim = (substs+indels)/(nchar(from.cdr3nt.aln)-indels)) %>%
          filter(dissim < dissim_threshold, 
                 substs <= substs_threshold, 
                 indels <= indel_threshold)
      } else {
        print(paste0("running mir... (", nrow(.dt.cdr3), " seqs)"))
        dt.mo = run_mir(.dt.cdr3, filename = filename, output_path = output_path, mir_path = mir_path) %>%
          mutate(dissim = (substs+indels)/(nchar(from.cdr3nt.aln)-indels)) %>%
          filter(dissim < dissim_threshold, 
                 substs <= substs_threshold, 
                 indels <= indel_threshold)
      }
      
      if (nrow(dt.mo) > 0){
        .dt.clns <- dt.clones %>%
          filter(donor==donor_id, v==vt, j==jt)
        
        .clusters = process_mir_out(dt.mo, .dt.cdr3, .dt.clns)
        
        # additional cluster merging
        if (merge_by_root){
          if (length(unique(.clusters$cluster)) > 1){
            print("merging clusters...")
            .clusters = merge_roots(.clusters, ncores = ncores)
          }
        }
        
        .clusters <- .clusters %>%
          group_by(cluster) %>%
          mutate(cluster.size = n()) %>%
          ungroup %>%
          filter(cluster.size > 1) %>%
          select(-cdr3.id)
        
      } else {
        print("single cdr3")
        .clusters <- dt.clones %>%
          filter(donor==donor_id, v==vt, j==jt) %>%
          mutate(cluster = 1, cluster.size = n())
      }
      clusters <- rbind(clusters, .clusters)
    }
  }
  
  clusters %>%
    mutate(donor = donor_id)
}

define_clusters <- function(clone_table_path, output_path, mir_path, 
                            ncores, save_mir_files = F, merge_by_root = F){
  dt.clones <- readRDS(clone_table_path)
  
  dt.cdr3 <- dt.clones %>%
    group_by(donor, v, j, cdr3nt=nSeqCDR3, cdr3aa=aaSeqCDR3) %>%
    summarise(count = n()) %>%
    mutate(freq = count/sum(count), d="IGHD",
           VEnd = -1, DStart = -1, DEnd = -1, JStart = -1) %>%
    as.data.table
  
  if (!file.exists(paste0(output_path, "mir"))){
    system(paste0("mkdir ", output_path, "/mir"))
    system(paste0("mkdir ", output_path, "/mir/output"))
    system(paste0("mkdir ", output_path, "/mir/input"))
  }
  
  clusters <- lapply(unique(dt.cdr3$donor), 
                     function(x) define_clusters_donor(donor_id = x, dt.clones, dt.cdr3, 
                                                       output_path = output_path, mir_path = mir_path,
                                                       merge_by_root = merge_by_root, ncores = ncores)) %>%
    rbindlist %>%
    mutate(cluster_id = paste(donor, v, j, cluster, sep = "_")) %>%
    select(-cluster)
  
  saveRDS(clusters, file=paste0(output_path, "clusters.Rds"))
  
  if (!save_mir_files){
    system(paste0("rm -r ", output_path, "/mir"))
  }
}