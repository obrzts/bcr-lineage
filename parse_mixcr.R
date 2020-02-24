suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(parallel))
suppressMessages(library(Biostrings))


def_alleles <- function(subs, segments, 
                        count_threshold = 5, 
                        read_threshold = 0.4,
                        clonotype_threshold = 0.4){
  
  segment.counts <- segments %>%
    group_by(donor, segment.name) %>%
    summarise(segment.count.clonotypes = n(), 
              segment.count.reads = sum(cloneCount)) %>%
    ungroup
  
  alleles <- subs %>%
    group_by(donor, segment.name, pos.nt, from.nt, to.nt) %>%
    summarise(count.clonotypes = n(), 
              count.read = sum(cloneCount)) %>%
    ungroup %>%
    merge(segment.counts) %>%
    mutate(freq.reads = count.read/segment.count.reads,
           freq.clonotypes = count.clonotypes/segment.count.clonotypes,
           is.allele = count.clonotypes >= count_threshold &
             freq.clonotypes > clonotype_threshold &
             freq.reads > read_threshold)
  
  subs %>%
    merge(alleles %>% select(donor, pos.nt, from.nt, to.nt, is.allele, segment.name))
}

transform_refpoints = function(data) {
  splited_refs = as.data.table(str_split_fixed(data$refPoints, ":", n = 22)[,c(10,12,13,16,17,6,7,8,9,19,20)])
  
  colnames(splited_refs) = c("cdr3Start", "vEnd", "dStart", "dEnd", "jStart",
                             "cdr1Start", "fr2Start", "cdr2Start", "fr3Start", "fr4Start", "fr4End")
  
  tmp$cdr3Start = as.integer(tmp$cdr3Start)
  tmp$vEnd = as.integer(tmp$vEnd)
  tmp$dStart = as.integer(tmp$dStart)
  tmp$dEnd = as.integer(tmp$dEnd)
  tmp$jStart = as.integer(tmp$jStart)
  tmp$cdr1Start = as.integer(tmp$cdr1Start)
  tmp$fr2Start = as.integer(tmp$fr2Start)
  tmp$cdr2Start = as.integer(tmp$cdr2Start)
  tmp$fr3Start = as.integer(tmp$fr3Start)
  tmp$fr4Start = as.integer(tmp$fr4Start)
  tmp$fr4End = as.integer(tmp$fr4End)
  
  data = cbind(as.data.table(data), tmp)
  data
}

def_region <- function(x){
  # assign region to mutation by its position
  # x = table row, x[1] = mutation position, x[2:7] - region borders
  regions <- c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
  pos.nt <- as.vector(x[1])
  points <- as.vector(x[2:7])
  p <- which.min(abs(points - pos.nt))
  ifelse(pos.nt < points[p], regions[p], regions[p+1])
}

merge_indels <- function(dt){
  # concatenate insertions\deletions that are next to each other
  
  dt_new = dt[1,]
  for (i in 2:nrow(dt)){
    n = nrow(dt_new)
    
    if (dt$type[i] == "I"){
      f = (dt$cloneId[i] == dt$cloneId[i-1]) & (dt$pos.nt[i] == dt$pos.nt[i-1]) & (dt$type[i] == dt$type[i-1])
    } else {
      f = (dt$cloneId[i] == dt$cloneId[i-1]) & ((dt$pos.nt[i]-1) == dt$pos.nt[i-1]) & (dt$type[i] == dt$type[i-1])
    }
    
    if (f){
      dt_new$seq[n] = paste0(dt_new$seq[n], dt$seq[i])
    } else {
      dt_new <- rbind(dt_new, dt[i,])
    }
  }
  dt_new
}

read_mutations_detailed <- function(dt.clones){
  muts <- dt.clones %>%
    select(sample, donor, cloneId, starts_with("mutationsDetailed")) %>%
    melt(id.vars = c("sample", "donor", "cloneId"))
  muts$mutation <- sapply(muts$value, function(x) str_split(x, ",")[[1]])
  muts <- muts %>%
    unnest %>%
    filter(mutation != "", mutation != "-")
  
  muts <- cbind(muts, str_split(muts$mutation, ":") %>%
                  unlist %>%
                  matrix(ncol=3, byrow=T) %>%
                  as.data.frame %>%
                  rename(ntMut=V1, aaMutInd=V2, aaMutCumul=V3))
  
  segments <- dt.clones %>%
    select(cloneId, sample, donor, v, j, cloneCount) %>%
    melt(id.vars=c("cloneId", "sample", "donor", "cloneCount")) %>%
    rename(segment.name = value, segment = variable) %>%
    mutate(segment = as.character(toupper(segment)))
  
  subs <- muts %>%
    filter(str_detect(ntMut, "S")) %>%
    mutate(from.nt = str_sub(ntMut,2,2),
           to.nt = str_sub(ntMut,-1,-1),
           pos.nt = as.integer(str_sub(ntMut, 3,-2)),
           from.aa = str_sub(aaMutCumul,2,2),
           to.aa = str_sub(aaMutCumul,-1,-1),
           pos.aa = as.integer(str_sub(aaMutCumul, 3,-2)),
           type = ifelse(from.aa == to.aa, "S", "R"),
           segment = str_sub(variable, -7, -7)) %>%
    select(-value, -variable) %>%
    merge(segments) %>%
    def_alleles(segments) %>%
    merge(dt.clones %>% 
            select(sample, cloneId, refPoints)) %>%
    transform_refpoints()
  
  subs$region <- unlist(apply(select(subs, pos.nt, cdr1Start, fr2Start,
                                     cdr2Start, fr3Start, cdr3Start, fr4Start), 1, def_region))
  
  
  indels <- rbind(muts %>%
                    filter(str_detect(ntMut, "I")) %>%
                    mutate(seq = str_sub(ntMut, -1,-1),
                           pos.nt = as.integer(str_sub(ntMut, 2,-2)),
                           type = "I",
                           segment = str_sub(variable, -7, -7)) %>%
                    merge(segments),
                  muts %>%
                    filter(str_detect(ntMut, "D")) %>%
                    mutate(seq = str_sub(ntMut, 2,2),
                           pos.nt = as.integer(str_sub(ntMut, 3)),
                           type = "D",
                           segment = str_sub(variable, -7, -7)) %>%
                    merge(segments)) %>%
    merge_indels %>%
    mutate(len = nchar(seq))
  
  return(list(subs = subs, indels = indels))
}

transform_alignment = function(data){
  data$vMutations <- sapply(data$allVAlignments, function(x) str_match(x, "\\d+\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|([\\w\\d>:]*)")[1,2])
  data$jMutations <- sapply(data$allJAlignments, function(x) str_match(x, "\\d+\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|([\\w\\d>:]*)")[1,2])
  
  data$vAlignTargetStart <- sapply(data$allVAlignments, function(x) str_match(x, "(\\d+)|")[1])
  data$jAlignTargetStart <- sapply(data$allJAlignments, function(x) str_match(x, "(\\d+)|")[1])
  
  data$vAlignTargetEnd <- sapply(data$allVAlignments, function(x) str_match(x, "\\d+\\|(\\d+)\\|")[2])
  data$jAlignTargetEnd <- sapply(data$allJAlignments, function(x) str_match(x, "\\d+\\|(\\d+)\\|")[2])
  
  data$vAlignQueryStart <- sapply(data$allVAlignments, function(x) str_match(x, "\\d+\\|\\d+\\|\\d+\\|(\\d+)\\|")[2])
  data$jAlignQueryStart <- sapply(data$allJAlignments, function(x) str_match(x, "\\d+\\|\\d+\\|\\d+\\|(\\d+)\\|")[2])
  
  data$allVAlignments <- NULL
  data$allJAlignments <- NULL
  data %>%
    mutate_at(c("vAlignTargetStart", "jAlignTargetStart", 
                "vAlignTargetEnd", "jAlignTargetEnd",
                "vAlignQueryStart", "jAlignQueryStart"), as.integer)
}

mutate_seq <- function(seq, idx, char){
  s = str_split(seq, "")[[1]]
  correct <- idx <= length(s)
  idx <- idx[correct]
  char <- char[correct]
  s[idx] <- char
  return(paste0(s, collapse=""))
}

read_substitutions <- function(dt.clones, ncores){
  dt.clones$cloneId_sample <- dt.clones$cloneId # not unique between samples
  dt.clonesId <- 1:nrow(dt.clones)
  
  v.sub <- mclapply(dt.clones$vMutations, function(x) str_match_all(x, "S([ACGT])([0-9]+)([ACTG])")[[1]],  mc.cores=ncores)
  j.sub <- mclapply(dt.clones$jMutations, function(x) str_match_all(x, "S([ACGT])([0-9]+)([ACTG])")[[1]],  mc.cores=ncores)
  
  names(v.sub) <- as.character(dt.clones$cloneId)
  names(j.sub) <- as.character(dt.clones$cloneId)
  
  dt.clones$v.with.subs <- mclapply(1:nrow(dt.clones),
                                    function(x) mutate_seq(dt.clones$v.germline[x], 
                                                           as.integer(v.sub[[x]][,3]) + 1, 
                                                           v.sub[[x]][,4]),
                                    mc.cores=ncores) %>%
    unlist
  dt.clones$j.with.subs <- mclapply(1:nrow(dt.clones), 
                                    function(x) mutate_seq(dt.clones$j.germline[x], 
                                                           as.integer(j.sub[[x]][,3]) + 1, 
                                                           j.sub[[x]][,4]),
                                    mc.cores=ncores) %>%
    unlist
  
  v.sub.list <- dt.clones$cloneId[str_detect(dt.clones$vMutations, "S")]
  v.sub.list = v.sub.list[!is.na(v.sub.list)]
  j.sub.list <- dt.clones$cloneId[str_detect(dt.clones$jMutations, "S")]
  j.sub.list = j.sub.list[!is.na(j.sub.list)]
  
  v.sub.df <- lapply(v.sub.list, function(x) as.data.frame(v.sub[[as.character(x)]]) %>% mutate(cloneId = x)) %>%
    rbindlist() %>%
    rename(shm=V1, from.nt=V2, pos.nt=V3, to.nt=V4) %>%
    merge(dt.clones %>% 
            rename(segment.name=v, segm.germline=v.germline, segm.with.subs=v.with.subs) %>%
            select(-j, -j.germline, -j.with.subs), by="cloneId") %>%
    mutate(pos.nt.target = as.integer(as.character(pos.nt)) + 1, # 1-based for R
           pos.nt = as.integer(as.character(pos.nt)) - vAlignTargetStart + vAlignQueryStart, # 0-based for plots
           segment = "V", shift = vAlignTargetStart)
  
  j.sub.df <- lapply(j.sub.list, function(x) as.data.frame(j.sub[[as.character(x)]]) %>% mutate(cloneId = x)) %>%
    rbindlist() %>%
    rename(shm=V1, from.nt=V2, pos.nt=V3, to.nt=V4) %>%
    merge(dt.clones %>% 
            rename(segment.name=j, segm.germline=j.germline, segm.with.subs=j.with.subs) %>%
            select(-v, -v.germline, -v.with.subs), by="cloneId") %>%
    mutate(pos.nt.target = as.integer(as.character(pos.nt)) + 1, # 1-based for R
           pos.nt = as.integer(as.character(pos.nt)) - jAlignTargetStart + jAlignQueryStart, # 0-based for plots
           segment = "J", shift = -(jAlignTargetStart+1)%%3,
           shift = replace(shift, shift == 0, 1),
           shift = replace(shift, shift == -1, 0),
           shift = replace(shift, shift == -2, 2))
  
  sub <- rbind(v.sub.df, j.sub.df) %>%
    mutate(pos.nt.target.shift = pos.nt.target - shift,
           segm.germline = str_sub(segm.germline, shift + 1),
           segm.with.subs = str_sub(segm.with.subs, shift + 1),
           pos.aa.target = ifelse(pos.nt.target.shift %% 3 == 0, 
                                  (pos.nt.target.shift %/% 3), 
                                  pos.nt.target.shift %/% 3 + 1), # 1-based
           from.codon = str_sub(segm.germline, (pos.aa.target-1)*3 + 1, pos.aa.target*3),
           to.codon = str_sub(segm.with.subs, (pos.aa.target-1)*3 + 1, pos.aa.target*3),
           pos.aa = ifelse((pos.nt+1) %% 3 == 0, ((pos.nt+1) %/% 3)-1, (pos.nt+1) %/% 3)) # 0-based
  sub <- sub %>%
    filter(!is.na(from.codon), !is.na(to.codon))
  
  codons <- data.frame(codon = unique(c(sub$from.codon, sub$to.codon))) %>%
    mutate(codon = as.character(codon)) %>%
    filter(nchar(codon) == 3)
  codons$aa = as.character(translate(DNAStringSet(codons$codon)))
  
  sub <- sub %>% 
    merge(codons %>% select(from.aa = aa, from.codon = codon)) %>%
    merge(codons %>% select(to.aa = aa, to.codon = codon))
  sub$type = ifelse(sub$from.aa == sub$to.aa, "S", "R")
  sub <- sub %>%
    mutate(all_na = is.na(FR1) + is.na(FR2) + is.na(FR3) + is.na(FR4) + is.na(CDR1) + is.na(CDR2) + is.na(CDR3)) %>%
    filter(all_na < 7)
  
  sub$region <- unlist(apply(select(sub, pos.nt, cdr1Start, fr2Start, cdr2Start, fr3Start, cdr3Start, fr4Start), 1, def_region))
  sub <- sub %>% 
    select(cloneId, segment, segment.name, region, 
           pos.nt.germ=pos.nt.target, pos.nt, from.nt, to.nt, pos.aa, from.aa, to.aa, type) %>%
    merge(dt.clones %>% select(cloneId, cloneId_sample, sample), by="cloneId") %>%
    select(-cloneId) %>%
    rename(cloneId = cloneId_sample)
}

read_indels <- function(dt.clones, ncores){
  dt.clones$cloneId_sample <- dt.clones$cloneId # not unique between samples
  dt.clonesId <- 1:nrow(dt.clones)
  
  v.ind <- mclapply(dt.clones$vMutations, function(x) str_match_all(x, "([ID])([ACGT]*)([0-9]+)([ACTG]*)")[[1]], mc.cores=ncores)
  j.ind <- mclapply(dt.clones$jMutations, function(x) str_match_all(x, "([ID])([ACGT]*)([0-9]+)([ACTG]*)")[[1]], mc.cores=ncores)
  
  names(v.ind) <- as.character(dt.clones$cloneId)
  names(j.ind) <- as.character(dt.clones$cloneId)
  
  v.ind.list <- dt.clones$cloneId[str_detect(dt.clones$vMutations, "[ID]")]
  v.ind.list = v.ind.list[!is.na(v.ind.list)]
  j.ind.list <- dt.clones$cloneId[str_detect(dt.clones$jMutations, "[ID]")]
  j.ind.list = j.ind.list[!is.na(j.ind.list)]
  
  v.ind.df <- mclapply(v.ind.list, 
                       function(x) as.data.frame(v.ind[[as.character(x)]]) %>% mutate(cloneId = x),
                       mc.cores=ncores) %>%
    rbindlist() %>%
    rename(shm=V1, type=V2, from.nt=V3, pos.nt=V4, to.nt=V5) %>%
    merge(dt.clones %>% rename(segment.name=v) %>% select(-j), by="cloneId") %>%
    mutate(pos.nt.germ = as.integer(as.character(pos.nt)) + 1, #1-based 
           pos.nt = as.integer(as.character(pos.nt)) - vAlignTargetStart + vAlignQueryStart, # 0-based
           segment = "V")
  
  j.ind.df <- mclapply(j.ind.list, 
                       function(x) as.data.frame(j.ind[[as.character(x)]]) %>% mutate(cloneId = x),
                       mc.cores=ncores) %>%
    rbindlist()
  
  if (nrow(j.ind.df) > 0){
    j.ind.df <- j.ind.df %>%
      rename(shm=V1, type=V2, from.nt=V3, pos.nt=V4, to.nt=V5) %>%
      merge(dt.clones %>% rename(segment.name=j) %>% select(-v), by="cloneId") %>%
      mutate(pos.nt.germ = as.integer(as.character(pos.nt)) + 1, # 1-based 
             pos.nt = as.integer(as.character(pos.nt)) - jAlignTargetStart + jAlignQueryStart, # 0-based
             segment = "J")
    indels <- rbind(v.ind.df, j.ind.df)
  } else {
    indels <- v.ind.df
  }
  
  indels <- indels %>%
    mutate(from.nt = as.character(from.nt), to.nt = as.character(to.nt), type = as.character(type))
  indels <- rbind(indels %>% filter(type == "I") %>% mutate(seq = to.nt),
                  indels %>% filter(type == "D") %>% mutate(seq = from.nt)) %>%
    filter(!is.na(pos.nt))
  indels <- merge_indels(indels)
  
  indels <- indels %>%
    mutate(all_na = is.na(FR1) + is.na(FR2) + is.na(FR3) + is.na(FR4) + is.na(CDR1) + is.na(CDR2) + is.na(CDR3)) %>%
    filter(all_na < 7)
  indels$region <- unlist(apply(select(indels, pos.nt, cdr1Start, fr2Start, cdr2Start, fr3Start, cdr3Start, fr4Start), 1, def_region))
  
  indels <- rbind(indels %>% filter(type == "I") %>% mutate(len.nt = nchar(seq)),
                  indels %>% filter(type == "D") %>% mutate(len.nt = -nchar(seq)))
  indels <- indels %>% 
    select(cloneId, type, segment, segment.name, region, pos.nt, pos.nt.germ, seq.nt=seq, len.nt) %>%
    merge(dt.clones %>% select(cloneId, cloneId_sample, sample), by="cloneId") %>%
    select(-cloneId) %>%
    rename(cloneId = cloneId_sample)
}

read_mutations <- function(dt.clones, germline_v_path, germline_j_path, ncores){
  germline.v <- fread(germline_v_path) %>%
    mutate(segment.name = str_split_fixed(name, "\\*", 2)[,1],
           segment.nt = paste0(UTR5, L1, L2, VRegion), # change if not aligned to VTranscript
           segment = "V")
  germline.j <- fread(germline_j_path) %>%
    mutate(segment.name = str_split_fixed(name, "\\*", 2)[,1], 
           segment.nt = paste0(JRegion),
           segment = "J")
  
  dt.clones <- dt.clones %>%
    transform_refpoints %>%
    transform_alignment %>%
    merge(germline.v %>% select(v = segment.name, v.germline = segment.nt), by="v") %>%
    merge(germline.j %>% select(j = segment.name, j.germline = segment.nt), by="j") %>%
    mutate(v.germline = as.character(v.germline), j.germline = as.character(j.germline)) %>%
    mutate(FR1 = as.integer(cdr1Start - 1),
           CDR1 = fr2Start - cdr1Start,
           FR2 = cdr2Start - fr2Start,
           CDR2 = fr3Start - cdr2Start,
           FR3 = cdr3Start - fr3Start,
           CDR3 = fr4Start - cdr3Start,
           FR4 = as.integer(nchar(targetSequences) - fr4Start + 1))
  
  return(list(subs = read_substitutions(dt.clones, ncores = ncores),
              indels = read_indels(dt.clones, ncores = ncores)))
}

read_mixcr <- function(mixcr_out_folder, germline_v_path = "", germline_j_path = "", ncores = 4){
  files = list.files(mixcr_out_folder, full.name=T)
  
  clones <- lapply(files, function(i) fread(i) %>% mutate(sample = i)) %>%
    rbindlist() %>%
    mutate(v = str_split_fixed(allVHitsWithScore, "\\*", 2)[,1],
           j = str_split_fixed(allJHitsWithScore, "\\*", 2)[,1], 
           donor = sample) %>% # ! should be changed if any donor is represented by multiple samples
    ungroup
  
  if (sum(str_detect(colnames(clones), "mutationsDetailed")) > 0){
    muts <- read_mutations_detailed(clones)
  } else {
    muts <- read_mutations(clones, germline_v_path, germline_j_path, ncores=ncores)
  }
  
  return(list(clones = clones, subs = muts[["subs"]], indels = muts[["indels"]]))
}
