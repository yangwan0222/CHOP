if(! "sesame" %in% tolower((.packages()))){
  library("sesame")
  (.packages())
}

if(! "matrixStats" %in% tolower((.packages()))){
  library("matrixStats")
  (.packages())
}

if(! "dplyr" %in% tolower((.packages()))){
  library("dplyr")
  (.packages())
}

#print("Path example: /home/wany/scr1_wany/IDAT_file/GSE129364\n")
p1 <- "/home/wany/scr1_wany/IDAT_file/GSE129364"
p2 <- "/home/wany/scr1_wany/IDAT_file/GSE35069"

probeID_hm450 <- read.table("/home/wany/scr1_wany/probeID_hm450")
seq_A_hm450 <- read.table("/home/wany/scr1_wany/seq_A_hm450")
hm450 <- cbind(probeID_hm450,seq_A_hm450)
setwd("/home/wany/scr1_wany/temp")

ssets_p1 <- lapply(searchIDATprefixes(p1),readIDATpair)
ssets_p2 <- lapply(searchIDATprefixes(p2),readIDATpair)

# use all oobg values to find mean and std, not using median z oobg
for (i in (1:2))
{
  ssets <- ssets_p1
  p3 <- "/home/wany/scr1_wany/IDAT_file/temp/t1_IG.csv"
  if (i == 2) 
  {
    ssets <- ssets_p2
    p3 <- "/home/wany/scr1_wany/IDAT_file/temp/t2_IG.csv"
  }

# calculate the mean of total intensities for each different probe
  a <- IG(ssets[[1]]) 
  for (i in 2:length(ssets))
  {
    b <- IG(ssets[[i]])
    a <- a+b
  }
  g_mean <- a/length(ssets)

# calculate the std of total intensities for each different probe  
  a1 <- (rowSums(IG(ssets[[1]])) - rowSums(g_mean))^2
  for (i in 2:length(ssets))
  {
    b <- rowSums(IG(ssets[[i]])) - rowSums(g_mean)
    a1 <- a1+b^2
  }
  g_std <- (a1/(length(ssets)-1))^0.5
  g_mean <- rowSums(g_mean)

# calculate the z-scores for each different probe
  a2 <- (rowSums(IG(ssets[[1]])) - g_mean) / g_std
  for (i in 2:length(ssets))
  {
    a2 <- cbind(a2, (rowSums(IG(ssets[[i]])) - g_mean) / g_std)
  }
  x1 <- IG(ssets[[1]])[,0]
  x <- cbind(x1, rowMedians(a2))

# lapply
  df_x = as.data.frame(x)
  newdf<- df_x %>% tibble::rownames_to_column("probeID")
  abc <- hm450[-1,]
  colnames(abc) <- c("probeID",'c_number')
  x <- merge(newdf,abc,by="probeID")
  names(x)[names(x) == 'V1'] <- 'z_scores'
  
  write.csv(x,p3,row.names = TRUE)
}


for (i in (1:2))
{
  ssets <- ssets_p1
  p3 <- "/home/wany/scr1_wany/IDAT_file/temp/t1_oobG.csv"
  if (i == 2) 
  {
    ssets <- ssets_p2
    p3 <- "/home/wany/scr1_wany/IDAT_file/temp/t2_oobG.csv"
  }
  
  # calculate the mean of total intensities for each different probe
  a <- oobG(ssets[[1]]) 
  for (i in 2:length(ssets))
  {
    b <- oobG(ssets[[i]])
    a <- a+b
  }
  g_mean <- a/length(ssets)
  
  # calculate the std of total intensities for each different probe  
  a1 <- (rowSums(oobG(ssets[[1]])) - rowSums(g_mean))^2
  for (i in 2:length(ssets))
  {
    b <- rowSums(oobG(ssets[[i]])) - rowSums(g_mean)
    a1 <- a1+b^2
  }
  g_std <- (a1/(length(ssets)-1))^0.5
  g_mean <- rowSums(g_mean)
  
  # calculate the z-scores for each different probe
  a2 <- (rowSums(oobG(ssets[[1]])) - g_mean) / g_std
  for (i in 2:length(ssets))
  {
    a2 <- cbind(a2, (rowSums(oobG(ssets[[i]])) - g_mean) / g_std)
  }
  x1 <- oobG(ssets[[1]])[,0]
  x <- cbind(x1, rowMedians(a2))
  
  # lapply
  df_x = as.data.frame(x)
  newdf<- df_x %>% tibble::rownames_to_column("probeID")
  abc <- hm450[-1,]
  colnames(abc) <- c("probeID",'c_number')
  x <- merge(newdf,abc,by="probeID")
  names(x)[names(x) == 'V1'] <- 'z_scores'
  
  write.csv(x,p3,row.names = TRUE)
}

t1_IG <- read.csv("/home/wany/scr1_wany/IDAT_file/temp/t1_IG.csv")
t1_oobG <- read.csv("/home/wany/scr1_wany/IDAT_file/temp/t1_oobG.csv")

t2_IG <- read.csv("/home/wany/scr1_wany/IDAT_file/temp/t2_IG.csv")
t2_oobG <- read.csv("/home/wany/scr1_wany/IDAT_file/temp/t2_oobG.csv")

part2 <- pnorm(t1_IG$z_scores,mean = mean(t1_oobG$z_scores),sd = sqrt(var(t1_oobG$z_scores)), lower.tail=FALSE)
t1_IG <- cbind(t1_IG,part2)
part2 <- pnorm(t2_IG$z_scores,mean = mean(t2_oobG$z_scores),sd = sqrt(var(t2_oobG$z_scores)), lower.tail=FALSE)
t2_IG <- cbind(t2_IG,part2)


