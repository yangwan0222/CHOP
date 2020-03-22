t1_IG <- read.csv("~/Desktop/t1_IG.csv")
x1 <- t1_IG$c_number
df1 <- data.frame(c_number=integer(),
                 prob=double(),
                 stringsAsFactors=FALSE)
for (i in min(x1):max(x1))
{
  df1[nrow(df1) + 1,] = c(i,sum(x1>i)/length(x1))
}
newdf1 <- merge(t1_IG,df1,by = 'c_number')
newdf1 <- newdf1[order(newdf1$probeID),]


t2_IG <- read.csv("~/Desktop/t2_IG.csv")
x2 <- t2_IG$c_number
df2 <- data.frame(c_number=integer(),
                  prob=double(),
                  stringsAsFactors=FALSE)

for (i in min(x2):max(x2))
{
  df2[nrow(df2) + 1,] = c(i,sum(x2>i)/length(x2))
}

newdf2 <- merge(t2_IG,df2,by = 'c_number')
newdf2 <- newdf2[order(newdf2$probeID),]

dp1 <- newdf1$prob*0.5/(newdf1$prob*0.5 + newdf1$part2*0.5)
dp2 <- newdf2$prob*0.5/(newdf2$prob*0.5 + newdf2$part2*0.5)