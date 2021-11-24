#---------------------------------#
## step_01 download GES dataset  ##
#---------------------------------#

# 下载数据集
# 下载GSE数据
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84402

# 下载依赖包 
# curl包  https://github.com/jeroen/curl   https://cran.r-project.org/web/packages/curl/index.html
# sudo apt-get install -y libcurl-dev
# install.packages("curl") 

# openssl包 https://github.com/jeroen/openssl   https://cran.r-project.org/web/packages/openssl/index.html
# sudo apt-get install -y libssl-dev
# install.packages("openssl")

# xml2包  https://github.com/r-lib/xml2   https://cran.r-project.org/web/packages/xml2/index.html 
# sudo apt-get install -y libxml2-dev
# install.packages("xml2") 
# or
# install.packages("devtools")
# devtools::install_github("r-lib/xml2")

# httr包  https://github.com/r-lib/httr  https://cran.r-project.org/web/packages/httr/index.html
# install.packages("httr")
# or
# install.packages("devtools")
# devtools::install_github("r-lib/httr")

install.packages("nlme")
install.packages("spatial")


#### 通过GEOquery下载GSE数据
# 通过BiocManager安装GEOquery包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery") 

# 加载依赖包
install.packages('R.utils')
library('R.methodsS3')
library('R.oo')

library('R.utils')
library(BiocGenerics)
library(Biobase)


# 加载GEOquery包
library(GEOquery)

eSet <- getGEO("GSE84402", destdir=".", getGPL = F) #下载在当前目录，并且不需要平台信息

# 查看eSet的类型
class(eSet)

# eSet对象里包含着各种各样的信息：表达矩阵，芯片是如何设计的，样本如何分组等。

# 需要从eSet中提取出表达矩阵，再进行后续操作

# 使用list取子集的方法，提取eSet的第一个元素：eSet[[1]];并使用exprs函数把它转化成矩阵
exp <- exprs(eSet[[1]])



#  标准化
 


# 归一化：将每个样本的特征值（在转录组中，特征值就是表达量）转换到同一量纲下，
# 把表达量映射到特定的区间内，区间的上下限由表达量的极值决定，这种区间缩放法是归一化的常用方法。
# 标准化：按照表达矩阵中的一个基因在不同样本中的表达量处理数据，
# 每个样本点都能对标准化产生影响，通过求z-score值，
# 转换为标准正态分布，经过处理的数据的均值为0，标准差为1，因此z-score也称为零-均值规范化。

# 取log对表达量的影响
# 原始的raw counts矩阵是一个离散型的变量，离散程度很高。有的基因表达丰度比较高，
# counts数为10000，有些低表达的基因counts可能10，甚至在有些样本中为0。
# 即使经过了RPKM/FPKM等方法抵消了一些测序技术误差的影响，但高低丰度基因的表达量的差距依然很大。
# 如果对表达量去一下log10，发现10000变成了4，10变成了1，这样之前离散程度很大的数据就被集中了。

# 判断GEO芯片数据表达矩阵是否需要log2转换

# ex <- exp
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# 
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exp <- log2(ex)
# print("log2 transform finished")}else{print("log2 transform not needed")}
# 
# View(exp)


# log2 transform
ex <- exp
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exp <- log2(ex) }


#------------------------------#
## step_02 探针id转换为基因id ##
#------------------------------#

# ID转换的第一步，必须要加载特定的R包，需要下在哪个包，根据GPL来定
# 通过eSet，可以看到，GPL号是GPL570

# 一般有三种方法可以得到芯片探针与gene的对应关系。
# 金标准当然是去基因芯片的厂商的官网直接去下载啦！！！
# 一种是直接用bioconductor的包
# 一种是从NCBI里面下载文件来解析好！https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84402

# 用R的bioconductor包来批量得到芯片探针与gene的对应关系！
# 一般重要的芯片在R的bioconductor里面都是有包的，用一个R包可以批量获取有注释信息的芯片平台

# 可以通过BioConductor官网搜索平台信息，下在相应平台的Bioc_package 
#
# 使用BiocManager安装hgu133plus2.db
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")

library(stats4)
library(S4Vectors)
library(IRanges)
library(AnnotationDbi)
library(hgu133plus2.db)


# 通过命令：
# ls("package:hgu133plus2.db")
# 我们可以看到这个包里面有很多数据集，想要得到probe_id和symbol的对应关系要用hgu133plus2SYMBOL数据集，用toTable函数提取数据集里面的信息：

ls("package:hgu133plus2.db")

ids <- toTable(hgu133plus2SYMBOL)

head(ids)
View(ids)

# 查看基因数量 20862个
# unique()函数的作用是：Extract Unique Elements 去除重复的symbol，只提取不同的元素；length()函数统计去重之后还有多少个基因
length(unique(ids$symbol))

# 再查看每个基因对应多少个探针
# table()函数可以生成频数统计表，这里是统计每个基因symbol出现的次数，然后将其表格化；sort()函数将symbol出现的频率从达到小进行排序；tail()取最后6个即出现频率最大的6个
tail(sort(table(ids$symbol)))

# table统计可以看出有多少基因设计了多少探针，从这个数据可以看出，有9818个基因设计了1个探针，有5312个基因设计了2个探针...，也就是说大部分基因只设计了1个或两个探针
table(sort(table(ids$symbol)))

# x %in% y表示 x 的元素在y中吗？然后返回逻辑值。
# rownames(exp)即表达矩阵exp的行名是文章数据中用到的所有探针ID（probe_id）；ids$probe_id是具有对应基因的所有探针。
# 所以返回的TRUE就是文章数据中有对应基因的探针数。
# 可以看到，有11530个探针没有对应基因名；43145个探针有对应的基因名
table(rownames(exp) %in% ids$probe_id)

# 对探针进行过滤，把没有对应基因名的探针过滤掉
# 过滤的本质就是矩阵取子集，如：matrix[2,]意思就是取矩阵matrix的第2行和所有的列。
# 同理，我们这里exp[rownames(exp) %in% ids$probe_id,]就是取具有对应基因的所有探针（行），和所有的列。
exp <- exp[rownames(exp) %in% ids$probe_id,]

View(exp)
dim(exp)

# 然后，我们使用match函数把ids里的探针顺序改一下，使ids里探针顺序和我们表达矩阵的顺序完全一样：
# match()函数返回的是一个位置向量，该向量记录着第一个参数中每个元素在第二个参数中的位置。所以，此时ids里的探针顺序与表达矩阵exp的探针顺序一一对应：
ids <- ids[match(rownames(exp),ids$probe_id),]


# ids里的探针顺序与表达矩阵exp的探针顺序一一对应后，再通过probe_id将表达矩阵exp进行分组，
# 将同一个symbol所对应的多个探针分成不同的组，并对每组探针进行统计：计算每组中每行探针表达量的平均值（也就是每个探针在6个样本中表达量的均值rowMeans(x)），
# 再取平均值最大的那个探针作为该symbol所对应的唯一探针，该组中的其它探针过滤掉，这样每个symbol就对应一个探针了

# 学习by()函数如何完成以上操作的。《R语言实战》这本书上是这样描述的
# 使用by()分组计算描述性统计量，它可以一次返回若干个统计量。格式为：
# by(data, INDICES, FUN)
# 其中data是一个数据框或矩阵；INDICES是一个因子或因子组成的列表，定义了分组；FUN是任意函数。

# --------------------------------------
# 简单一句话理解就是：by()函数就是根据因子将整个data分成几个小的data.frame，然后进行运算处理。
# 同理，我们这里：
# by(exp, ids$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
# 第二个参数ids$symbol定义了分组，将第一参数—exp表达矩阵分成了若干个小矩阵，每个小矩阵里存放着同一个symbol所对应的所有探针。
# 第三个参数是我们自己定义的函数：计算每个小矩阵中每行探针表达量的平均值（也就是每个探针在28个样本中表达量的均值rowMeans(x)），
# 再取平均值最大的那个探针作为该symbol所对应的唯一探针which.max(rowMeans(x))。
# by()函数就可以返回每个分组里的统计结果，即每个symbol所对应的唯一探针IDprobe_id。
#
# a <- matrix(c(1:6),nrow = 3)
# rowMeans(a)
# ----------------------------------

# by()函数在这里发挥的功能就是将表达矩阵exp中的探针分组，同一个symbol所对应的多个探针分到一组，
# 并对每组探针进行统计得到symbol所对应的唯一探针，如果一个symbol对应的有多个探针，那么
# 所以tmp里放着by()函数的统计结果即每个symbol所对应的唯一探针IDprobe_id，
# 用probes = as.character(tmp)将结果变身为纯字符型向量：
tmp <- by(exp,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # 过滤有多个探针的基因

dim(exp)

# 这时，探针ID和基因symbol就一一对应了，将表达矩阵探针ID(即exp表达矩阵的行名rownames(exp))换为基因symbol:
rownames(exp) <- ids[match(rownames(exp),ids$probe_id),2]

#-------
# 此时，我们已经将探针ID转化成基因symbol了。
# 在转换ID中最重要的是根据GPL平台号找到所对应的R注释包，可是如果找不到GPL平台对应的R注释包怎么办呢？
# 答：我们不用GEO号进行下载，而是下载平台信息(GPL)，从平台信息中选择我们想要的列：探针名、基因名....
# GPL里面的信息量特别大，下载特别考验网速。
# gpl <- getGEO('GPL570', destdir = ".")
# colnames(Table(gpl))
# head(Table(gpl)[,c(1,6,7)]) #看gpl对象中哪一列是我们想要的取出来，发现1/6/7列是我们想要的
# write.csv(Table(gpl)[,c(1,6,7)],"GPL570.csv") #把我们想要的部分即探针名对应的基因名....存起来
#------------



# 获取分组信息 group_list
# 分组信息就是告诉我们哪些组是control；哪些组是tumor。
# 使用pData函数获取分组信息—group_list：
# 
# pd <- pData(eSet[[1]])   
# pData()函数可以得到每个样本的描述信息，一般来说数据框的第一列(title列)描述了哪些是control；哪些是treatment。
pd <- pData(eSet[[1]])
View(pd)
save(exp,file = "DEGinput.Rdata")

# 根据第一列所描述的信息我们自己创建分组信息group_list：
# 方法一：使用stringr函数
library(stringr)
# stringr包用于字符串的处理，str_detect是该包里的函数，用来确定一个字符向量能否匹配一种模式。它返回一个与输入向量具有同样长度的逻辑向量：
# str_detect(pd$title,"cancertissue")
# [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE

# 这里的输入向量是数据框pd的第一列pd$title内容，即由28个元素组成的字符型向量。
# str_detect()函数会自动判断cancertissue，是否存在于pd$title向量的每一个元素中，存在返回TRUE，否则返回FALSE。
# str_detect函数处理后我们再使用 ifelse生成符合要求的分组信息group_list
group_list <- ifelse(str_detect(pd$title,"cancertissue")==TRUE,"cancertissue","noncanceroustissue")
group_list

# 方法二：自己造一个
# 我们已经知道了1，3，5，7...27是cancertissue，2,4,6,8...28是non-canceroustissue，那就自己生成一个符合要求的分组信息：
# group_list <- rep(c(rep("cancertissue",times=1), rep("non-canceroustissue",times=1)) ,times = 14)


#------------------------------#
## step_03 检查表达矩阵 ##
#------------------------------#



# 得到表达矩阵就是描述了某个基因在某个样本的表达量。有了这个表达矩阵我们可以做后面的分析，
# 第一步就是确定我们得到的表达矩阵是否正确：
# 
# 查看管家基因的表达量
# 检测分组之间是否有差异：PCA图、热图和hclust图等等


# 3.1   检验常见基因的表达量
# 查看典型管家基因（如：GAPDH、ACTB）的表达量，如果表达量高于正常值，说明我们数据没问题。
# 此时表达矩阵exp_1的行名已经由探针ID转换成基因名了，所以我们使用exp_l['GAPDH',]来提取该基因在所有样品中的表达量。

exp['GAPDH',]
# 我们可以看到我们数据中两个管家基因的表达量都偏高，符合预期。为什么知道它偏高呢？
# 画一个整体样本所有基因的表达量的boxplot：boxplot(exp)
boxplot(exp)

# 发现大部分基因的表达量都在5-6，而GAPDH、ACTB在13-14左右，所以是偏高的。
# 假如，我们发现管家基因表达量特别低，那我们就要思考是不是在提取表达矩阵的时候哪里出了问题：
# 比如ID转换的时候转换错了等等....



# 从图中可以看到两个分组control和treat基本在一条线上，这样的数据说明可以进行后续比较，
# 如果不在一条线上说明有批次效应batch infect，
# 需要用limma包内置函数normalizeBetweenArrays人工校正一下(Normalization)：
library(limma) 
exp = normalizeBetweenArrays(exp)
boxplot(exp)



# 3.2   看表达矩阵的分布图—画图看各个样本的表达量
# 使用ggplot2画各个样本表达量的boxplot图

# 准备画图所需数据exp_L
# https://cran.r-project.org/web/packages/plyr/index.html
library(plyr)
# https://cran.r-project.org/web/packages/reshape2/index.html
library(reshape2)

# exp_L矩阵是这样分布的：每个基因在第一个样本中的value值，
# 每个基因在第二个样本中的value值....以此类推一共有28个样本。

# 难点攻克：如何得到这样的exp_L矩阵呢？？？使用reshape2包
#reshape2包是一套重构和整合数据集的绝妙的万能工具。

# 大致用法就是，需要首先将数据融合（melt），以使每一行都是唯一的标识符-变量组合。
# 然后将数据重塑（cast）为你想要的任何形状。
# 在重铸过程中，你可以使用任何函数对数据进行整合。
# 我们这里只用到这个包里的数据融合（melt）功能。
# 数据集的融合（melt）是将它重构为这样一种格式：每个测量变量（每个基因在每个样本中的表达量）独占一行，
# 行中带有要唯一确定这个测量所需的标识符变量（基因symbol和样本sample）。
# 注意，必须指定要唯一确定每个测量所需的变量（也就是说基因symbol和样本sample必须对应唯一的表达量），
# 而表示测量变量名的变量将由程序为你自动创建（即表达量独占一行后程序会自动创建表达量所对应的symbol和sample）。
# 
# 说成人话就是，以前exp矩阵是一个基因在28个样本中的表达量占一行，
# melt后就会将表达量独占一行。
# 一个表达量的值需要有两个定语才能唯一指定，即这个表达量是哪个样本（sample）中的哪个基因（symbol）的。

head(exp)
exp_L = melt(exp)
head(exp_L)
colnames(exp_L)=c('symbol','sample','value')
head(exp_L)

# 获得分组信息
library(stringr)
group_list = ifelse(str_detect(pd$title,"cancertissue")==TRUE,"cancertissue","noncanceroustissue")
group_list
exp_L$group = rep(group_list,each=nrow(exp))
head(exp_L)

# ggplot2画图 
library(ggplot2)
p = ggplot(exp_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)


##boxplot图精修版
p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
# p=p+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)


# 关于画样本表达量的分布图，除了上面介绍的boxplot，ggplot2还可以画别的，
# 看情况使用就好，不同的图有不同的展现方式但都在展现同一个问题那就是各个样本的表达量，看自己喜欢用就好：
# p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+geom_violin()
# print(p)
# 
# p=ggplot(exp_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
# print(p)





# 检查样本分组信息
# 检查样本分组信息，一般看PCA图，hclust图

# 更改表达矩阵列名
head(exp)
colnames(exp) = paste(group_list,1:28,sep='')
head(exp)

# hclust图
# 定义nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# 聚类
hc=hclust(dist(t(exp)))
par(mar=c(5,5,5,10)) 
# 绘图
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)


# PCA

# https://cran.r-project.org/web/packages/ggfortify/index.html
library(ggfortify)
# 互换行和列，再dim一下
df=as.data.frame(t(exp))
# 不要view df，列太多，软件会卡住；
dim(df)
dim(exp)

exp[1:6,1:6]
df[1:6,1:6]

df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
save(exp,group_list,file = "step2output.Rdata")






#------------------------------#
## step_04 差异分析及可视化 ##
#------------------------------#

# 芯片数据做差异分析最常用的就是limma包
# 使用这个包需要三个数据：
# 
# 表达矩阵(exp)
# 分组矩阵(design)
# 差异比较矩阵(contrast.matrix)
# 下面我们开始准备这三个输入数据：
# 表达矩阵(exp)我们早就得到了，不用再制作了；
# 我们也得到了存放分组信息的向量group_list，用它来制作我们的分组矩阵


# --------------
# 4.1 limma包做差异分析输入数据的准备

# 输入数据—分组矩阵
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = "step2output.Rdata")
dim(exp)
library(limma)
# 做分组矩阵 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp)
design  #得到的分组矩阵



# 输入数据—差异比较矩阵
# contrast.matrix 这个矩阵声明，我们要把treat组和contor组进行差异分析比较：
# -1和1的意思是contorl是用来被比的，treat是来比的即：treat/contorl

# contrast.matrix <- makeContrasts(paste0(c("noncanceroustissue","cancertissue"),collapse = "-"),levels = design)
contrast.matrix <- makeContrasts(paste0(c("cancertissue","noncanceroustissue"),collapse = "-"),levels = design)
contrast.matrix

# 到此，做差异分析所需要的三个矩阵就做好了：表达矩阵(exp)、分组矩阵(design)、差异比较矩阵(contrast.matrix)
# 我们已经制作好了必要的输入数据，下面开始讲如何使用limma包来进行差异分析


# ----------------
# 4.2 limma包做差异分析
# 只有三个步骤：
# 
# lmFit
# eBayes
# topTable

##step1
fit <- lmFit(exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)

head(nrDEG)
save(nrDEG,file = "DEGoutput.Rdata")

# 此时我们就得到差异分析矩阵（nrDEG），重点看logFC和P值：
# 差异分析就是对每个基因都进行检验，检验基因的logFG是多大、平均表达量是多少、p.value是否显著等...

#---------------
# 4.3 差异表达基因的可视化

# 用limma包得到差异分析表达矩阵后作图检查差异基因是否真的很差异

# 画热图
# 选差异最显著的前25个基因画热图，查看差异是否真的很显著
##热图
rm(list = ls())  ## 魔幻操作，一键清空~

options(stringsAsFactors = F)
load(file = "DEGoutput.Rdata")
load(file = "DEGinput.Rdata")

library(pheatmap)
choose_gene = head(rownames(nrDEG),25)
choose_matrix = exp[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix)

# 火山图
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = "DEGoutput.Rdata")
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))

DEG=nrDEG
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)





# 对比文章的结果
# GPC3
tmp<-(which(row.names(DEG)=="PEG10"))
PEG10 <- DEG[tmp,]
PEG10






