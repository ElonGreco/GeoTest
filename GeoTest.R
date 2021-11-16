#---------------------------------#
## step_01 download GES dataset  ##
#---------------------------------#

# 下载数据集
# 下载GSE数据
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84402

#### 通过GEOquery下载GSE数据
# 通过BiocManager安装GEOquery包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

# 加载依赖包
install.packages('R.utils')
library('R.oo')
library('R.methodsS3')
library('R.utils')
library(Biobase)
library(BiocGenerics)

# 加载GEOquery包
library(GEOquery)

eSet <- getGEO("GSE84402", destdir=".", getGPL = F) #下载在当前目录，并且不需要平台信息

# 查看eSet的类型
class(eSet)

# eSet对象里包含着各种各样的信息：表达矩阵，芯片是如何设计的，样本如何分组等。

# 需要从eSet中提取出表达矩阵，再进行后续操作

# 使用list取子集的方法，提取eSet的第一个元素：eSet[[1]];并使用exprs函数把它转化成矩阵
exp <- exprs(eSet[[1]])


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


library(AnnotationDbi)
library(stats4)
library(IRanges)
library(S4Vectors)
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

# 根据第一列所描述的信息我们自己创建分组信息group_list：
# 方法一：使用stringr函数
library(stringr)
# stringr包用于字符串的处理，str_detect是该包里的函数，用来确定一个字符向量能否匹配一种模式。它返回一个与输入向量具有同样长度的逻辑向量：
# str_detect(pd$title,"cancertissue")
# [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE

# 这里的输入向量是数据框pd的第一列pd$title内容，即由28个元素组成的字符型向量。
# str_detect()函数会自动判断cancertissue，是否存在于pd$title向量的每一个元素中，存在返回TRUE，否则返回FALSE。
# str_detect函数处理后我们再使用 ifelse生成符合要求的分组信息group_list
group_list <- ifelse(str_detect(pd$title,"cancertissue")==TRUE,"canceroustissue","non-canceroustissue")
group_list

# 方法二：自己造一个
# 我们已经知道了1，3，5，7...27是cancertissue，2,4,6,8...28是non-canceroustissue，那就自己生成一个符合要求的分组信息：
# group_list <- rep(c(rep("cancertissue",times=1), rep("non-canceroustissue",times=1)) ,times = 14)


#------------------------------#
## step_03 检查表达矩阵 ##
#------------------------------#
















