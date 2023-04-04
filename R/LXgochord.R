
#-------------------------------------------------------------------

LXgochord <- function(gene_file,species,term_number){

#list all the packages that have been installed
all_packages <- data.frame(installed.packages())

#To judge whether a package was installed. If not, it will be installed.
pack <- data.frame(c("devtools","BiocManager","ggnewscale","R.utils",
                     "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych",
                     "ggplot2","ggrepel","RColorBrewer", "ggthemes","rticles",
                     "grid","patchwork","Hmisc","pak","GOplot","ggplotify") )
colnames(pack) <- c("pakages_name")

# To judge whether a package was included in the all_packages: %in%
pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){
  if(pack[i,2]==FALSE)
    install.packages(pack[i,1],update = F,ask = F)
}
rm(i)

# 批量library
packages <- as.character(pack[,1])

for(i in packages){
  library(i, character.only = T)
}
rm(i)


#-----------------

if("tidyverse" %in% all_packages$Package==FALSE)
  pak::pak("tidyverse/tidyverse")
library(tidyverse)


#-----------------
BiocManager_pack <- data.frame(c("DOSE","clusterProfiler","do","enrichplot",
                                 "pathview","BiocParallel","GO.db","KEGGREST",
                                 "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"))
# human: "org.Hs.eg.db"
# mouse: "org.Mm.eg.db"
# rat: "org.Rn.eg.db"

BiocManager_pack$type <- BiocManager_pack[,1] %in% all_packages$Package

for (i in 1:nrow(BiocManager_pack)){
  if(BiocManager_pack[i,2]==FALSE)
    BiocManager::install(BiocManager_pack[i,1],update = F,ask = F)
}

# 批量library
Bio_packages <- as.character(BiocManager_pack[,1])
for(i in Bio_packages){
  library(i, character.only = T)
}
rm(i)


#---------------------------------------------------------------
if(dir.exists("analysis result"))
   unlink(x = "analysis result",recursive = T) # 删除文件夹

 dir.create("analysis result") # 新建文件夹

gene_df <- read.xlsx(gene_file) # 读取数据
gene_df <- na.omit(gene_df)     # 去除NA项
gene_df[,2] <- as.numeric(gene_df[,2]) # 转换成数据类型
gene_df[,3] <- as.numeric(gene_df[,3])  # 转换成数据类型

table(is.infinite(gene_df[,2])) #判断是否存在无穷值Inf
gene_df <- dplyr::filter(gene_df,!gene_df[,2]==Inf) #去除Inf项
write.xlsx(gene_df,"analysis result/all_genes.xlsx")

if(max(gene_df[,2])>0){
      up_gene <- dplyr::filter(gene_df,gene_df[,2]>0)
      write.xlsx(up_gene,"analysis result/up_genes.xlsx")
     }

if((min(gene_df[,2])<0)){
        down_gene <- dplyr::filter(gene_df,gene_df[,2]<0)
        write.xlsx(down_gene,"analysis result/down_genes.xlsx")
      }

#--------------------- GOchord --------------------------------
df <- gene_df

log <- colnames(df)[2] %>% tolower() # 转换成小写字母

log <- gsub("\\s", "", log) # 去掉所有的空格，包括前、后和中间的空格
# log <- gsub("[[:space:]]", "", log) # 去掉所有的空格，包括前后和中间的空格

if (str_sub(log,1,4)=="log2"){
    colnames(df)[2] <- "log2FC"
    df$log2FC <- as.numeric(df$log2FC)
    df$log2FC <- 2^(df$log2FC)
    df$log2FC <- log(df$log2FC)
    colnames(df)[2] <- "logFC"
   }

colnames(df) <- c("gene_symbol","logFC","pvalue")
table(duplicated(df$gene_symbol))

df <- distinct(.data = df,gene_symbol,.keep_all = T)

df <- dplyr::filter(df,abs(logFC)>0.30103) # 选出 FC>2的项，(log10(2)=0.30103
df <- dplyr::filter(df,pvalue<0.05) # 筛选p<0.05

spe <- c("human","mouse","rat")
species_input <- tolower(species)  # 转换成小写字母
species_input <- gsub("\\s", "", species_input) # 去掉所有的空格，包括前、后和中间的空格

if(!species_input %in% spe)
{text <- paste0(species," was not the species of human, rat, or mouse.Please check it.")
stop(text)
}

# hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
gene_human <- df
colnames(gene_human)[1] <- "gene_symbol"
gene_human$gene_symbol <- toupper(gene_human$gene_symbol)

# mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
gene_animal <- df
colnames(gene_animal)[1] <- "gene_symbol"

#------------enrichGO分析--------------------------------------

if(species_input== "human"){
  go <- enrichGO(gene=gene_human$gene_symbol,
                OrgDb=org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "ALL",  #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", #one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                # universe,    #background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
                # qvalueCutoff = 0.2,
                minGSSize = 2, # minimal size of genes annotated by Ontology term for testing.
                maxGSSize = 2000, # maximal size of genes annotated for testing
                readable = T, #GO富集的结果的gene以gene symbol形式出现。whether mapping gene ID to gene Name
                pool = FALSE  # If ont='ALL', whether pool 3 GO sub-ontologies
                )

  go_result <- data.frame(go)
  }


if(species_input== "mouse"){
  go <- enrichGO(gene=gene_animal$gene_symbol,
                 OrgDb=org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont = "ALL",  #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", #one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                 # universe,    #background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
                 # qvalueCutoff = 0.2,
                 minGSSize = 2, # minimal size of genes annotated by Ontology term for testing.
                 maxGSSize = 2000, # maximal size of genes annotated for testing
                 readable = T, #GO富集的结果的gene以gene symbol形式出现。whether mapping gene ID to gene Name
                 pool = FALSE  # If ont='ALL', whether pool 3 GO sub-ontologies
                )

  go_result <- data.frame(go)
  }


if(species_input== "rat"){
  go <- enrichGO(gene=gene_animal$gene_symbol,
                 OrgDb=org.Rn.eg.db,
                 keyType = "SYMBOL",
                 ont = "ALL",  #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", #one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                 # universe,    #background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
                 # qvalueCutoff = 0.2,
                 minGSSize = 2, # minimal size of genes annotated by Ontology term for testing.
                 maxGSSize = 2000, # maximal size of genes annotated for testing
                 readable = T, #GO富集的结果的gene以gene symbol形式出现。whether mapping gene ID to gene Name
                 pool = FALSE  # If ont='ALL', whether pool 3 GO sub-ontologies
                 )

  go_result <- data.frame(go)
  }

go_result$geneID <- str_replace_all(go_result$geneID,"/",",") #把/换成,

write.xlsx(go_result,"analysis result/GO_analysis_result.xlsx",rowNames=TRUE)

process <- go_result[,c(1,3,10)]
colnames(process) <- c("Category","Term","gene_count")
process$term_number <- paste0("T",c(1:nrow(process)))

write.xlsx(process,"analysis result/term number.xlsx",rowNames=T)

#-------------准备go数据格式---------------------
go_df <- go_result[,c(1:3,9,7)]
colnames(go_df) <- c("Category","ID","Term","Genes","adj_pvalue")

#-------------准备process（这里是 Term）---------
term_number <- toupper(term_number) #转换成大写

term_list <- go_df[,c(1,3)]
term_list$term_number <- paste0("T",c(1:nrow(term_list)))
term <- term_list$Term[term_list$term_number %in% term_number]  # 抽取出包含在 term_number 的term
term

#------------准备与term相对应的 target_gene ------------
term_df <- go_df[,c(3,4)]

target_gene <- term_df$Genes[term_df$Term %in% term]

target_gene

target_gene <- str_split(string = target_gene,pattern = ",") # 拆分字符串

target_gene <- unlist(target_gene) %>%  # 把列表变为字符串
              dplyr::tibble() # 变成表

colnames(target_gene) <- c("Genes")

table(duplicated(target_gene$Genes)) #查看重复

target_gene <- distinct(.data = target_gene,Genes,.keep_all = T) #去重

df_gene <- df[,c(1:2)]

target_gene <- dplyr::filter(df_gene,gene_symbol %in% target_gene$Genes)


#---------------GOplot分析和可视化-------------------

# Building the circ object
go_circ <- go_df
go_circ$Genes <- str_replace_all(go_circ$Genes,"/",",") #把/换成,
go_circ$Genes <- toupper(go_circ$Genes)

gene_circle <- df[,c(1:2)]
colnames(gene_circle) <- c("ID","logFC")
table(duplicated(gene_circle$ID)) #查看重复

circ <- circle_dat(go_circ,gene_circle)
#table(is.na(circ))

# Generating the binary matrix
gene_chord <- target_gene
colnames(gene_chord) <- c("ID","logFC")
gene_chord$ID <- toupper(gene_chord$ID)

process_chord <- term

chord<-chord_dat(circ,gene_chord,process_chord)

class(chord)

chord_data <- data.frame(chord)
colnames(chord_data) <- c(term,"logFC")

write.xlsx(chord_data,"analysis result/chord_data.xlsx",rowNames=T,colNaems=T)

# Creating the chord plot
# GOChord(chord)

# Excluding process with less than 5 assigned genes
# GOChord(chord, limit = c(0,5))

# Creating the chord plot genes ordered by logFC and a different logFC color scale
lfc <- circ$logFC
table(is.infinite(lfc)) # 无穷大的数据
lfc[is.infinite(lfc)]=NA # 把无穷大的数据转换成NA
lfc <- na.omit(lfc) # 去掉NA

#min <- ceiling(min(lfc))-1
#max <- ceiling(max(lfc))+1

my_col <- brewer.pal(length(term),"Set2")

GOChord <- GOChord(chord,
       # title="The chord diagram of GO enrichment",
        space=0.02,
        gene.order='logFC',
        lfc.col=c('firebrick3', 'white','royalblue3'), #上调下调颜色设置
    #    lfc.min=min,
    #    lfc.max=max,
        gene.size=2.5,  # 字体大小
        gene.space=0.25, # gene和环之间的距离
        #border.size=0.1
    #    ribbon.col=c("gray"), # The background color of the ribbons
        ribbon.col=my_col,#GO term 颜色设置
        process.label=10) # legend 字体大小


GOChord

my_theme <- theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), # 去除边框
        panel.background = element_blank(),#去除背景
        legend.title = element_blank(),
        axis.text.x =element_blank(), #删除x轴标签
        axis.text.y =element_blank(), #删除y轴标签
        axis.ticks = element_blank(), #删除x，y轴标签和刻度
        plot.title = element_text(size = 16, face = "bold", # labs(title) 字体大小
                                  hjust = 0.5), #标题居中
        plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=0.5),"cm"),) #边距

ggplot_chord <- as.ggplot(plot=GOChord) +my_theme+
                labs(title = "The chord diagram of GO enrichment")+
                xlab(label = "") + # 去掉 x 标签
                ylab(label = "") # 去掉 y 标签

ggsave("analysis result/GO_chord.png",ggplot_chord,width=1800, height =1500, dpi=180,units = "px")

print("The GO chord diagram can be found in the folder of analysis result")

ggplot_chord

 }

