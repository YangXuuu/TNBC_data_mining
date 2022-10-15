library(GOplot)
library(GOSim)
library(AnnotationDbi)
library(hgu133plus2.db)
library("GO.db")
library(plyr)
library(AnnotationDbi)
library(hgu133plus2.db)
anno<- AnnotationDbi::select(hgu133plus2.db,
                                  keys = probesets,
                                  columns = c("SYMBOL", "GENENAME"),
                                  keytype = "PROBEID")
anno2 <- subset(anno, !is.na(SYMBOL))
anno_grouped <- group_by(anno2, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

anno_multiple<- filter(anno_summarized, no_of_matches > 1)
ids_to_exlude <- (anno2$PROBEID %in% anno_multiple$PROBEID)
anno_final <- subset(anno2, !ids_to_exlude)

limma1<- read.csv('dataset1b_topgenes.csv')
limma2<- read.csv('dataset2b_topgenes.csv')
gfs1<- read.csv('GSE76275GFS_DEGs.csv')
colnames(gfs1)[1]<- 'PROBEID'
gfs1b<- merge(anno_final, gfs1, by.x = "PROBEID")
colnames(gfs1)[2]<-'ID'
gfs2<- read.csv('GSE43358_DEG.csv')
colnames(gfs2)[1]<- 'PROBEID'
gfs2b<- merge(anno_final, gfs2, by.x = "PROBEID")
colnames(gfs2)[2]<-'ID'


limma1_ora <- read.csv("ora_limma_dataset1.csv")
limma2_ora <- read.csv("ora_limma_dataset2.csv")
gfs1_ora<- read.csv("ora_gfs_dataset1.csv")
gfs2_ora<- read.csv("ora_gfs_dataset2.csv")


limma1b<- limma1[,-c(1,3)]
limma2b<- limma2[,-c(1,3)]
colnames(limma1b)[1]<-"ID"
colnames(limma2b)[1]<-"ID"
#colnames(limma1b)<- c("ID", "BP", "CC", "MF", "logFC","AveExpr","t" ,"P.Value","adj.P.Val", "B")
#limma1c<-  reshape2::melt(limma1b, id=c("ID","logFC","AveExpr","t" ,"P.Value","adj.P.Val", "B"))
limma1c <-ddply(limma1b, 'ID', dplyr::summarize, logFC =mean(logFC), AveExpr=head(AveExpr,1), t=head(t,1), P.Value=head(P.Value,1), adj.P.Val = head(adj.P.Val,1), B=head(B,1))
limma2c <-ddply(limma2b, 'ID', dplyr::summarize, logFC =mean(logFC), AveExpr=head(AveExpr,1), t=head(t,1), P.Value=head(P.Value,1), adj.P.Val = head(adj.P.Val,1), B=head(B,1))

limma1b_ora<- limma1_ora %>%
mutate(Category = "BP") %>%
dplyr::select(Category, term_id,term_name, intersection, p_value) 
colnames(limma1b_ora)<- c("Category","ID", "Term", "Genes", "adj_pval")
circlimma1<- circle_dat(limma1b_ora,limma1c )
chordlimma1 <- chord_dat(data = circlimma1, genes= limma1c[,1:2], process =limma1b_ora$Term)

limma2b_ora<- limma2_ora %>%
mutate(Category = "BP") %>%
dplyr::select(Category, term_id,term_name, intersection, p_value) 
colnames(limma2b_ora)<- c("Category","ID", "Term", "Genes", "adj_pval")
circlimma2<- circle_dat(limma2b_ora,limma2c )
chordlimma2 <- chord_dat(data = circlimma2, genes= limma2c[,1:2], process =limma2b_ora$Term)

gfs2b_ora<- gfs2_ora %>%
mutate(Category = "BP") %>%
dplyr::select(Category, term_id,term_name, intersection, p_value) 
colnames(gfs2b_ora)<- c("Category","ID", "Term", "Genes", "adj_pval")
circgfs2<- circle_dat(gfs2b_ora,gfs2c )
chordgfs2 <- chord_dat(data = circgfs2, genes= gfs2c[,1:2], process =gfs2b_ora$Term)


distance <- dist(chordlimma1 )
rownames(distance) <- row.names(chordlimma1 )
cluster <- hclust(distance)
M <- dim(chordlimma1)[2]
nterm <- M - 1
tmp <- NULL
for(r in 1:nrow(chordlimma1)){
tmp <- c(tmp, as.numeric(gsub(1, chordlimma1[r, (nterm + 1)], chordlimma1[r, 1:nterm])))
}
df <- data.frame(x = rep(cluster$order, each = nterm), y = rep(colnames(chordlimma1[,1:nterm]), length(rownames(chordlimma1))), z = tmp, lab = rep(rownames(chordlimma1), each = nterm))
df_o <- df[order(df$x),]
g <- ggplot(data = df_o, aes(x =lab, y = y)) + 
        geom_tile(aes( fill = z))+
		scale_fill_distiller(palette = "YlGnBu")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.y = element_text(size = 14), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
ggsave("heat1.png", width = 20, height=8)

distance <- dist(chordlimma2 )
rownames(distance) <- row.names(chordlimma2 )
cluster <- hclust(distance)
M <- dim(chordlimma2)[2]
nterm <- M - 1
tmp <- NULL
for(r in 1:nrow(chordlimma2)){
tmp <- c(tmp, as.numeric(gsub(1, chordlimma2[r, (nterm + 1)], chordlimma2[r, 1:nterm])))
}
df <- data.frame(x = rep(cluster$order, each = nterm), y = rep(colnames(chordlimma2[,1:nterm]), length(rownames(chordlimma2))), z = tmp, lab = rep(rownames(chordlimma2), each = nterm))
df_o <- df[order(df$x),]
g <- ggplot(data = df_o, aes(x =lab, y = y)) + 
        geom_tile(aes( fill = z))+
		scale_fill_distiller(palette = "YlGnBu")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.y = element_text(size = 14), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
ggsave("heat2.png", width = 10, height=3)




####IGNORE BELOW

library('org.Hs.eg.db')
library(GOxploreR)
library(GO.db)


limma1d<- cbind(limma1c,mapIds(org.Hs.eg.db,limma1c$ID, "ENTREZID", "SYMBOL"))
colnames(limma1d)[8]<- "Entrez"
bp<-Gene2GOTermAndLevel(genes = limma1d$Entrez, organism = "Homo sapiens", domain = "BP")
cc<-Gene2GOTermAndLevel(genes = limma1d$Entrez, organism = "Homo sapiens", domain = "CC")
mf<-Gene2GOTermAndLevel(genes = limma1d$Entrez, organism = "Homo sapiens", domain = "MF")
combine<-do.call("rbind", list(bp,cc,mf))
colnames(combine)<- c("Entrez", "GO", "Domain"   ,     "Level" )
GOIDs <- combine$GO
select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")


limma2d<- cbind(limma2c,mapIds(org.Hs.eg.db,limma2c$ID, "ENTREZID", "SYMBOL"))
colnames(limma2d)[8]<- "Entrez"


limma1b_ora<- limma1_ora %>%
mutate(Category = "BP") %>%
select(Category, term_id,term_name, intersection, p_value) 
colnames(limma1b_ora)<- c("Category","ID", "Term", "Genes", "adj_pval")
circlimma1<- circle_dat(limma1b_ora,limma1c )
chordlimma1 <- chord_dat(data = circlimma1, genes= limma1c[,1:2], process =limma1b_ora$Term)
M <- dim(chordlimma1)[2]
nterm <- M - 1
tmp <- NULL
for(r in 1:nrow(chordlimma1)){
tmp <- c(tmp, as.numeric(gsub(1, chordlimma1[r, (nterm + 1)], chordlimma1[r, 1:nterm])))
}
df <- data.frame(x = rep(cluster$order, each = nterm), y = rep(colnames(chordlimma1[,1:nterm]), length(rownames(chordlimma1))), z = tmp, lab = rep(rownames(chordlimma1), each = nterm))
df_o <- df[order(df$x),]
g <- ggplot(data = df_o, aes(x =lab, y = y)) + 
        geom_tile(aes( fill = z))+
		scale_fill_distiller(palette = "YlGnBu")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.y = element_text(size = 14), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
ggsave("heat1.png", width = 10, height=8)


limma2b<- limma2[,-c(1,2,3,4,5,6,7,8,9,12,13)]
colnames(limma1b)[1]<- "ID"
limma1c <-ddply(limma1b, 'ID', summarize, logFC =mean(logFC), AveExpr=head(AveExpr,1), t=head(t,1), P.Value=head(P.Value,1), adj.P.Val = head(adj.P.Val,1), B=head(B,1))
limma1b_ora<- limma1_ora %>%
mutate(Category = "BP") %>%
select(Category, term_id,term_name, intersection, p_value) 
colnames(limma1b_ora)<- c("Category","ID", "Term", "Genes", "adj_pval")
circlimma1<- circle_dat(limma1b_ora,limma1c )
chordlimma1 <- chord_dat(data = circlimma1, genes= limma1c[,1:2], process =limma1b_ora$Term)
M <- dim(chordlimma1)[2]
nterm <- M - 1
tmp <- NULL
for(r in 1:nrow(chordlimma1)){
tmp <- c(tmp, as.numeric(gsub(1, chordlimma1[r, (nterm + 1)], chordlimma1[r, 1:nterm])))
}
df <- data.frame(x = rep(cluster$order, each = nterm), y = rep(colnames(chordlimma1[,1:nterm]), length(rownames(chordlimma1))), z = tmp, lab = rep(rownames(chordlimma1), each = nterm))
df_o <- df[order(df$x),]
g <- ggplot(data = df_o, aes(x =lab, y = y)) + 
        geom_tile(aes( fill = z))+
		scale_fill_distiller(palette = "YlGnBu")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.y = element_text(size = 14), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
ggsave("heat1.png", width = 10, height=8)