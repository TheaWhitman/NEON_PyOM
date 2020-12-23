library(phyloseq)
library(wesanderson)
library(ggplot2)
library(vegan)
library(dplyr)
library(corncob)
library(reshape)
library(RColorBrewer)

ps=readRDS("../data/Cornell16S/ps.16S")
ps.norm=readRDS("../data/Cornell16S/ps.16S.norm")
ps

#ps.class = tax_glom(ps,"Class")
#ps.norm.class = tax_glom(ps.norm,"Class")
#ps = ps.class
#ps.norm = ps.norm.class


# Most interested in the final timepoint, for OM and for PyOM

Factors = expand.grid(Soil_Trtmt=c("Hawaii","Alaska","Utah","New York","Florida"))
Factors


# Make function to test for differential abundance
# controlling for differential variance
da_analysis = function(Factors){
  Soil_Trtmt=paste(Factors["Soil_Trtmt"])
  ps.DA = prune_samples(sample_data(ps)$Soil_Trtmt==Soil_Trtmt & sample_data(ps)$Day %in% c("1","26") & (sample_data(ps)$Amdmt == "Soil"),ps)
  ps.DA = prune_taxa(taxa_sums(ps.DA)>0,ps.DA)
  print(ps.DA)
  AbundTaxa = taxa_names(filter_taxa(ps.norm, function(x) mean(x) > 0.0001, TRUE))
  ps.DA = prune_taxa(AbundTaxa,ps.DA)
  print(ps.DA)
  dT = differentialTest(formula = ~ Day,
                        phi.formula = ~ Day,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Day,
                        test = "Wald", boot = FALSE,
                        data = ps.DA,
                        fdr_cutoff = 0.05)
  results=list(dT,Soil_Trtmt)
  return(results)
}

results = apply(Factors,1,da_analysis)
# This produces a list that holds the results of the differential abundance test

cutoff=0.05
r = lapply(results, function(l) l[[1]])

sigOTUs = c(
  row.names(tax_table(r[[1]]$data)[!is.na(r[[1]]$p_fdr) & r[[1]]$p_fdr<cutoff]),
  row.names(tax_table(r[[2]]$data)[!is.na(r[[2]]$p_fdr) & r[[2]]$p_fdr<cutoff]),
  row.names(tax_table(r[[3]]$data)[!is.na(r[[3]]$p_fdr) & r[[3]]$p_fdr<cutoff]),
  row.names(tax_table(r[[4]]$data)[!is.na(r[[4]]$p_fdr) & r[[4]]$p_fdr<cutoff]),
  row.names(tax_table(r[[5]]$data)[!is.na(r[[5]]$p_fdr) & r[[5]]$p_fdr<cutoff]))
levels(as.factor(sigOTUs))

# Which are significant, in how many soils?
data.frame(sigOTUs) %>%
  group_by(sigOTUs)%>%
  summarize(Count=n())%>%
  arrange(-Count)

# Example plot
p = plot(r[[1]])
p = p + theme(axis.title.y=element_blank(),
              #axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
p
# Lots of taxa increased over time

# Check out a single model to confirm structure and coefficients
results[[1]][[1]]$significant_models[[1]]

results[[1]][[1]]$significant_models[[1]]$coefficients[1:2,1]


mu = data.frame(t(as.matrix(results[[1]][[1]]$significant_models[[1]]$coefficients[1:2,1])))
mu

df = data.frame()
for (i in c(1:5)){
  r = results[[i]][[1]]
  Soil = try(results[[i]][[2]])
  for (j in 1:length(r$significant_taxa)){
    sig_models = try(r$significant_models[[j]],silent=TRUE)
    mu = try(data.frame(t(as.matrix(sig_models$coefficients[1:2,1]))),silent=TRUE)
    se = try(data.frame(t(as.matrix(sig_models$coefficients[1:2,2]))),silent=TRUE)
    mu = cbind(mu,se)
    p_fdr = try(r$p_fdr[r$significant_taxa][j],silent=TRUE)
    mu$p_fdr = try(p_fdr,silent=TRUE)
    mu$Soil = try(Soil,silent=TRUE)
    row.names(mu)= try(paste(row.names(try(data.frame(p_fdr))),"_",i,sep=""),silent=TRUE)
    colnames(mu) = c("mu.Intercept","mu.Day26","se.Intercept","se.Day26","p_fdr","Soil")
    df = try(rbind(df,mu))
  }
}

head(df)


df$OTU = row.names(df)
df$OTU = sub("_[0-9]","",df$OTU)
head(df$OTU)


# Calculate the fold-change value using corncob's inverse logit function and our estimates
df$FC = corncob::invlogit(df$mu.Intercept+df$mu.Day26)/corncob::invlogit(df$mu.Intercept)


# Get the log2-fold change
df$log2FC = log(df$FC,base=2)
# Get the baseline relative abundance (of that OTU in bulk soil)
df$Relabund = (corncob::invlogit(df$mu.Intercept))


SigOTUs = levels(as.factor(df$OTU))
pruned = prune_taxa(SigOTUs,ps.norm)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
joined = merge(taxtab,df,by=c("OTU"))

ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium")

#Burkholderia-Caballeronia-Paraburkholderia
#Colstridium senu stricto 10
#Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
#Plactomycetales bacterium Ellin 6207

joined$Phylum = ifelse(joined$Phylum=="WPS-2","Eremiobacterota (WPS-2)",paste(joined$Phylum))

joined = joined %>%
  mutate(Name = ifelse(Genus %in% ignoreList |is.na(Genus),ifelse(Family %in% ignoreList |is.na(Family),ifelse(Class %in% ignoreList |is.na(Class),ifelse(Phylum %in% ignoreList |is.na(Phylum),paste(OTU),
                                                                                                                                                          paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))%>%
  mutate(Name = ifelse(Name == "Burkholderia-Caballeronia-Paraburkholderia","*Burkholderia",Name))%>%
  mutate(Name = ifelse(Name == "Clostridium sensu stricto 10","Clostridium",Name))%>%
  mutate(Name = ifelse(Name == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","*Rhizobium",Name))%>%
  mutate(Name = ifelse(Name == "Planctomycetales bacterium Ellin6207","Planctomycetales",Name))%>%
  mutate(Name = ifelse(Name == "uncultured thaumarchaeote","Thaumarchaeota",Name))
head(joined$Name)

TaxonOrder = joined %>%
  dplyr::select(Phylum,Genus,OTU,Name)%>%
  dplyr::arrange(Phylum)
OTUOrder = unique(TaxonOrder$OTU)
GenusOrder = unique(TaxonOrder$Genus)
NameOrder = unique(TaxonOrder$Name)

joined$OTU = factor(joined$OTU, levels = OTUOrder)
joined$Genus = factor(joined$Genus, levels = GenusOrder)
joined$Name = factor(joined$Name, levels = NameOrder)
joined$Soil_Trtmt = factor(joined$Soil, levels = c("Hawaii","Alaska","Utah","New York","Florida"))

levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Hawaii"] = "Hydrudand"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="New York"] = "Fragiudept"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Alaska"] = "Cryaquept"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Utah"] = "Haplocalcid"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Florida"] = "Quartzipsamment"

joined.plot = joined %>%
  dplyr::filter(log2FC>0)

joined.plot2 = joined %>%
  dplyr::filter(log2FC<0)

p = ggplot(joined.plot,aes(y=log2FC,color=Phylum,x=Name))
p = p + theme_bw()
p = p + geom_point() + facet_grid(~Soil_Trtmt,scales="free_x",space="free_x")
p = p + theme(axis.text.x=element_text(angle=90,size=7,face="italic",vjust=0,hjust=1))
p = p + ylab("log2-fold change in relative abundance with amendment") + xlab("")
p


# Summarizing response for positive responders

summary1=joined.plot%>% 
  distinct(OTU)%>%
  summarize(Count=n())%>%
  arrange(-Count)

summary2=joined.plot%>%
  filter(!(Genus %in% ignoreList))%>%
  group_by(Genus,Soil_Trtmt)%>%
  summarize(Count=n())%>%
  arrange(Genus,-Count)%>%
  group_by(Genus)%>%
  summarize(Sum=sum(Count),Soils=n())%>%
  arrange(-Sum)

summary1
summary2


############################# Running Amdmt analysis at class level

# Most interested in the final timepoint, for OM and for PyOM

ntaxa(ps.norm)

Factors = expand.grid(Soil_Trtmt=c("Hawaii","Alaska","Utah","New York","Florida"),Amdmt=c("OM","PyOM"))
Factors

# Make function to test for differential abundance
# controlling for differential variance
da_analysis = function(Factors){
  Soil_Trtmt=paste(Factors["Soil_Trtmt"])
  Amdmt=paste(Factors["Amdmt"])
  ps.DA = prune_samples(sample_data(ps)$Soil_Trtmt==Soil_Trtmt & sample_data(ps)$Day %in% c("10","26") & (sample_data(ps)$Amdmt == Amdmt | sample_data(ps)$Amdmt == "Soil"),ps)
  ps.DA = prune_taxa(taxa_sums(ps.DA)>0,ps.DA)
  print(ps.DA)
  AbundTaxa = taxa_names(filter_taxa(ps.norm, function(x) mean(x) > 0.0001, TRUE))
  ps.DA = prune_taxa(AbundTaxa,ps.DA)
  print(ps.DA)
  dT = differentialTest(formula = ~ Day+Amdmt,
                        phi.formula = ~ Day+Amdmt,
                        formula_null = ~ Day,
                        phi.formula_null = ~ Day+Amdmt,
                        test = "Wald", boot = FALSE,
                        data = ps.DA,
                        fdr_cutoff = 0.05)
  results=list(dT,Amdmt,Soil_Trtmt)
  return(results)
}

results = apply(Factors[,],1,da_analysis)
# This produces a list that holds the results of the differential abundance test

cutoff=0.05
r = lapply(results, function(l) l[[1]])
sigOTUs = c(
  row.names(tax_table(r[[1]]$data)[!is.na(r[[1]]$p_fdr) & r[[1]]$p_fdr<cutoff]),
  row.names(tax_table(r[[2]]$data)[!is.na(r[[2]]$p_fdr) & r[[2]]$p_fdr<cutoff]),
  row.names(tax_table(r[[3]]$data)[!is.na(r[[3]]$p_fdr) & r[[3]]$p_fdr<cutoff]),
  row.names(tax_table(r[[4]]$data)[!is.na(r[[4]]$p_fdr) & r[[4]]$p_fdr<cutoff]),
  row.names(tax_table(r[[5]]$data)[!is.na(r[[5]]$p_fdr) & r[[5]]$p_fdr<cutoff]),
  row.names(tax_table(r[[6]]$data)[!is.na(r[[6]]$p_fdr) & r[[6]]$p_fdr<cutoff]),
  row.names(tax_table(r[[7]]$data)[!is.na(r[[7]]$p_fdr) & r[[7]]$p_fdr<cutoff]),
  row.names(tax_table(r[[8]]$data)[!is.na(r[[8]]$p_fdr) & r[[8]]$p_fdr<cutoff]),
  row.names(tax_table(r[[9]]$data)[!is.na(r[[9]]$p_fdr) & r[[9]]$p_fdr<cutoff]),
  row.names(tax_table(r[[10]]$data)[!is.na(r[[10]]$p_fdr) & r[[10]]$p_fdr<cutoff]))
levels(as.factor(sigOTUs))

data.frame(sigOTUs) %>%
  group_by(sigOTUs)%>%
  summarize(Count=n())%>%
  arrange(-Count)


# Check out a single model to confirm structure and coefficients
results[[10]][[1]]$significant_models[[1]]

results[[1]][[1]]$significant_models[[1]]$coefficients[1:3,1]

mu = data.frame(t(as.matrix(results[[1]][[1]]$significant_models[[1]]$coefficients[1:3,1])))
mu


results[[8]][[1]]$significant_taxa

# For each list, and for each significant taxon in that list, extract the mu coefficients and the p_fdr
# Currently not running for 6, because no significant taxa

df = data.frame()
#for (i in c(1:5,7,9:10)){ # For classs
for (i in c(1:5,7:10)){  
 r = results[[i]][[1]]
  Amdmt = try(results[[i]][[2]])
  Soil = try(results[[i]][[3]])
  for (j in 1:length(r$significant_taxa)){
    sig_models = try(r$significant_models[[j]],silent=TRUE)
    mu = try(data.frame(t(as.matrix(sig_models$coefficients[1:3,1]))),silent=TRUE)
    se = try(data.frame(t(as.matrix(sig_models$coefficients[1:3,2]))),silent=TRUE)
    mu = cbind(mu,se)
    p_fdr = try(r$p_fdr[r$significant_taxa][j],silent=TRUE)
    mu$p_fdr = try(p_fdr,silent=TRUE)
    mu$Amdmt = try(Amdmt,silent=TRUE)
    mu$Soil = try(Soil,silent=TRUE)
    row.names(mu)= try(paste(row.names(try(data.frame(p_fdr))),"_",i,sep=""),silent=TRUE)
    colnames(mu) = c("mu.Intercept","mu.Day26","mu.Amdmt","se.Intercept","se.Day26","se.Amdmt","p_fdr","Amdmt","Soil")
    df = try(rbind(df,mu))
  }
}

df$OTU = row.names(df)
df$OTU = sub("_[0-9]","",df$OTU)
head(df$OTU)

# Calculate the fold-change value using corncob's inverse logit function and our estimates
df$FC.10 = corncob::invlogit(df$mu.Intercept+df$mu.Amdmt)/corncob::invlogit(df$mu.Intercept)
df$FC.26 = corncob::invlogit(df$mu.Intercept+df$mu.Day26+df$mu.Amdmt)/corncob::invlogit(df$mu.Intercept+df$mu.Day26)
df$FC = (df$FC.10+df$FC.26)/2

# Get the log2-fold change
df$log2FC = log(df$FC,base=2)
# Get the baseline relative abundance (of that OTU in bulk soil)
df$Relabund = (corncob::invlogit(df$mu.Intercept)+corncob::invlogit(df$mu.Intercept+df$mu.Day26))/2
df=df[,c(7:10,13:15)]

SigOTUs = levels(as.factor(df$OTU))
pruned = prune_taxa(SigOTUs,ps.norm)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
joined = merge(taxtab,df,by=c("OTU"))


ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium")

#Burkholderia-Caballeronia-Paraburkholderia
#Colstridium senu stricto 10
#Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
#Plactomycetales bacterium Ellin 6207

joined$Phylum = ifelse(joined$Phylum=="WPS-2","Eremiobacterota (WPS-2)",paste(joined$Phylum))

joined = joined %>%
  mutate(Name = ifelse(Genus %in% ignoreList |is.na(Genus),ifelse(Family %in% ignoreList |is.na(Family),ifelse(Class %in% ignoreList |is.na(Class),ifelse(Phylum %in% ignoreList |is.na(Phylum),paste(OTU),
                                                                                                                                                          paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))%>%
  mutate(Name = ifelse(Name == "Burkholderia-Caballeronia-Paraburkholderia","*Burkholderia",Name))%>%
  mutate(Name = ifelse(Name == "Clostridium sensu stricto 10","Clostridium",Name))%>%
  mutate(Name = ifelse(Name == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","*Rhizobium",Name))%>%
  mutate(Name = ifelse(Name == "Planctomycetales bacterium Ellin6207","Planctomycetales",Name))%>%
  mutate(Name = ifelse(Name == "uncultured thaumarchaeote","Thaumarchaeota",Name))
head(joined$Name)


TaxonOrder = joined %>%
  dplyr::select(Phylum,Genus,OTU,Name)%>%
  dplyr::arrange(Phylum)
OTUOrder = unique(TaxonOrder$OTU)
GenusOrder = unique(TaxonOrder$Genus)
NameOrder = unique(TaxonOrder$Name)


joined$OTU = factor(joined$OTU, levels = OTUOrder)
joined$Genus = factor(joined$Genus, levels = GenusOrder)
joined$Name = factor(joined$Name, levels = NameOrder)
joined$Soil_Trtmt = factor(joined$Soil, levels = c("Hawaii","Alaska","Utah","New York","Florida"))

levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Hawaii"] = "Hydrudand"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="New York"] = "Fragiudept"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Alaska"] = "Cryaquept"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Utah"] = "Haplocalcid"
levels(joined$Soil_Trtmt)[levels(joined$Soil_Trtmt)=="Florida"] = "Quartzipsamment"

FCcutoff = 2

joined.plot = joined %>%
  dplyr::filter(log2FC>2)

joined.plot2 = joined %>%
  dplyr::filter(log2FC<0)



p = ggplot(joined.plot,aes(y=log2FC,color=Phylum,x=Class))
p = p + theme_bw()
p = p + geom_jitter(width=0.1) + facet_grid(~Amdmt~Soil_Trtmt,scales="free_x",space="free_x")
#p = p + geom_point() + facet_wrap(~Soil_Trtmt~Amdmt,scales="free_x")
#p = p + geom_errorbar(aes(ymax=log2FC+1.96*se,ymin=log2FC-1.96*se),width=0.2)
p = p + theme(axis.text.x=element_text(angle=90,size=10,face="italic",vjust=0,hjust=1))
#p = p + theme(strip.text.x=element_text(angle=90,size=10))
p = p + ylab("log2-fold change in relative abundance with amendment") + xlab("")
p = p + ylim(c(0,12.7))
#p = p + scale_color_manual(values=palette)
p

# Summarizing response for positive responders



# Want to report total unique OTUs that were positive responders

summary1=joined.plot%>% 
  group_by(Amdmt)%>%
  distinct(OTU)%>%
  summarize(Count=n())%>%
  arrange(-Count)

summary3=joined.plot%>% 
  group_by(Soil_Trtmt,Amdmt)%>%
  distinct(OTU)%>%
  group_by(Soil_Trtmt,OTU)%>%
  summarize(Count=n())%>%
  filter(Count==2)%>%
  summarize(n())

summary2=joined.plot%>%
  #filter(!(Genus %in% ignoreList))%>%
  group_by(Class,Amdmt,Soil_Trtmt)%>%
  summarize(Count=n())%>%
  arrange(Class,-Count)%>%
  group_by(Class,Amdmt)%>%
  summarize(Sum=sum(Count),Soils=n())%>%
  arrange(Class,Amdmt)

summary1
summary3
summary2


### Correls

# Want to make analogous figures to ISME paper
# So, need to match up OTUs across days
# Goal: OTU taxonomy, Soil, Day, Response to PyOM, Response to OM

d.cast.l2FC = cast(joined, Soil_Trtmt+OTU+Phylum~Amdmt,value=c("log2FC"))
d.cast.padj = cast(joined, Soil_Trtmt+OTU+Phylum~Amdmt,value=c("p_fdr"))
d.cast = merge(d.cast.l2FC,d.cast.padj,by=c("Soil_Trtmt","OTU","Phylum"))

#colnames(d.cast) = c(colnames(d.cast)[1:3],"Relabund","logFC.OM","logFC.PyOM","Relabund.2","logFC.padj.OM","logFC.padj.PyOM")
colnames(d.cast) = c(colnames(d.cast)[1:3],"logFC.OM","logFC.PyOM","logFC.padj.OM","logFC.padj.PyOM")
#d.cast=d.cast[,c(1:6,8:9)]
head(d.cast)



d.plot = d.cast%>%
  dplyr::filter(!is.na(logFC.padj.PyOM) & !is.na(logFC.padj.OM))

d.plot$Soil_Name=d.plot$Soil_Trtmt

d.plot$Soil_Name = ordered(d.plot$Soil_Name, levels = c("Hydrudand", "Cryaquept", "Haplocalcid","Fragiudept","Quartzipsamment"))

p = ggplot(d.plot,aes(x=logFC.PyOM,y=logFC.OM,shape=Soil_Name,color=Phylum,fill=Soil_Name))
p = p + geom_point(size=3,alpha=0.7)
p = p + geom_hline(yintercept = 0.0, linetype=2) + theme_bw()
p = p + geom_vline(xintercept = 0.0, linetype=2) + theme_bw()
p = p + xlab("Response to PyOM (log2-fold change)")
p = p + ylab("Response to OM (log2-fold change)")

#palette = c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2")[2:4])
#palette = palette[c(2,1,4,3,5,6)]
#p = p + scale_fill_manual(values=palette)+ scale_color_manual(values=palette) + scale_shape_manual(values=c(21:25))
p = p + guides(shape = guide_legend(title="Soil Type"), color=guide_legend(title="Soil Type"), fill=guide_legend(title="Soil Type"))
p = p + facet_wrap(~Soil_Name)
p

Relabunds = subset_taxa(ps.norm, taxa_names(ps.norm) %in% d.cast$OTU)
Relabunds = psmelt(Relabunds)
Relabunds = Relabunds %>%
  #filter(Amdmt=="Soil")%>%
  filter(Day!=1)%>%
  group_by(Soil_Trtmt,OTU)%>%
  summarize(MeanAbund = mean(Abundance))%>%
  filter(MeanAbund>0)

levels(Relabunds$Soil_Trtmt) = c("Hydrudand", "Cryaquept", "Haplocalcid","Fragiudept","Quartzipsamment")


d.cast.2 = merge(d.cast,Relabunds,by=c("Soil_Trtmt","OTU"))
dim(d.cast.2)
dim(d.cast)


d.plot = d.cast.2%>%
  filter(!is.na(logFC.padj.PyOM) & !is.na(logFC.padj.OM))%>%
  filter(!(Phylum %in% c("BRC1", "Cyanobacteria","Thaumarchaeota")))

d.plot$Soil_Name=d.plot$Soil_Trtmt

d.plot$Soil_Name = ordered(d.plot$Soil_Name, levels = c("Hydrudand", "Cryaquept", "Haplocalcid","Fragiudept","Quartzipsamment"))

p = ggplot()
p = p + geom_abline(slope=m,intercept=b,colour="grey")
?geom_abline
p = p + theme_bw()
p = p + geom_point(data=d.plot,aes(x=logFC.PyOM,y=logFC.OM,shape=Soil_Name,color=Soil_Name,fill=Soil_Name,size=MeanAbund), alpha=0.5)
p = p + geom_hline(yintercept = 0.0, linetype=2) + theme_bw()
p = p + geom_vline(xintercept = 0.0, linetype=2) + theme_bw()
p = p + xlab("Response to PyOM (log2-fold change)")
p = p + ylab("Response to OM (log2-fold change)")
p = p + theme(strip.text = element_text(face="italic"))
palette = c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2")[2:4])
palette = palette[c(2,1,4,3,5,6)]
p = p + scale_fill_manual(values=palette)+ scale_color_manual(values=palette) + scale_shape_manual(values=c(21:25))
p = p + guides(size = guide_legend(title="Mean Relative\nAbundance"), shape = guide_legend(title="Soil Type"), color=guide_legend(title="Soil Type"), fill=guide_legend(title="Soil Type"))
p = p + facet_wrap(~Phylum)
p


########### Fitting fit line

d.plot.lm = lm(data=d.plot,logFC.OM~logFC.PyOM)
summary(d.plot.lm)

b = d.plot.lm$coefficients[1]
m = d.plot.lm$coefficients[2]
b
m

