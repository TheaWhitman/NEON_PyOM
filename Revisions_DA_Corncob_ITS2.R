library(phyloseq)
library(wesanderson)
library(ggplot2)
library(vegan)
library(dplyr)
library(corncob)
library(reshape)
library(RColorBrewer)

ps=readRDS("../data/CornellITS2/ps.ITS2")
ps.norm=readRDS("../data/CornellITS2/ps.ITS2.norm")
ps

# If trying with genus-glommed taxa
ps.norm.genus = tax_glom(ps.norm,"Genus")
ps.genus = tax_glom(ps,"Genus")
ps.norm.fam = tax_glom(ps.norm,"Family")
ps.fam = tax_glom(ps,"Family")


ps = ps.fam
ps.norm = ps.norm.fam



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

results[[5]][[1]]$significant_models

sigOTUs = c(
  row.names(tax_table(r[[1]]$data)[!is.na(r[[1]]$p_fdr) & r[[1]]$p_fdr<cutoff]),
  #row.names(tax_table(r[[2]]$data)[!is.na(r[[2]]$p_fdr) & r[[2]]$p_fdr<cutoff]),
  #row.names(tax_table(r[[3]]$data)[!is.na(r[[3]]$p_fdr) & r[[3]]$p_fdr<cutoff]),
  #row.names(tax_table(r[[4]]$data)[!is.na(r[[4]]$p_fdr) & r[[4]]$p_fdr<cutoff]),
  row.names(tax_table(r[[5]]$data)[!is.na(r[[5]]$p_fdr) & r[[5]]$p_fdr<cutoff]))
levels(as.factor(sigOTUs))

# Which are significant, in how many soils?
data.frame(sigOTUs) %>%
  group_by(sigOTUs)%>%
  summarize(Count=n())%>%
  arrange(-Count)

# Example plot
p = plot(r[[5]])
p = p + theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
p
# Onlye one taxon increased over time

# Check out a single model to confirm structure and coefficients
results[[5]][[1]]$significant_models[[1]]

results[[5]][[1]]$significant_models[[1]]$coefficients[1:2,1]


mu = data.frame(t(as.matrix(results[[1]][[1]]$significant_models[[1]]$coefficients[1:2,1])))
mu

df = data.frame()
for (i in c(1,5)){
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





##### Running with glommed genus

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
results[[1]][[1]]$significant_models[[1]]

results[[1]][[1]]$significant_models[[1]]$coefficients[1:3,1]

mu = data.frame(t(as.matrix(results[[1]][[1]]$significant_models[[1]]$coefficients[1:3,1])))
mu


results[[10]][[1]]$significant_taxa

# For each list, and for each significant taxon in that list, extract the mu coefficients and the p_fdr
# Currently not running for 6, because no significant taxa

df = data.frame()
#for (i in c(1:5,7,9:10)){ # For classs
for (i in c(1,3:5,8:10)){  
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


#write.csv(joined,"../data/CornellITS2/TableS10.ITS2.Responders.Revised.csv")

FCcutoff = 2

joined.plot = joined %>%
  dplyr::filter(log2FC>0)

joined.plot2 = joined %>%
  dplyr::filter(log2FC<0)



p = ggplot(joined.plot,aes(y=log2FC,color=Phylum,x=Family))
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

summary1=joined.plot2%>% 
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

###### 
