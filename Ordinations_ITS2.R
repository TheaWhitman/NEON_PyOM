library(phyloseq)
library(wesanderson)
library(ggplot2)
library(vegan)
library(dplyr)

ps.hell=readRDS("../data/CornellITS2/ps.ITS2.hell")
ps.hell

ps.norm = transform_sample_counts(ps.hell,function(x) x^2)
ps.norm.genus = tax_glom(ps.norm,"Genus")
ps.hell.genus = transform_sample_counts(ps.norm.genus,function(x) x^0.5)

ps.norm.class = tax_glom(ps.norm,"Class")
ps.hell.class = transform_sample_counts(ps.norm.class,function(x) x^0.5)

ps.hell = ps.hell.class

Sums = taxa_sums(ps.hell)
Genus = data.frame(as(tax_table(ps.hell),"matrix"))$Genus
OTUs = names(taxa_sums(ps.hell))

df = data.frame(Sums,OTUs, Genus)%>%
  arrange(-Sums)
head(df)
df

head(tax_table(ps.hell))

# Full dataset analysis

ord = ordinate(ps.hell, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(ps.hell), "data.frame")
d = phyloseq::distance(ps.hell, method = "bray")

d.adonis = adonis(d~sample_data(ps.hell)$Soil_Trtmt+sample_data(ps.hell)$Amdmt+sample_data(ps.hell)$Day+
                    sample_data(ps.hell)$Soil_Trtmt*sample_data(ps.hell)$Amdmt+
                    sample_data(ps.hell)$Soil_Trtmt*sample_data(ps.hell)$Day,
                  df)
d.adonis

p = plot_ordination(ps.hell, ord, axes=c(1,2), type="samples", color="Day",shape="Soil_Trtmt")
#palette = c(wes_palette("Darjeeling1"))[c(2,1,4,3,5)]
palette = c(wes_palette("Cavalcanti1"))[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + scale_shape_manual(values=c(15,16,17,18,3))
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title="Day"), color=guide_legend(title="Amendment"))
p


### Additional adonis

ps.unamended = subset_samples(ps.hell,Amdmt=="Soil")
ps.unamended = subset_samples(ps.unamended,Day=="1")


d = phyloseq::distance(ps.unamended, method = "bray")
df = as(sample_data(ps.unamended), "data.frame")

d.adonis = adonis(d~sample_data(ps.unamended)$ph+sample_data(ps.unamended)$CEC+sample_data(ps.unamended)$Ca+sample_data(ps.unamended)$Mg+sample_data(ps.unamended)$Na+sample_data(ps.unamended)$K+sample_data(ps.unamended)$totalCarbon+sample_data(ps.unamended)$totalNitrogen,
                  df)
d.adonis

### Hawaii
Soil = "Hawaii"
physeq = prune_samples(sample_data(ps.hell)$Soil_Trtmt==Soil,ps.hell)
ord = ordinate(physeq, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray")

d.adonis = adonis(d~sample_data(physeq)$Amdmt+
                    sample_data(physeq)$Day,
                  df)
d.adonis

p = plot_ordination(physeq, ord, axes=c(1,2), type="samples", color="Amdmt",shape="Day")
palette = c(wes_palette("Cavalcanti1"))
palette = palette[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title=""), color=guide_legend(title=""))
p = p + facet_wrap(~Soil_Name)
p

### Alaska 
Soil = "Alaska"
physeq = prune_samples(sample_data(ps.hell)$Soil_Trtmt==Soil,ps.hell)
ord = ordinate(physeq, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray")

d.adonis = adonis(d~sample_data(physeq)$Amdmt+
                    sample_data(physeq)$Day,
                  df)
d.adonis

options(repr.plot.width = 3, repr.plot.height = 2)

p = plot_ordination(physeq, ord, axes=c(1,2), type="samples", color="Amdmt",shape="Day")
palette = c(wes_palette("Cavalcanti1"))
palette = palette[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title=""), color=guide_legend(title=""))
p = p + facet_wrap(~Soil_Name)
p

## Utah

Soil = "Utah"
physeq = prune_samples(sample_data(ps.hell)$Soil_Trtmt==Soil,ps.hell)
ord = ordinate(physeq, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray")

d.adonis = adonis(d~sample_data(physeq)$Amdmt+
                    sample_data(physeq)$Day,
                  df)
d.adonis

options(repr.plot.width = 3, repr.plot.height = 2)

p = plot_ordination(physeq, ord, axes=c(1,2), type="samples", color="Amdmt",shape="Day")
palette = c(wes_palette("Cavalcanti1"))
palette = palette[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title=""), color=guide_legend(title=""))
p = p + facet_wrap(~Soil_Name)
p


## New York
Soil = "New York"
physeq = prune_samples(sample_data(ps.hell)$Soil_Trtmt==Soil,ps.hell)
ord = ordinate(physeq, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray")

d.adonis = adonis(d~sample_data(physeq)$Amdmt+
                    sample_data(physeq)$Day,
                  df)
d.adonis

options(repr.plot.width = 3, repr.plot.height = 2)

p = plot_ordination(physeq, ord, axes=c(1,2), type="samples", color="Amdmt",shape="Day")
palette = c(wes_palette("Cavalcanti1"))
palette = palette[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title=""), color=guide_legend(title=""))
p = p + facet_wrap(~Soil_Name)
p

## Florida 

Soil = "Florida"
physeq = prune_samples(sample_data(ps.hell)$Soil_Trtmt==Soil,ps.hell)
ord = ordinate(physeq, "NMDS", "bray", k=2, trymax=500)
ord
df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray")

d.adonis = adonis(d~sample_data(physeq)$Amdmt+
                    sample_data(physeq)$Day,
                  df)
d.adonis

options(repr.plot.width = 3, repr.plot.height = 2)

p = plot_ordination(physeq, ord, axes=c(1,2), type="samples", color="Amdmt",shape="Day")
palette = c(wes_palette("Cavalcanti1"))
palette = palette[c(3,1,2)]
p = p + theme_bw()
p = p + scale_color_manual(values=palette)
p = p + geom_point(size=3)
p = p + guides(shape = guide_legend(title=""), color=guide_legend(title=""))
p = p + facet_wrap(~Soil_Name)
p

