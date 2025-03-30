# Population Selection Analysis
This is about calculating and visualizing genetic diversity metrics such as π (nucleotide diversity), Tajima’s D (neutrality test), and Fst (population differentiation).  
Required Tools and Files:    
vcftools: For calculating π, Tajima’s D, and Fst.  
R: For visualizing the results using ggplot2 and ggridges.  
Data Files:  
VCF File: data/all.vcf  
Population Lists: data/pop_wild.list and data/pop_cultivated.list (lists of individuals in wild and cultivated populations).  

# Step 1: Calculate π Values

In this step, we calculate nucleotide diversity (π) for both wild and cultivated populations using a sliding window approach.
```
# Parameter settings
VCF="data/all.vcf"
WINDOW=100000 # 100 kb window
STEP=10000 # 10 kb sliding step

# Wild population
vcftools --vcf $VCF \
--window-pi $WINDOW \
--window-pi-step $STEP \
--keep data/pop_wild.list \
--out results/pi/wild_population

# Cultivated population
vcftools --vcf $VCF \
--window-pi $WINDOW \
--window-pi-step $STEP \
--keep data/pop_cultivated.list \
--out results/pi/cultivated_population
```
# Step 2: Visualize π Values

After calculating π, we can visualize the values across the genome with a Manhattan plot.  
In R
```
library(ggplot2)
library(dplyr)

# Read data
wild_pi <- read.table("results/pi/wild_population.windowed.pi", header=T)
cultivated_pi <- read.table("results/pi/cultivated_population.windowed.pi", header=T)

# Merge data
combined <- bind_rows(
mutate(wild_pi, Population = "Wild"),
mutate(cultivated_pi, Population = "Cultivated")
)

# Generate Manhattan plot
png("figures/pi_manhattan.png", width=12, height=6, units="in", res=300)
ggplot(combined, aes(x=BIN_START/1e6, y=PI, color=CHROM)) +
geom_point(alpha=0.6) +
facet_grid(Population ~ .) +
labs(x="Position (Mb)", y="π Value") +
theme_bw() +
scale_color_viridis_d() # Colorblind-friendly palette
dev.off()
```
# Step 3: Calculate Tajima’s D
We calculate Tajima’s D to test for deviations from neutrality in both populations.
```
# Wild population
vcftools --vcf $VCF \
--TajimaD $WINDOW \
--keep data/pop_wild.list \
--out results/tajimaD/wild_population

# Cultivated population
vcftools --vcf $VCF \
--TajimaD $WINDOW \
--keep data/pop_cultivated.list \
--out results/tajimaD/cultivated_population
```
# Step 4: Visualize Tajima’s D

After calculating Tajima’s D, we can visualize its distribution in both populations using a density plot.  
In R
```
library(ggridges)

# Read data
wild_tajima <- read.table("results/tajimaD/wild_population.TajimaD", header=T)
cultivated_tajima <- read.table("results/tajimaD/cultivated_population.TajimaD", header=T)

# Merge data
combined <- bind_rows(
mutate(wild_tajima, Population = "Wild"),
mutate(cultivated_tajima, Population = "Cultivated")
)

# Generate density plot
png("figures/tajimad_density.png", width=8, height=6, units="in", res=300)
ggplot(combined, aes(x=TajimaD, y=Population, fill=Population)) +
geom_density_ridges(alpha=0.6, scale=0.9) +
labs(x="Tajima's D", y="") +
theme_minimal() +
scale_fill_manual(values=c("#1b9e77", "#d95f02"))
dev.off()
```
# Step 5: Calculate Fst
We calculate Fst to measure the genetic differentiation between the wild and cultivated populations using a sliding window approach.
```
vcftools --vcf $VCF \
--fst-window-size $WINDOW \
--fst-window-step $STEP \
--weir-fst-pop data/pop_wild.list \
--weir-fst-pop data/pop_cultivated.list \
--out results/fst/wild_cultivated
```
# Step 6: Visualize Fst
We can visualize Fst values across the genome using a genome track plot.   
In R
```
library(ggbio)
library(GenomicRanges)

# Read data
fst_data <- read.table("results/fst/wild_cultivated.windowed.weir.fst", header=T)

# Create genomic ranges
gr <- GRanges(
seqnames = fst_data$CHROM,
ranges = IRanges(start=fst_data$BIN_START, end=fst_data$BIN_END),
FST = fst_data$WEIGHTED_FST
)

# Generate genome track plot
png("figures/fst_genome_track.png", width=15, height=4, units="in", res=300)
autoplot(gr, aes(y=FST)) +
geom_line(color="#7570b3", size=0.3) +
facet_grid(. ~ seqnames, scales="free_x") +
labs(title="Wild vs Cultivated Fst") +
theme_bw() +
ylim(0, 0.5) # Set Fst axis range
dev.off()
```

