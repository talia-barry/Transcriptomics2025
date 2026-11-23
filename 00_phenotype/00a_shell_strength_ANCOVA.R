# ===================================================================
# 00a_shell_strength_ANCOVA.R
# Purpose: Shell strength of the 5-month-old oysters (2.5 mm start size)
#          used for transcriptomics: picks the best size covariate (AIC),
#          runs the ANCOVA (Table 1) and plots strength / whole dry
#          weight by treatment (Figure 4).
# Input:   data/phenotype/shell_strength_5month_2.5mm.csv
#          (the 180 oysters from the Chapter 3 master data sheet,
#           "Inducible Defenses End 8-26 Master.csv", with
#           Age.Months == 5 and Start.Size..mm. == 2.5)
# Output:  results/00_phenotype/
# Run from the repository root (open Transcriptomics2025.Rproj).
# ===================================================================
dir.create("results/00_phenotype", recursive = TRUE, showWarnings = FALSE)


ID_5_2.5 <- read.csv("data/phenotype/shell_strength_5month_2.5mm.csv", header = TRUE)
dim(ID_5_2.5)   # 180 oysters (90 control, 90 treatment)

ID_5_2.5$treatment        <- as.factor(ID_5_2.5$treatment)
ID_5_2.5$strength.right   <- as.numeric(ID_5_2.5$strength.right)
ID_5_2.5$dry.weight.whole <- as.numeric(ID_5_2.5$dry.weight.whole)

# ---- Which size variable best explains strength? (lowest AIC wins) ----
m1 <- lm(strength.right ~ height,           data = ID_5_2.5)
m2 <- lm(strength.right ~ dry.weight.whole, data = ID_5_2.5)
m3 <- lm(strength.right ~ length,           data = ID_5_2.5)
m4 <- lm(strength.right ~ dry.weight.right, data = ID_5_2.5)
AIC(m1, m2, m3, m4)

# Average size at crushing
mean(ID_5_2.5$height)

# ---- ANCOVA: strength ~ treatment * whole dry weight (Table 1) ----
ancova_5_2.5 <- lm(strength.right ~ treatment * dry.weight.whole, data = ID_5_2.5)
summary(ancova_5_2.5)
table1 <- anova(ancova_5_2.5)
print(table1)
write.csv(as.data.frame(table1), "results/00_phenotype/Table1_ANCOVA_strength.csv")

# ---- Figure 4: strength / whole dry weight by treatment (mean +/- SE) ----
strength.dryweightwhole <- ID_5_2.5$strength.right / ID_5_2.5$dry.weight.whole

means <- tapply(strength.dryweightwhole, ID_5_2.5$treatment, mean)
sds   <- tapply(strength.dryweightwhole, ID_5_2.5$treatment, sd)
ns    <- tapply(strength.dryweightwhole, ID_5_2.5$treatment, length)
ses   <- sds / sqrt(ns)
print(data.frame(mean = means, se = ses, n = ns))

png("results/00_phenotype/Figure4_strength_by_dry_weight.png", width = 1200, height = 1200, res = 200)
b <- barplot(means, ylim = c(min(pretty(means - ses)), max(pretty(means + ses))),
             xpd = FALSE, ylab = "strength by whole dry weight")
arrows(b, means + ses, b, means - ses, angle = 90, code = 3)
dev.off()
