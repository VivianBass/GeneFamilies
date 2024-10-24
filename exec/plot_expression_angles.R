
require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)
library(RColorBrewer)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/")

load(file.path(output_data_dir, "exp.angle.diag.RData"))                       
 
# plotting the expression angle data
plot.df <- expr.angle.diag.df[!is.nan(expr.angle.diag.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = c("ortholog", "paralog"))
pdf(file.path(output_data_dir, "expressionAngleToDiagonalBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(angle.diag ~ gene.type, data = plot.df, xlab = "Type of Gene", 
    ylab = "relative tissue specificity", border = pushed.colors, col = addAlpha(pushed.colors))
dev.off()

plot.df$rel.vers <- 1 - plot.df$angle.diag
pdf(file.path(output_data_dir, "relativeExpressionVersatilityBoxplot.pdf"))
boxplot(rel.vers ~ gene.type, data = plot.df, xlab = "Type of Gene", ylab = "relative tissue versatility", 
    border = pushed.colors, col = addAlpha(pushed.colors), outline = FALSE)
dev.off()

message("DONE")
