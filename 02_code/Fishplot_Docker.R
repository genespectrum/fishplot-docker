library(fishplot)
library(clonevol)
library(sciClone)

clonevol.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')
input <- read.csv(file = "01_data/AML_07.txt", header = TRUE, dec = ".", stringsAsFactors = FALSE, sep = ",") 
input_names <- c("Time_1", "Time_2", "Time_3", "Time_4", "Time_5")

colnames_input <- colnames(input)
#input_number <- sum(startsWith(colnames_input, "var"), na.rm = TRUE)
input_number <- length(input_names)

# ==============================================================================
# Splitting Input-File
# ==============================================================================
#chr, pos, ref_reads, var_reads, vaf

tmp <- matrix(ncol = 5, nrow = nrow(input))
tmp <- as.data.frame(tmp)
colnames(tmp) <- c("chr", "pos", "ref_reads", "var_reads", "vaf")
input_annotation <- matrix(ncol = 3, nrow = nrow(input))
colnames(input_annotation) <- c("chr", "pos", "gene_name")

tmp[,1] <- input[,1]
tmp[,2] <- as.numeric(input[,2])

for (i in 1:input_number)
{
    tmp[,3] <- as.numeric(input[[paste0("ref_",i)]])
    tmp[,4] <- as.numeric(input[[paste0("var_",i)]])
    tmp[,5] <- as.numeric(input[[paste0("vaf_",i)]])
    assign(paste0("sample_",i), tmp)
}

input_annotation[,1] <- input[,1]
input_annotation[,2] <- input[,2]
input_annotation[,3] <- input[,ncol(input)]


# ==============================================================================
# SciClone Plots
# ==============================================================================

sample_list <- vector("list", input_number)

for (i in 1:input_number)
{
    sample_list[[i]] <- get(paste0("sample_",i))
}

sc = sciClone(
    vafs = sample_list,
    sampleNames = input_names,
    annotation = input_annotation,
    useSexChrs = FALSE,
    doClusteringAlongMargins = FALSE
)

sc.plot1d(sc, "03_output/sciclone_cluster1d.pdf")
sc.plot2d(sc, "03_output/sciclone_cluster2d.pdf")  
writeClusterTable(sc, "03_output/sciclone_clusterinfo.tsv")
sciclone_out <- read.table(file = "03_output/sciclone_clusterinfo.tsv", header = TRUE, sep = "\t")

# ==============================================================================
# ClonEvol
# ==============================================================================

vaf.col.names <- grep('.vaf', colnames(sciclone_out), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
sciclone_out[, sample.names] <- sciclone_out[, vaf.col.names] 
vaf.col.names <- sample.names

row.has.na <- apply(sciclone_out, 1, function(x){any(is.na(x))})
sciclone_out <- sciclone_out[!row.has.na,]
sciclone_out <- sciclone_out[order(sciclone_out$cluster),]

clonevol_input <- matrix(ncol = (2 + input_number), nrow = nrow(sciclone_out))
clonevol_input <- as.data.frame(clonevol_input, stringsAsFactors = FALSE)
colnames(clonevol_input) <- c("cluster", input_names, "gene")
clonevol_input[,1] <- as.numeric(sciclone_out$cluster)
for (i in 1:input_number)
{
    clonevol_input[,as.numeric(i + 1)] <- as.numeric(sciclone_out[,which(colnames(sciclone_out) == input_names[i])])
}
clonevol_input[,ncol(clonevol_input)] <- as.character(sciclone_out$gene_name)

clonevol_model <- infer.clonal.models(
    variants = clonevol_input,
    cluster.col.name = "cluster",
    model = "monoclonal",
    vaf.col.names = vaf.col.names,
    subclonal.test = "bootstrap",
    subclonal.test.model = "non-parametric",
    sample.groups = NULL,
    cluster.center = "mean",
    num.boots = 1000,
    founding.cluster = '1',
    min.cluster.vaf = 0.0001,
    p.value.cutoff = 0.0001,
    # alpha level in confidence interval estimate for CCF(clone)
    alpha = 0.01,
    random.seed = 63108,
    clone.colors = clonevol.colors,
    score.model.by = "probability",
    vaf.in.percent = TRUE
)

clonevol_model <- transfer.events.to.consensus.trees(
    clonevol_model,
    clonevol_input,
    cluster.col.name = 'cluster',
    event.col.name = 'gene')

clonevol_model <- convert.consensus.tree.clone.to.branch(
    clonevol_model, 
    branch.scale = 'sqrt')

plot.clonal.models(
    clonevol_model,
    # box plot parameters
    box.plot = TRUE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'NULL',
    fancy.variant.boxplot.highlight.shape = 21,
    fancy.variant.boxplot.highlight.fill.color = 'red',
    fancy.variant.boxplot.highlight.color = 'black',
    #fancy.variant.boxplot.highlight.note.col.name = 'is.driver',
    fancy.variant.boxplot.highlight.note.color = 'blue',
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.jitter.size = 4,
    fancy.variant.boxplot.base_size = 12,
    #size letter axis boxplot
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    
    # bell plot parameters
    clone.shape = 'bell',
    bell.event = TRUE,
    bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    
    
    # node-based consensus tree parameters
    merged.tree.plot = TRUE,
    tree.node.label.split.character = ',',
    tree.node.shape = 'circle',
    tree.node.size = 40,
    tree.node.text.size = 1,
    merged.tree.node.size.scale = 2,
    merged.tree.node.text.size.scale = 1,
    merged.tree.cell.frac.ci = FALSE,
    
    # branch-based consensus tree parameters
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 0.7,
    #size of the genes
    mtcab.branch.width = 0.5,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 0.5,
    #title in branched tree
    
    # cellular population parameters
    cell.plot = TRUE,
    num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'random',
    #horizontal, vertical, random
    
    #meta-parameters
    scale.monoclonal.cell.frac = FALSE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = TRUE,
    max.num.models.to.plot = 100,
    # output figure parameters
    out.dir = "03_output",
    out.format = 'pdf',
    overwrite.output = TRUE,
    width = 24,
    height = 16,
    
    # vector of width scales for each panel from left to right
    panel.widths = c(2, 1, 0.6, 1, 1.5)
)


# ==============================================================================
# FishPlot Plots
# ==============================================================================

fishplot_input <- generateFishplotInputs(clonevol_model, rescale = TRUE)

vafs = clonevol:::estimate.clone.vaf(clonevol_model$variants, vaf.col.names=vaf.col.names)
scales = vafs[vafs$cluster == 1, vaf.col.names]/max(vafs[vafs$cluster == 1, vaf.col.names])
for (i in 1:length(fishplot_input$cell.fractions)){
fishplot_input$cell.fractions[[i]] = (fishplot_input$cell.fractions[[i]] *
as.matrix(scales[rep(1,nrow(fishplot_input$cell.fractions[[i]])),]))
}

pdf("03_output/fishplot.pdf", width = 8, height = 4)
for (i in 1:length(fishplot_input$cell.fractions))
{
    
    fishplot_objects <- createFishObject(fishplot_input$cell.fractions[[i]], fishplot_input$parents[[i]])
    fishplot_layout = layoutClones(fishplot_objects)
    fishplot_layout = setCol(fishplot_layout, fishplot_input$clonevol.clone.colors)
    fishPlot(
        fishplot_layout,
        shape = "spline",
        cex.title = 0.7,
        vlines = seq(1, length(input_names)),
        vlab = input_names,
        pad.left = 0.5
    )
}
dev <- dev.off()

