#!/usr/bin/env Rscript

require(ggplot2)
require(doBy)
require(RColorBrewer)
require(scales)
require(plyr)

args = commandArgs(trailingOnly = TRUE)
class_in = args[1]
out_plot = args[2]

s = read.table(class_in, comment.char = '', quote = '', sep = '\t')
colnames(s) = c('Sample', 'Condition', 'Timepoint', 'Relative Abundance', 'Reads', 'Direct Reads', 'Taxid', 'Tax')

#remove human
s = s[s$Taxid != 'Homo_sapiens',]
s = s[s$Taxid != 'Homo',]

s$`Relative Abundance` = ddply(s, c('Sample', 'Condition', 'Timepoint', 'Taxid'),
	function(x) x$Reads/sum(s$Reads[s$Sample == x[1,1]]))[,5]

toplot = s
abundance.threshold = sort(toplot$`Relative Abundance`, decreasing = T)[11]

toplot = cbind(toplot, toplot$`Relative Abundance` > abundance.threshold)
toplot$Taxid = factor(toplot$Taxid)
toplot$Taxid = factor(toplot$Taxid, levels = unique(toplot$Taxid[order(toplot$`Relative Abundance`, decreasing = T)]))

colors = c(brewer.pal(12, 'Set3'), rep(c('#C0C0C0', '#DCDCDC'), length(toplot$Taxid)))
colors[9] = '#CD6155'
ggplot(toplot) + geom_bar(aes(x = Condition, y = `Relative Abundance`, fill = Taxid), stat = 'identity') +
	scale_fill_manual(values = colors, breaks = levels(toplot$Taxid)[1:12], name = '') +
	scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
	facet_grid(.~Timepoint) + ylab('Relative Abundance') +
	theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(out_plot, width = 8.5, height = 4)
