# ------------ themes.r ------------------
# a file containing the figure themes used in R in the lab, 
# to maintain consistent visuals in papers and presentations
# ------------------------------

# The themes -----

# The size of the font is relative to the size of the figure. You can change it here.
fontsize <- 14

# theme for pie chart
paper_figs_pie_theme <- 
  theme_bw()+
  theme(plot.title=element_text(size=14, face="bold", hjust = 0.5, color = "#666666"),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "cm", data = NULL),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(), 
        legend.text=element_text(size=fontsize, color='black'),
        legend.title=element_text(size=fontsize, color='black'))

paper_figs_theme <- 
  theme_bw()+
  theme(plot.title=element_text(size=14, face="bold", hjust = 0.5, color = "#666666"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=fontsize, color='black'),
        axis.title = element_text(size=fontsize, color='black'),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=fontsize-2, color='black'),
        legend.title=element_text(size=fontsize, color='black'))


paper_figs_theme_no_legend <- 
  paper_figs_theme +
  theme(legend.position = 'none')
