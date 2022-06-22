library(flipped)
library(ggplot2)

nheads=2
nflips=10
congruent <- find_congruent_models(nheads=nheads, nflips=nflips, slopes=0)
congruent$probabilities_per_flip$Model <- gsub("Slope of 0", "Constant Binomial", congruent$probabilities_per_flip$Model)
pdf(file="flipped_models.pdf", width=5, height=5)
print(ggplot(congruent$probabilities_per_flip, aes(x=Flips, y=Probability_of_heads_this_flip, colour=Model, group=Model)) + geom_line() + geom_point() + ylab("Probability of heads") + scale_x_continuous(name="Flip", breaks=seq(from=0, to=nflips, by=1))  + scale_colour_manual(values=c("black", "red", "blue")) + theme_light() + theme(legend.position = c(0.8, 0.8)))
dev.off()
system("open flipped_models.pdf")
