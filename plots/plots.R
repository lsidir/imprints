# on-demand package autoinstaller hack
(function(lp) {
np <- lp[!(lp %in% installed.packages()[,"Package"])]
if(length(np)) install.packages(np,repos=c("http://cran.rstudio.com/"))
x <- lapply(lp,function(x){library(x,character.only=TRUE)}) 
})(c("dplyr", "ggplot2", "ggthemes", "ggrepel", "stringr"))


# plot style setup
theme <- theme_few(base_size = 24) + 
theme(axis.title.y=element_text(vjust=0.9), 
  axis.title.x=element_text(vjust=-0.1),
  axis.ticks.x=element_blank(),
  text=element_text(family="serif"))

cachefile <- "all.log.rds"

if (!file.exists(cachefile)) {
	# parsing log format into data frame
	logfile <- grep("^type .* not supported", readLines("all.log"), invert=T, value=T)
	endoff <- (1:length(logfile))[grepl("^end of run", logfile)]
	startoff <- c(1, endoff[1:(length(endoff)-1)]+1)
	impres <- data.frame()
	for (i in 1:length(startoff)) {
		chunk <- paste(logfile[(startoff[i]):(endoff[i])], collapse="\n")
		label <- str_split(chunk, " ")[[1]][[1]]

		meta <- str_match(chunk, "filesize (\\d+) type (\\w+) [^\\n]+ records (\\d+)")
		filesize <- as.integer(meta[, 2])
		typename <- meta[, 3]
		records <- as.integer(meta[, 4])

		res <- str_match_all(chunk, "imprints (\\w+) creation time=(\\d+), (\\d+) usec per [^\\n]+\n[^\\n]+imprints data_size=(\\d+)\\(bytes\\) typesize=(\\d+)[^\\n]+bins=(\\d+) blocksize=(\\d+)\\(bytes\\) imprintsize=(\\d+) values_per_bloc=(\\d+)\n")[[1]]
		resdf <- as.data.frame(matrix(as.integer(res[, 3:10]),nrow=nrow(res)))
		names(resdf) <- c("time_creation_usec", "time_per_value_usec", "data_size_bytes", "type_size_bytes", "bins", "block_size_bytes", "imprints_size_bits", "values_per_block")
		resdf$label <- label
		resdf$exp <- res[, 2]
		resdf$filesize <- filesize
		resdf$typename <- typename
		resdf$records <- records

		impres <- rbind(impres, resdf)
	}
	saveRDS(impres, cachefile)
} else {	
	impres <- readRDS(cachefile)
}

impres <- impres %>% filter(time_creation_usec > 1000)

# plots

print(head(impres))

summary(impres$filesize)

scalebins1 <- impres  %>% select(exp, bins, time_per_value_usec) %>% group_by(exp, bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

print(scalebins1)

breaks <- sort(unique(scalebins1$bins))

pdf("scalebins1.pdf", width=10, height=5)
ggplot(scalebins1, aes(y=mean, x=bins, group=exp)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none") + xlab("Bins / Imprint size") + ylab("Time per 1000 values (usec)") + geom_text_repel(data = scalebins1 %>% filter((exp == "scalar" & bins==64) | (exp == "simd" & bins==256)), aes(label = exp), size = 8,  nudge_x = .4, segment.color = NA, family="serif") + scale_x_continuous(breaks=breaks, labels=as.character(breaks))
dev.off()



typewidth1 <- impres %>% filter(exp=="simd", bins >=64) %>% select(type_size_bytes, bins, time_per_value_usec) %>% group_by(type_size_bytes, bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

typewidth1$type_size_bytes <- as.character(typewidth1$type_size_bytes)
print(typewidth1)

pdf("typewidth1.pdf", width=10, height=5)
ggplot(typewidth1, aes(y=mean, x=bins, group=type_size_bytes)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none")  + xlab("Bins / Imprint size") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks=breaks, labels=as.character(breaks)) + geom_text_repel(data = typewidth1 %>% filter(bins==256), aes(label = type_size_bytes), size = 8,  nudge_x = .4, segment.color = NA, family="serif")
dev.off()



valuesperblock1 <- impres %>% filter(exp=="simd", bins == 256) %>% select(values_per_block, time_per_value_usec) %>% group_by(values_per_block) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

#valuesperblock1$values_per_block <- as.character(valuesperblock1$values_per_block)
print(valuesperblock1)

breaks <- sort(unique(valuesperblock1$values_per_block))


pdf("valuesperblock1.pdf", width=10, height=5)
ggplot(valuesperblock1, aes(y=mean, x=values_per_block)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none")  + xlab("Values per imprint") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks=breaks, labels=as.character(breaks)) #+ geom_text_repel(data = typewidth1 %>% filter(bins==256), aes(label = type_size_bytes), size = 8,  nudge_x = .4, segment.color = NA, family="serif")
dev.off()


# values per block






# + geom_line(size=1.5) + geom_point(size=3) + scale_y_log10() + scale_x_log10(breaks=c(0.1,1,10,100), limits=c(0.1, 500)) + theme  + xlab("Latency (ms, log)") + ylab("Wall clock time (s, log)") +  theme(legend.position = "none") +  geom_text_repel(data = dd5 %>% filter(is.na(throughput), latency==100 | (system=="DBMS X" & latency==10)), aes(label = system), size = 8,  nudge_x = .3, segment.color = NA, family="serif") 

