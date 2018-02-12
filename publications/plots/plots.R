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

cachefile1 <- "imprints_create.rds"
cachefile2 <- "imprints_query.rds"

if (!file.exists(cachefile1)) {
	# parsing log format into data frame
	logfile <- grep("^type .* not supported", readLines("all2.log"), invert=T, value=T)
	endoff <- (1:length(logfile))[grepl("^end of run", logfile)]
	startoff <- c(1, endoff[1:(length(endoff)-1)]+1)
	impres <- data.frame()
	queryres <- data.frame()

	for (i in 1:length(startoff)) {
		chunk <- paste(logfile[(startoff[i]):(endoff[i])], collapse="\n")
		label <- str_split(chunk, " ")[[1]][[1]]

		meta <- str_match(chunk, "filesize (\\d+) type (\\w+) [^\\n]+ records (\\d+)")
		filesize <- as.integer(meta[, 2])
		typename <- meta[, 3]
		records <- as.integer(meta[, 4])

		res <- str_match_all(chunk, "imprints (\\w+) creation time=(\\d+), (\\d+) usec per [^\\n]+\n[^\\n]+imprints data_size=(\\d+)\\(bytes\\) typesize=(\\d+) imprints_size=(\\d+)\\(bytes\\)[^\\n]+bins=(\\d+) blocksize=(\\d+)\\(bytes\\) imprintsize=(\\d+) values_per_bloc=(\\d+)\n")[[1]]
		resdf <- as.data.frame(matrix(as.integer(res[, 3:11]),nrow=nrow(res)))
		
		names(resdf) <- c("time_creation_usec", "time_per_value_usec", "data_size_bytes", "type_size_bytes", "imprint_index_size_bytes", "bins", "block_size_bytes", "imprints_size_bits", "values_per_block")
		resdf$label <- label
		resdf$exp <- res[, 2]
		resdf$filesize <- filesize
		resdf$typename <- typename
		resdf$records <- records
		impres <- rbind(impres, resdf)
		res <- str_match_all(chunk, "\n([\\w.]+)\\s+query\\[(\\d+)\\]=\\s+(\\d+)\\s+selectivity=([0-9.]+)%\\s+bins = (\\d+)\\s+blocksize = (\\d+)\\(bytes\\)\\s+simd_imprints = (\\d+)\\(usec\\)\n")[[1]]

		resdf <- data.frame(exp=res[, 2], query=as.integer(res[, 3]), result_size=as.integer(res[, 4]), selectivity = as.numeric(res[, 5]), bins= as.integer(res[, 6]), blocksize=as.integer(res[, 7]), query_time= as.integer(res[, 8]))
		resdf$filesize <- filesize
		resdf$typename <- typename
		resdf$records <- records
		queryres <- rbind(queryres, resdf)
	}
	saveRDS(impres, cachefile1)
	saveRDS(queryres, cachefile2)

} else {	
	impres <- readRDS(cachefile1)
	queryres <- readRDS(cachefile2)
}

print(head(impres))
print(head(queryres))

impres <- impres %>% filter(time_creation_usec > 1000)

queryres <- queryres %>% filter(exp %in% unique(impres$label))
queryres$typelen <- as.integer(queryres$filesize/queryres$records)
queryres$values_per_block <- as.integer(queryres$blocksize/queryres$typelen)

# plots

scalebins1 <- impres  %>% select(exp, bins, time_per_value_usec) %>% group_by(exp, bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))


breaks <- sort(unique(scalebins1$bins))

pdf("scalebins1.pdf", width=10, height=5)
ggplot(scalebins1, aes(y=mean, x=bins, group=exp)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none") + xlab("Bins / Imprint size (Bits)") + ylab("Time per 1000 values (usec)") + geom_text_repel(data = scalebins1 %>% filter((exp == "scalar" & bins==64) | (exp == "simd" & bins==256)), aes(label = exp), size = 8,  nudge_x = .4, segment.color = NA, family="serif") + scale_x_continuous(breaks=breaks, labels=as.character(breaks))
dev.off()



typewidth1 <- impres %>% filter(exp=="simd", bins >=64) %>% select(type_size_bytes, bins, time_per_value_usec) %>% group_by(type_size_bytes, bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

print(typewidth1)

typewidth1$type_size_bytes <- as.character(typewidth1$type_size_bytes)

pdf("typewidth1.pdf", width=10, height=5)
ggplot(typewidth1, aes(y=mean, x=bins, group=type_size_bytes)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none")  + xlab("Bins / Imprint size (Bits)") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks=breaks, labels=as.character(breaks)) + geom_text_repel(data = typewidth1 %>% filter(bins==256), aes(label = type_size_bytes), size = 8,  nudge_x = .4, segment.color = NA, family="serif")
dev.off()



valuesperblock1 <- impres %>% filter(exp=="simd", bins == 256) %>% select(values_per_block, time_per_value_usec) %>% group_by(values_per_block) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

breaks <- sort(unique(valuesperblock1$values_per_block))


pdf("valuesperblock1.pdf", width=10, height=5)
ggplot(valuesperblock1, aes(y=mean, x=values_per_block)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none")  + xlab("Values per imprint") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks=breaks, labels=as.character(breaks))
dev.off()



block_size_bytes1 <- impres %>% filter(exp=="simd", bins == 256) %>% select(block_size_bytes, time_per_value_usec) %>% group_by(block_size_bytes) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))


breaks <- sort(unique(block_size_bytes1$block_size_bytes))


pdf("block_size_bytes1.pdf", width=10, height=5)
ggplot(block_size_bytes1, aes(y=mean, x=block_size_bytes)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme + theme(legend.position = "none")  + xlab("Block size (Bytes)") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks = breaks, labels = as.character(breaks)) 
dev.off()



# # querying

# stupid me, using different names for the same things above
queryres$label <- queryres$exp
queryres$exp <- NULL

impres$blocksize <- impres$block_size_bytes
impres$block_size_bytes <- NULL

# dplyr magic
queryres_joined <- left_join(queryres, impres %>% select(label, bins, blocksize, imprint_index_size_bytes, type_size_bytes), by=c("label", "bins", "blocksize")) %>% mutate(imprint_entries = (imprint_index_size_bytes/(bins/8))) %>% filter(imprint_entries > 1000) %>% mutate(query_time_1k_imprints = query_time / (imprint_entries/1000), values_per_block = blocksize/type_size_bytes)



# moar data into imprint -> crappier results but way better with larger imprints
querybins1 <- queryres_joined %>% filter(bins %in% c(8,256)) %>% select(bins, blocksize, query_time_1k_imprints) %>% group_by(bins, blocksize) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))


	breaks <- sort(unique(querybins1$blocksize))
	querybins1$bins <- as.character(querybins1$bins)

	pdf("qblocksize1.pdf", width=10, height=5)
	print(ggplot(querybins1, aes(y=mean, x=blocksize, group=bins)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme  + xlab("Block size (Bytes)") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks = breaks, labels = as.character(breaks), limits=c(64,265)) + theme(legend.position = "none") + geom_text(data = querybins1 %>% filter(blocksize==256), aes(x=260,label = bins), size = 8, hjust=0, family="serif")) 
	dev.off()




# for (bi in c(64,128,256)) {

# 	querybins1 <- queryres_joined %>% filter(blocksize==bi)  %>% select(bins, query_time_1k_imprints) %>% group_by(bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

# 	breaks <- sort(unique(querybins1$bins))


# 	pdf(paste0("queryblocksize-bi-",bi,".pdf"), width=10, height=5)
# 	print(ggplot(querybins1, aes(y=mean, x=bins)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme  + xlab("Bins") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks = breaks, labels = as.character(breaks)))

# 	 # + geom_text_repel(data = typewidth1 %>% filter(bins==256), aes(label = type_size_bytes), size = 8,  nudge_x = .4, segment.color = NA, family="serif")
# 	dev.off()


# }

	querybins1 <- queryres_joined %>% filter(blocksize %in% c(64, 256))  %>% select(blocksize, bins, query_time_1k_imprints) %>% group_by(blocksize, bins) %>% summarise_each(funs(mean,se=sd(.)/sqrt(n())))

	breaks <- sort(unique(querybins1$bins))
	querybins1$blocksize <- as.character(querybins1$blocksize)


	pdf("qblocksize2.pdf", width=10, height=5)
	print(ggplot(querybins1, aes(y=mean, x=bins, group=blocksize)) + geom_line(size=1.1) + geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=1, size=1) + theme  + xlab("Bins / Imprint size (Bits)") + ylab("Time per 1000 values (usec)") + scale_x_continuous(breaks = breaks, labels = as.character(breaks), limits=c(8,265)) + geom_text(data = querybins1 %>% filter(bins==256), aes(x=260,label = blocksize), size = 8, hjust=0, family="serif"))

	 # + geom_text_repel(data = typewidth1 %>% filter(bins==256), aes(label = type_size_bytes), size = 8,  nudge_x = .4, segment.color = NA, family="serif")
	dev.off()


