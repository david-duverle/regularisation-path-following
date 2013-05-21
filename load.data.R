source("itemset.coxpath.R");

# first load data using one of the example below, then run itemset.coxpath() or run.cv()

### Lymphoma:

load.data.lymphoma.2002 <- function(location = 'data/lymphoma.2002/', filter.na = F, average.gene.exp = T) {
	survival = read.csv(paste(location, "survival.csv", sep=''), header = T)

	gene.exps = t(as.matrix(read.csv(paste(location, "exprs.preprocessed.csv", sep=''), header=FALSE)))
	

	platform = read.csv(paste(location, "geneName.csv", sep=''), header = T)
	gene.symbols = cbind(platform[,1], apply(platform, 1, function(line) { elems = strsplit(as.character(line[2]), "|", fixed=T)[[1]]; if(elems[2] != '') elems[2] else elems[4]; }))
	colnames(gene.symbols) = c("ID", "Gene.Symbol")
	
	empty.sym = (gene.symbols[,'Gene.Symbol'] == '')
	gene.symbols[empty.sym, 'Gene.Symbol'] = paste('pb', gene.symbols[empty.sym, 'ID'], sep='')
	
	if(average.gene.exp) {
		unique.gene.symbols = unique(gene.symbols[,'Gene.Symbol']);
		unique.gene.exps = matrix(0, nrow(gene.exps), length(unique.gene.symbols))
		for(i in 1:length(unique.gene.symbols)) {
			lines = which(gene.symbols[, 'Gene.Symbol'] == unique.gene.symbols[i])
			unique.gene.exps[, i] = apply(gene.exps[,lines, drop=F], 1, mean)
		}
		gene.exps = unique.gene.exps
		gene.symbols = cbind(1:length(unique.gene.symbols), unique.gene.symbols)
		colnames(gene.symbols) = c("ID", "Gene.Symbol")
	}
	gene.symbols = as.data.frame(gene.symbols)
	
	data.lymphoma.2002 = data.frame(x = I(gene.exps), time = as.numeric(survival$survival), status = as.numeric(survival$event))
	
	return (list(data.lymphoma.2002, gene.symbols))
}


### Breast cancer:


load.data.bc.2002 <- function(location = 'data/bc.2002/', filter.na = F,  average.gene.exp = T) {

	survival = read.csv(paste(location, "survival.csv", sep=''), header = T)

	gene.exps = t(as.matrix(read.csv(paste(location, "exprs.preprocessed.csv", sep=''), header=FALSE)))
	
	platform = read.csv(paste(location, "geneName.csv", sep=''), header = F, col.names = c("ID", "Gene.Symbol"))
	gene.symbols = cbind(platform[,1], apply(platform, 1, function(line) { if(! is.na(line[2]) && line[2] != '') line[2] else line[1]; }))
	colnames(gene.symbols) = c("ID", "Gene.Symbol")
	
	if(average.gene.exp) {
		unique.gene.symbols = unique(gene.symbols[,'Gene.Symbol']);
		unique.gene.exps = matrix(0, nrow(gene.exps), length(unique.gene.symbols))
		for(i in 1:length(unique.gene.symbols)) {
			lines = which(gene.symbols[, 'Gene.Symbol'] == unique.gene.symbols[i])
			unique.gene.exps[, i] = apply(gene.exps[,lines, drop=F], 1, mean)
		}
		gene.exps = unique.gene.exps
		gene.symbols = cbind(1:length(unique.gene.symbols), unique.gene.symbols)
		colnames(gene.symbols) = c("ID", "Gene.Symbol")
	}

	gene.symbols = as.data.frame(gene.symbols)

	data.bc.2002 = data.frame(x = I(gene.exps), time = as.numeric(survival$TIMEsurvival), status = as.numeric(survival$EVENTdeath))


	return(list(data.bc.2002, gene.symbols))
}

### Data05:

load.data.05 <- function(location = 'data/data05/', load.age.data = T, filter.na = F,  average.gene.exp = F) {
	
	gene.exps = t(as.matrix(read.csv(paste(location, "exprs.csv", sep=''), header=FALSE)))
	
	platform = read.table(paste(location, "GPL1875.annot", sep=''), header=TRUE, skip=22, sep="\t", quote = "", comment.char="")
	
	gene.symbols = platform[,c("GenBank.Accession", "Gene.symbol")]
	colnames(gene.symbols) = c("ID", "Gene.Symbol")
	gene.symbols$ID = as.character(gene.symbols$ID)
	gene.symbols$ID[which(gene.symbols$ID == "")] = paste("?", which(gene.symbols$ID == ""), sep="")
	# temp = as.character(gene.symbols$Gene.Symbol)
	# temp[temp == ""] = paste("pb", gene.symbols$ID[temp == ""], sep="")
	# gene.symbols$Gene.Symbol = temp
	
	if(average.gene.exp) {
		unique.gene.symbols = unique(gene.symbols[,'Gene.Symbol']);
		unique.gene.exps = matrix(0, nrow(gene.exps), length(unique.gene.symbols))
		for(i in 1:length(unique.gene.symbols)) {
			lines = which(gene.symbols[, 'Gene.Symbol'] == unique.gene.symbols[i])
			unique.gene.exps[, i] = apply(gene.exps[,lines, drop=F], 1, function(x) { mean(x, na.rm=T) })
		}
		gene.exps = unique.gene.exps
		gene.symbols = cbind(1:length(unique.gene.symbols), as.data.frame(unique.gene.symbols))
		colnames(gene.symbols) = c("ID", "Gene.Symbol")
	}

	gene.symbols = as.data.frame(gene.symbols)

	data.05 = data.frame(x = I(gene.exps), time = as.numeric(read.csv(paste(location, "survival.csv", sep=''), header = FALSE)$V1), status = as.numeric(read.csv(paste(location, "event.csv", sep=''), header = FALSE)$V1))
	
	if(load.age.data) {
		data.05$extra.features = cbind(read.table(paste(location, "withOrWithoutMYCNAmplification.csv", sep='')) > 0, read.table(paste(location, "age.csv", sep='')) >= 120, read.table(paste(location, "age.csv", sep='')) < 12)
		colnames(data.05$extra.features) = c("MYCN-amp", "over-10yo", "under-1yo")
	}
	else {
		data.05$extra.features = cbind(read.table(paste(location, "withOrWithoutMYCNAmplification.csv", sep='')) > 0)
		colnames(data.05$extra.features) = c("MYCN-amp")	
	}

	if(filter.na) {
		na.filtered = !probe.contains.na(data.05$x)
		data.05.na.filtered = list(x = data.05$x[, na.filtered], time = data.05$time, status = data.05$status, extra.features = data.05$extra.features)
		gene.symbols = gene.symbols[na.filtered, , drop = F]
		
		return(list(data.05.na.filtered, gene.symbols))
	}
	else
		return(list(data.05, gene.symbols))
}

# source('coxpath.R'); train.test.05 = training.testing.coxpath(data.05, max.steps=200, trace=T, deviation.threshold = 1.5);

### Data02:

load.data.02 <- function(location = 'data/data02/',  average.gene.exp = T) {

	gene.exps = t(as.matrix(read.csv(paste(location, "exprs.csv", sep=''), header=FALSE)))
	
	platform = read.csv(paste(location, "geneName.csv", sep=''), header=F)

	gene.symbols = cbind(paste("?", seq(length(platform$V1)), sep=""), as.character(platform$V1))
	colnames(gene.symbols) = c("ID", "Gene.Symbol")
	# gene.symbols[gene.symbols[,'Gene.Symbol'] == '', 'Gene.Symbol'] = gene.symbols[gene.symbols[,'Gene.Symbol'] == '', 'ID']
	
	if(average.gene.exp) {
		unique.gene.symbols = unique(gene.symbols[,'Gene.Symbol']);
		unique.gene.exps = matrix(0, nrow(gene.exps), length(unique.gene.symbols))
		for(i in 1:length(unique.gene.symbols)) {
			lines = which(gene.symbols[, 'Gene.Symbol'] == unique.gene.symbols[i])
			unique.gene.exps[, i] = apply(gene.exps[,lines, drop=F], 1, mean)
		}
		gene.exps = unique.gene.exps
		gene.symbols = cbind(1:length(unique.gene.symbols), unique.gene.symbols)
		colnames(gene.symbols) = c("ID", "Gene.Symbol")
	}

	gene.symbols = as.data.frame(gene.symbols)

	data.02 = data.frame(x = I(gene.exps), time = as.numeric(read.csv(paste(location, "survivalB.csv", sep=''), header = FALSE)$V1), status = as.numeric(read.csv(paste(location, "eventB.csv", sep=''), header = FALSE)$V1))

	data.02$extra.features = cbind(read.table(paste(location, "withOrWithoutMYCNAmplification.csv", sep='')) > 0, read.table(paste(location, "age.csv", sep='')) > 3650, read.table(paste(location, "age.csv", sep='')) <= 365)
	colnames(data.02$extra.features) = c("MYCN-amp", "over-10yo", "under-1yo")
	data.02$extra.features[is.na(data.02$extra.features)] = FALSE

	return(list(data.02, gene.symbols))
}
