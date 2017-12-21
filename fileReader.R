cellFile <- read.csv(cellLocation,sep="=",header=FALSE,stringsAsFactors = FALSE)
for (numVar in 1:length(cellFile$V1)){
  value <- (unlist(strsplit(cellFile$V2[numVar],"\t"))[1]); value <- (unlist(strsplit(value,"`"))[1])
  if (class(parse(text=value)) == "expression"){
    value <- eval(parse(text=value))
  }
  assign(cellFile[numVar,1],value)
}

enzymeFile <- read.csv(enzymeLocation,sep="=",header=FALSE,stringsAsFactors = FALSE)
for (numVar in 1:length(enzymeFile$V1)){
  value <- (unlist(strsplit(enzymeFile$V2[numVar],"\t"))[1]);value <- (unlist(strsplit(value,"`"))[1])
  if (class(parse(text=value)) == "expression"){
    value <- eval(parse(text=value))
  }
  assign(enzymeFile[numVar,1],value)
}

dosageFile <- read.csv(dosageLocation,sep="=",header=FALSE,stringsAsFactors = FALSE)
for (numVar in 1:length(dosageFile$V1)){
  value <- (unlist(strsplit(dosageFile$V2[numVar],"\t"))[1]);value <- (unlist(strsplit(value,"`"))[1])
  if (class(parse(text=value)) == "expression"){
    value <- eval(parse(text=value))
  }
  assign(dosageFile[numVar,1],value)
}
