
### Author: Carl de Boer
### Note: Please don't share this for the time being.


library(rvest)
library(Biostrings)
library(stringr)
library(textclean)
library(ggplot2)
library(amap)
library(reshape)
library(stringi)
source("https://raw.githubusercontent.com/Carldeboer/SingleCellR/master/ConvenienceFunctions.R")


getCode = function(html){
  code = html %>% html_nodes('pre.r') %>% html_nodes('code') %>% html_text()
  return(code)
}

getOutput = function(html){
  out = html %>% html_nodes('pre:not(pre.r)') %>% 
    html_nodes('code') %>% html_text()
  return(out)
}

getTextContents = function(html){
  code = html %>% html_nodes('p, h1, h2, h3, h4, h5, h6, h7, th, td') %>% html_text()
  return(code)
}


getImgSrc = function(html){
  code = html %>% html_nodes("img")  %>% html_attr('src')
  return(code)
}

getTags = function(html){
  code = html %>% html_nodes("*") %>% html_name()
  return(code)
  #return(str_extract_all(html_text, "</?[a-zA-Z]+( ?[^>]*)>"))
}


noScale = function(x,y){
  return(1);
}
pairwiseAlignText = function(x,y, mismatch=0, match=1, scoreScale=mean, scoreOnly=T, ...){
  if(nchar(x)==0 && nchar(y)==0){
    return(NA);
  }
  if(nchar(x)>45000){
    x=substr(x,0,45000);
    warning("x is longer than 45,000 characters - truncating")
  }
  if(nchar(y)>45000){
    y=substr(y,0,45000);
    warning("y is longer than 45,000 characters - truncating")
  }
  
  allChars = unique(c(unlist(strsplit(x,"")), unlist(strsplit(y,""))))
  subMat = matrix(mismatch, nrow=length(allChars), ncol=length(allChars), dimnames=list(allChars, allChars))
  diag(subMat)=match;
  gAlign=pairwiseAlignment(x,y, substitutionMatrix = subMat, scoreOnly=scoreOnly, ...)
  if(scoreOnly){
    return(gAlign/scoreScale(nchar(x),nchar(y)))
  }else{
    gAlign@score=gAlign@score/scoreScale(nchar(x),nchar(y))
    return(gAlign)
  }
  
}

fractionIdentical = function(x,y, scoreScale=min){
  return(sum(x %in% y)/scoreScale(length(x),length(y)))
}

cleanText = function(x){
  if(length(x)==1 && x==""){
    return(x)
  }
  x= gsub("\r\n","\n",x)
  x=unlist(strsplit(x, split="\n"))
  #x=replace_non_ascii(x)
  x=stri_enc_toascii(replace_curly_quote(x))
  x=paste(x, collapse = "\n")
  return(x)
}

loadCorpus = function(corpusDF){
  stopifnot("file" %in% names(corpusDF));
  stopifnot("ID" %in% names(corpusDF));
  stopifnot("assignment" %in% names(corpusDF));
  corpus=list();
  for (assignment in unique(corpusDF$assignment)){
    allCode = list()
    allSCode = list()
    allOutput = list()
    allWriting = list()
    allTags = list();
    uTags = c();
    allImages = list();
    allComments = list();
    allIDs = unique(corpusDF$ID[corpusDF$assignment==assignment]);
    for(id in allIDs){
      message(id);
      doc1 = read_html(corpusDF$file[corpusDF$assignment==assignment & corpusDF$ID == id])
      code = paste(getCode(doc1), collapse="\n")
      code = gsub("\r\n","\n",code)
      parsedCode = getParseData(parse(text = code))
      
      justCode = paste(parsedCode$text[parsedCode$token!="COMMENT"],collapse=" ")
      justComments = paste(parsedCode$text[parsedCode$token=="COMMENT"],collapse="\n")
      strippedCode = parsedCode[parsedCode$token!="COMMENT",]
      strippedCode$text[strippedCode$token=="SYMBOL"]="SYMBOL"
      strippedCode$text[strippedCode$token=="STR_CONST"]="STRING"
      strippedCode$text[strippedCode$token=="LEFT_ASSIGN"]="="
      strippedCode=paste(strippedCode$text, collapse=" ")
      allCode[id] = cleanText(justCode)
      allSCode[id] = cleanText(strippedCode)
      allComments[id] = cleanText(justComments)
      allOutput[id]=cleanText(getOutput(doc1))
      allWriting[id]=cleanText(getTextContents(doc1))
      allImages[id] = list(getImgSrc(doc1))
      allTags[id] = list(getTags(doc1))
      uTags = unique(c(uTags,allTags[[id]]));
      #save the components separately
      write.table(code, file = sprintf("%s.code.txt",corpusDF$file[corpusDF$assignment==assignment & corpusDF$ID == id]), row.names = F, quote=F, col.names = F)
      write.table(strippedCode, file = sprintf("%s.strippedcode.txt",corpusDF$file[corpusDF$assignment==assignment & corpusDF$ID == id]), row.names = F, quote=F, col.names = F)
      write.table(allOutput[[id]], file = sprintf("%s.output.txt",corpusDF$file[corpusDF$assignment==assignment & corpusDF$ID == id]), row.names = F, quote=F, col.names = F)
      write.table(allWriting[[id]], file = sprintf("%s.writing.txt",corpusDF$file[corpusDF$assignment==assignment & corpusDF$ID == id]), row.names = F, quote=F, col.names = F)
    }
    
    #encode html tags
    uTagsMap = c(LETTERS,tolower(LETTERS))[1:length(uTags)]
    names(uTagsMap)=uTags;
    for(id in allIDs){
      allTags[[id]] = paste(uTagsMap[allTags[[id]]],collapse="")
    }
    corpus[[assignment]] = list(ids = allIDs, code=allCode, scode = allSCode, output=allOutput, writing=allWriting, tags=allTags, images=allImages, comments=allComments)
  }
  return(corpus);
}

getPairwiseSimilarity = function(corpus, scoreScale=noScale, ...){
  similarityScores = data.frame();
  for(assignment in names(corpus)){
    #calculate similarity for different measures
    for(i in 1:(length(corpus[[assignment]]$ids)-1)){
      id1=corpus[[assignment]]$ids[i];
      for (j in (i+1):length(corpus[[assignment]]$ids)){
        id2=corpus[[assignment]]$ids[j];
        message(sprintf("%s - %s",id1,id2))
        similarityScores = rbind(similarityScores, data.frame(
          ID1 = id1, ID2=id2, assignment=assignment,
          codeScore = pairwiseAlignText(corpus[[assignment]]$code[[id1]],corpus[[assignment]]$code[[id2]], scoreScale=scoreScale, ...),
          scodeScore = pairwiseAlignText(corpus[[assignment]]$scode[[id1]],corpus[[assignment]]$scode[[id2]], scoreScale=scoreScale, ...),
          commentsScore = pairwiseAlignText(corpus[[assignment]]$comments[[id1]],corpus[[assignment]]$comments[[id2]], scoreScale=scoreScale, ...),
          outputScore = pairwiseAlignText(corpus[[assignment]]$output[[id1]],corpus[[assignment]]$output[[id2]], scoreScale=scoreScale, ...),
          writingScore = pairwiseAlignText(corpus[[assignment]]$writing[[id1]],corpus[[assignment]]$writing[[id2]], scoreScale=scoreScale, ...),
          imageScore = fractionIdentical(corpus[[assignment]]$images[[id1]],corpus[[assignment]]$images[[id2]], scoreScale=scoreScale),
          tagsScore = pairwiseAlignText(corpus[[assignment]]$tags[[id1]],corpus[[assignment]]$tags[[id2]], scoreScale=scoreScale, ...)
        ))
      }
    }
  }
  return(similarityScores)
}

makeSimilarityReport = function(corpus, corpusSimilarity, assignment, id1, id2, filePre, grade="", comment="", ...){
  curSim = corpusSimilarity[corpusSimilarity$assignment==assignment & corpusSimilarity$ID1 %in% c(id1,id2) & corpusSimilarity$ID2 %in% c(id1,id2),]
  curSim = curSim[1,] # if I have duplicated this, then just take the first.
  fileConn<-file(sprintf("%s.%s.%s.%s.txt",filePre,id1,id2,assignment),open = "w")
  spacer="##########################################################";
  writeLines(c(spacer,sprintf("#Similarity report for %s and %s, assignment %s", id1,id2,assignment),spacer), fileConn)
  writeLines(comment, fileConn)
  imageScore = fractionIdentical(corpus[[assignment]]$images[[id1]],corpus[[assignment]]$images[[id2]], scoreScale=noScale);
  writeLines(sprintf("#Of %i and %i images, %i were identical (%g percentile)", length(corpus[[assignment]]$images[[id1]]), length(corpus[[assignment]]$images[[id1]]), imageScore, mean(curSim$imageScore>=corpusSimilarity$imageScore, na.rm=T)*100), fileConn)
  if (grade!=""){
    writeLines(sprintf("#Grades were %g correlated (%g percentile)", curSim[[grade]], mean(curSim[[grade]]>=corpusSimilarity[[grade]], na.rm=T)*100), fileConn)
  }
  for(scoreType in c("scode","code","output","writing","tags","comments")){
    temp = pairwiseAlignText(corpus[[assignment]][[scoreType]][[id1]], corpus[[assignment]][[scoreType]][[id2]], scoreOnly=F, ...)
    #r = BStringSet( c( toString( subject(temp)), toString(pattern(temp))))
    writeLines(c("","",spacer,sprintf("#Alignment for %s, where score = %f (%g percentile)", scoreType, temp@score, mean(curSim[[sprintf("%sScore",scoreType)]]>=corpusSimilarity[[sprintf("%sScore",scoreType)]], na.rm=T)*100),spacer), fileConn)
    if (nchar(corpus[[assignment]][[scoreType]][[id1]])==0 || nchar(corpus[[assignment]][[scoreType]][[id2]])==0){
      writeLines("Nothing to align. One entry is blank.", fileConn)
    }else{
      tryCatch({writePairwiseAlignments(temp, file=fileConn, Matrix=NA, block.width=50)},
             error=function(X){writeLines("Unable to print alignment- bug in someone else's code", fileConn)})
    }
    
  }
  close(fileConn)
}

combineZStouffer = function(x){sum(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))}



similaritiesToZs = function(corpusSimilarity, meanFunc = median, sdFunc=mad){
  corpusZs= corpusSimilarity;
  for(assignment in unique(corpusSimilarity$assignment)){
    for(curScore in names(corpusSimilarity)[!(names(corpusSimilarity) %in% c("ID1","ID2","assignment"))]){
      corpusZs[corpusZs$assignment==assignment, curScore]= qnorm((rank(corpusSimilarity[corpusZs$assignment==assignment, curScore], na.last="keep")/sum(corpusZs$assignment==assignment)))
      corpusZs[is.infinite(corpusZs[[curScore]]),curScore]=NA;
      corpusZs[corpusZs$assignment==assignment, curScore]= corpusZs[corpusZs$assignment==assignment, curScore]-meanFunc(corpusZs[corpusZs$assignment==assignment, curScore], na.rm=T)
      corpusZs[corpusZs$assignment==assignment, curScore]=corpusZs[corpusZs$assignment==assignment, curScore]/sdFunc(corpusZs[corpusZs$assignment==assignment, curScore], na.rm=T)
    }
  }
  return(corpusZs)
}

####Analysis


fileLocation = "C:\\Users\\cdeboer\\VirtualBox VMs\\GenomeInformatics\\Shared\\A1_compiled"
outDir = sprintf("%s/../analysis_no_grades_for_A2/", fileLocation)
dir.create(outDir)

#A1
corpusDF = read.table(sprintf("%s/corpus.A1.txt",fileLocation), header = T)
corpusDF$file = sprintf("%s/%s",fileLocation, corpusDF$file)
#A2
corpusDF2 = data.frame(file=list.files(sprintf("%s/../A2_compiled/",fileLocation)))
corpusDF2$ID = gsub(".html$","",corpusDF2$file);
corpusDF2$assignment="A2";
corpusDF2$file = sprintf(sprintf("%s/../A2_compiled/%s",fileLocation, corpusDF2$file))
corpusDF2 = corpusDF2[grepl(".html$",corpusDF2$file),]

#grades for A1 and student info
studentInfo = read.table(sprintf("%s/studentInfo.txt",fileLocation), header = T, sep="\t")
studentGrades = read.table(sprintf("%s/studentGrades.txt",fileLocation), header = T, sep="\t")
studentGrades = merge(studentInfo[c("GitHub.UN","SIS.User.ID")], studentGrades, by.x="SIS.User.ID", by.y="Student.ID", all=T)
studentGrades$SIS.User.ID[is.na(studentGrades$Q1.1..3.)]
studentGrades[is.na(studentGrades$GitHub.UN),]
studentGrades = studentGrades[studentGrades$GitHub.UN!="",]
row.names(studentGrades)=studentGrades$GitHub.UN
studentGrades$SIS.User.ID=NULL; studentGrades$GitHub.UN=NULL;
studentGrades[is.na(studentGrades)]=0
maxPerSection =unlist(apply(studentGrades, 2, max))
studentGrades = scale(studentGrades,center=F, scale=as.numeric(maxPerSection))
gradeCor = cor(t(as.matrix(studentGrades)))
meltedGradeCor = melt(gradeCor)
names(meltedGradeCor) = c("ID1","ID2","r_grade")
p=ggplot(meltedGradeCor, aes(x=r_grade))+geom_density(); print(p)

#load the data
corpus = loadCorpus(corpusDF)
corpusA2 = loadCorpus(corpusDF2)

#compare the assignments
corpusSimilarityNoPenalty = getPairwiseSimilarity(corpus, gapOpen=0, gapExtension=0, match=1, mismatch=0, scoreScale=min)
corpusSimilarityNoPenaltyA2 = getPairwiseSimilarity(corpusA2, gapOpen=0, gapExtension=0, match=1, mismatch=0, scoreScale=min)


#merge with grade similarity
corpusSimilarityNoPenalty = merge(corpusSimilarityNoPenalty, meltedGradeCor, by=c("ID1","ID2"), all.x=T)
corpusSimilarityNoPenaltyA2$r_grade=NA;

corpusSimilarityBoth = rbind(corpusSimilarityNoPenalty, corpusSimilarityNoPenaltyA2)

#duplicate so that X-Y and Y-X are both present
corpusSimilarityBoth2=corpusSimilarityBoth;
corpusSimilarityBoth2$ID3 = corpusSimilarityBoth2$ID1; 
corpusSimilarityBoth2$ID1 = corpusSimilarityBoth2$ID2;
corpusSimilarityBoth2$ID2 = corpusSimilarityBoth2$ID3;
corpusSimilarityBoth2$ID3=NULL;
corpusSimilarityBoth = rbind(corpusSimilarityBoth, corpusSimilarityBoth2)


#cluster data by similarity, both assignments at once
idOrderImages = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="imageScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderImages)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderImages)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=imageScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_image_identity.pdf",outDir), width=19, height=9)

idOrderOutput = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="outputScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderOutput)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderOutput)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=outputScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_output_identity.pdf",outDir), width=19, height=9)


idOrderCode = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="codeScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderCode)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderCode)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=codeScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_code_identity.pdf",outDir), width=19, height=9)

idOrderSCode = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="scodeScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderSCode)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderSCode)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=scodeScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_strippedcode_identity.pdf",outDir), width=19, height=9)

idOrderWriting = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="writingScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderWriting)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderWriting)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=writingScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_writing_identity.pdf",outDir), width=19, height=9)

idOrderTags = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="tagsScore", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderTags)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderTags)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=tagsScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_HTML_tag_identity.pdf",outDir), width=19, height=9)

idOrderGrades = clusterDataFrame(corpusSimilarityBoth, ID1 ~ID2+assignment, value="r_grade", NA.set = 0)
corpusSimilarityBoth$ID1 = factor(as.character(corpusSimilarityBoth$ID1), levels=idOrderGrades)
corpusSimilarityBoth$ID2 = factor(as.character(corpusSimilarityBoth$ID2), levels=idOrderGrades)
p=ggplot(corpusSimilarityBoth, aes(x=ID1, y=ID2, fill=r_grade)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
ggsave(sprintf("%s/clustered_by_grade_r.pdf",outDir), width=19, height=9)


#separateky
for(assignment in c("A1","A2")){
  corpusSimilarityCur = corpusSimilarityBoth[corpusSimilarityBoth$assignment==assignment,]
  #cluster data by similarity, both assignments at once
  idOrderImages = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="imageScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderImages)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderImages)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=imageScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_image_identity.pdf",outDir, assignment), width=10, height=9)
  
  idOrderOutput = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="outputScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderOutput)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderOutput)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=outputScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_output_identity.pdf",outDir, assignment), width=10, height=9)
  
  idOrderCode = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="codeScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderCode)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderCode)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=codeScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_code_identity.pdf",outDir, assignment), width=10, height=9)
  
  idOrderSCode = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="scodeScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderSCode)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderSCode)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=scodeScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_strippedcode_identity.pdf",outDir, assignment), width=10, height=9)
  
  idOrderWriting = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="writingScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderWriting)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderWriting)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=writingScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_writing_identity.pdf",outDir, assignment), width=10, height=9)
  
  idOrderTags = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="tagsScore", NA.set = 0)
  corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderTags)
  corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderTags)
  p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=tagsScore)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
  ggsave(sprintf("%s/%s_clustered_by_HTML_tag_identity.pdf",outDir, assignment), width=10, height=9)
  
  if(assignment=="A1"){
    idOrderGrades = clusterDataFrame(corpusSimilarityCur, ID1 ~ID2+assignment, value="r_grade", NA.set = 0)
    corpusSimilarityCur$ID1 = factor(as.character(corpusSimilarityCur$ID1), levels=idOrderGrades)
    corpusSimilarityCur$ID2 = factor(as.character(corpusSimilarityCur$ID2), levels=idOrderGrades)
    p=ggplot(corpusSimilarityCur, aes(x=ID1, y=ID2, fill=r_grade)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust=1))+facet_grid(. ~ assignment); print(p);
    ggsave(sprintf("%s/%s_clustered_by_grade_r.pdf",outDir, assignment), width=10, height=9)
  }
}


corpusSimilarityBothZ = similaritiesToZs(corpusSimilarityBoth, meanFunc=median, sdFunc=sd)

corpusSimilarityBothZMelted = melt(corpusSimilarityBothZ, id.vars = c("ID1","ID2","assignment"))
names(corpusSimilarityBothZMelted)[(ncol(corpusSimilarityBothZMelted)-1):ncol(corpusSimilarityBothZMelted)]=c("criteria","Z")

p=ggplot(corpusSimilarityBothZMelted, aes(x=Z, colour=criteria)) + geom_density() + facet_grid(.~assignment)+theme_classic(); print(p)
ggsave(sprintf("%s/Similarity_Z_score_dist.pdf",outDir), width=1+3+3, height=3)


corpusSimilarityBothZMeltedFlagged = corpusSimilarityBothZMelted[naIsFalse(corpusSimilarityBothZMelted$Z>3),]
for(i in 1:nrow(corpusSimilarityBothZMeltedFlagged)){
  corpusSimilarityBothZMeltedFlagged[i,1:2]= sort(corpusSimilarityBothZMeltedFlagged[i,1:2])
}

reportDir = sprintf("%s/../analysis_no_grades_for_A2/reports/",fileLocation)
dir.create(reportDir)
corpusSimilarityBothZMeltedFlagged=unique(corpusSimilarityBothZMeltedFlagged)
corpusSimilarityBothZMeltedFlagged$reason = sprintf("Z(%s)=%g", corpusSimilarityBothZMeltedFlagged$criteria, corpusSimilarityBothZMeltedFlagged$Z)
flaggedPairs = cast(corpusSimilarityBothZMeltedFlagged, ID1+ID2 +assignment~., value="reason", fun.aggregate = function(x){paste(x, collapse="; ")})
names(flaggedPairs)[ncol(flaggedPairs)]="reason"
for(i in 1:nrow(flaggedPairs)){
  curReason=sprintf("%s and %s were flagged on %s because %s",as.character(flaggedPairs$ID1[i]), 
                        as.character(flaggedPairs$ID2[i]), flaggedPairs$assignment[i], flaggedPairs$reason[i])
  message(curReason)
  makeSimilarityReport(corpusBoth, corpusSimilarity = corpusSimilarityBoth, 
                       assignment=flaggedPairs$assignment[i],
                       id1=as.character(flaggedPairs$ID1[i]), id2=as.character(flaggedPairs$ID2[i]), 
                       file = sprintf("%s/simReport.test",reportDir), comment = curReason, grade = "r_grade", 
                       gapOpen=0, gapExtension=0, match=1, mismatch=0, scoreScale=noScale )
}

assignment=flaggedPairs$assignment[i]
id1=flaggedPairs$ID1[i]
id2=flaggedPairs$ID2[i]



temp = pairwiseAlignText(corpusBoth$A1$code$AlixSavard, corpusBoth$A1$code$ctcheune, gapOpen=0, gapExtension=0, match=1, mismatch=0, scoreScale=noScale, scoreOnly = F)
fileConn = file(sprintf("%s.temp.align",reportDir),open = "w")
writePairwiseAlignments(temp, file=fileConn, Matrix=NA, block.width=50)
close(fileConn)
makeSimilarityReport(corpusBoth, corpusSimilarity = corpusSimilarityBoth, assignment="A1",id1="arshaan9", id2="WalterHansWong", file = sprintf("%s/simReport.test",fileLocation), gapOpen=0, gapExtension=0, match=1, mismatch=0, scoreScale=noScale, grade="r_grade" )
