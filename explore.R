library(trena)
library(trenaViz)
library(colorspace)
library(annotate)
library(org.Hs.eg.db)
library(MotifDb)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
# make sure we are executing in the right directory, where all the needed Allen Brain Atalas data can
# be found
stopifnot(length(grep("H0", list.files("./data/", all.files=TRUE))) == 12)
#----------------------------------------------------------------------------------------------------
stopifnot(packageVersion("trena")    >= "0.99.161")
stopifnot(packageVersion("trenaViz") >= "0.99.18")
stopifnot(packageVersion("httpuv")   >= "1.3.5")
#----------------------------------------------------------------------------------------------------
if(!exists("trena"))
   trena <- Trena("hg38")

PORT.RANGE <- 8000:8020

if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=FALSE)
   setGenome(tv, "hg38")
   }


#----------------------------------------------------------------------------------------------------
#
#  The main components of the basal ganglia – as defined functionally – are the
#  striatum; both dorsal striatum (caudate nucleus and putamen) and ventral
#  striatum (nucleus accumbens and olfactory tubercle), globus pallidus, ventral
#  pallidum, substantia nigra, and subthalamic nucleus.[3] Each of these components
#  has a complex internal anatomical and neurochemical organization. The largest
#  component, the striatum (dorsal and ventral), receives input from many brain
#  areas beyond the basal ganglia, but only sends output to other components of the
#  basal ganglia. The pallidum receives input from the striatum, and sends
#  inhibitory output to a number of motor-related areas. The substantia nigra is
#  the source of the striatal input of the neurotransmitter dopamine, which plays
#  an important role in basal ganglia function. The subthalamic nucleus receives
#  input mainly from the striatum and cerebral cortex, and projects to the globus
#  pallidus.
#
#
#   Popular theories implicate the basal ganglia primarily in action selection –
#   in helping to decide which of several possible behaviors to execute at any
#   given time. In more specific terms, the basal ganglia's primary function is
#   likely to control and regulate activities of the motor and premotor cortical
#   areas so that voluntary movements can be performed smoothly.[1][4] Experimental
#   studies show that the basal ganglia exert an inhibitory influence on a number
#   of motor systems, and that a release of this inhibition permits a motor system
#   to become active. The "behavior switching" that takes place within the basal
#   ganglia is influenced by signals from many parts of the brain, including the
#   prefrontal cortex, which plays a key role in executive functions.[2][5]
#
#   The importance of these subcortical nuclei for normal brain function and behavior
#   is emphasized by the numerous and diverse neurological conditions associated with
#   basal ganglia dysfunction, which include: disorders of behavior control such as
#   Tourette syndrome, hemiballismus, and obsessive–compulsive disorder; dystonia;
#   addiction; and movement disorders, the most notable of which are Parkinson's
#   disease, which involves degeneration of the dopamine-producing cells in the
#   substantia nigra pars compacta, and Huntington's disease, which primarily involves
#   damage to the striatum.[1][3] The basal ganglia have a limbic sector whose components
#   are assigned distinct names: the nucleus accumbens, ventral pallidum, and ventral
#   tegmental area (VTA). There is considerable evidence that this limbic part plays
#   a central role in reward learning, particularly a pathway (mesolimbic pathway)
#   from the VTA to the nucleus accumbens that uses the neurotransmitter dopamine.
#   A number of highly addictive drugs, including cocaine, amphetamine, and nicotine,
#   are thought to work by increasing the efficacy of this dopamine signal. There is
#   also evidence implicating overactivity of the VTA dopaminergic projection in schizophrenia.[6]
#
read.expression.data <- function()
{
  subjects <- c("H0351.1009", "H0351.1012", "H0351.1015", "H0351.1016", "H0351.2001", "H0351.2002")

  data <- list()
  for(subject in subjects){
     data.filename <- sprintf("data/%s.RData", subject)
     metadata.filename <- sprintf("data/%s.metadata.csv", subject)
     load(data.filename)
     tbl.md <- read.table(metadata.filename, sep=",", nrow=-1, as.is=TRUE, header=TRUE)
     printf("--- %s mtx: %d %d, metadata: %d %d", subject, nrow(expr), ncol(expr), nrow(tbl.md), ncol(tbl.md))
     sampleNames <- sprintf("%s.%s", tbl.md$structure_acronym, subject)
     colnames(expr) <- sampleNames
     data[[subject]] <- list(data=expr, md=tbl.md)
     print(fivenum(expr))
     } # for subject

   invisible(data)

} # read.expression.data
#----------------------------------------------------------------------------------------------------
create.htt.model <- function(shoulder)
{
   target.gene <- "HTT"
   gene.tss <- 3074678

   #gene.tss <- 3074510
   # 1.5kb region covering tss of HTT and HTT-AS: chr4:3,073,917-3,075,476
   roi <- "chr4:3,073,917-3,075,476"
   showGenomicRegion(tv, roi)
   chrom.loc <- parseChromLocString(getGenomicRegion(tv))
   #shoulder <- 2000
   #shoulder <- 0
   loc.chrom <- chrom.loc$chrom
   loc.start <- chrom.loc$start - shoulder
   loc.end   <- chrom.loc$end   + shoulder

   expanded.region <- with(chrom.loc, sprintf("%s:%d-%d", chrom, start-shoulder, end+shoulder))
   showGenomicRegion(tv, expanded.region)

   db.hint.20 <- sprintf("postgres://%s", "whovian/brain_hint_20")
   db.hint.16 <- sprintf("postgres://%s", "whovian/brain_hint_16")
   db.wellington.20 <- sprintf("postgres://%s", "whovian/brain_wellington_20")
   db.wellington.16 <- sprintf("postgres://%s", "whovian/brain_wellington_16")

   sources <- c(dhs="encodeHumanDHS", hint20=db.hint.20, hint16=db.hint.16,
                w20=db.wellington.20, w16=db.wellington.16)

   tbls.regions <- getRegulatoryChromosomalRegions(trena, loc.chrom, loc.start, loc.end, as.character(sources),
                                                   target.gene, gene.tss)
   colors <- rainbow_hcl(length(sources))

   # in exploratory work, you may want, now and again, to remova all tracks but the first, Gencode.v24
   #    removeTracksByName(tv, getTrackNames(tv)[-1])

   for(i in seq_len(length(sources))){
      sourceName <- names(sources)[i]
      tbls.regions[[i]]$source <- sourceName   # use this for loop to add source column to each table
      tbl.regions.trimmed <- tbls.regions[[i]][, c(1,2,3,4,6)]  # chrom, start, end, name, score
      addBedTrackFromDataFrame(tv, sourceName, tbl.regions.trimmed, color=colors[i])
      }

   tbl.motifs <- do.call(rbind, tbls.regions)
   shortMotifs <- unlist(lapply(strsplit(tbl.motifs$motifName, "-"), function(tokens) return(tokens[length(tokens)])))
   tbl.motifs$shortMotif <- shortMotifs
   tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="TFClass", expand.rows=TRUE)
   motifs.without.tf <- which(is.na(tbl.motifs$geneSymbol))
   if(length(motifs.without.tf) > 0)
      tbl.motifs <- tbl.motifs[-motifs.without.tf,]
   length(unique(tbl.motifs$geneSymbol))

   mtx.rna <- all.data[[1]]$data
   htt.variance <- var(mtx.rna["HTT",])
   variance <- apply(all.data[[1]]$data, 1, var)
   quiet.genes <- as.integer(which(variance < (htt.variance/1.2)))
   length(quiet.genes)
   if(length(quiet.genes) > 1)
      mtx.rna <- mtx.rna[-quiet.genes,]
   print(dim(mtx.rna))
   stopifnot("HTT" %in% rownames(mtx.rna))

   tbl.geneModel <- createGeneModel(trena, target.gene,
                                    c("lasso", "randomForest", "pearson", "spearman"),
                                    tbl.motifs, mtx.rna)

   colnames(tbl.geneModel) <- gsub(".", "", colnames(tbl.geneModel), fixed=TRUE)
   return(list(model=tbl.geneModel, regions=tbl.motifs))


} # create.hdd.model
#----------------------------------------------------------------------------------------------------
run <- function()
{
   x1 <- create.htt.model(shoulder=1000)
   tbl.model.strong.1 <- subset(x1$model, pcaMax >= 1)
   tbl.regions.strong.1 <- subset(x1$regions, geneSymbol %in% tbl.model.strong.1$gene)
   study.1 <- list(model=tbl.model.strong.1, regions=tbl.regions.strong.1)

   x2 <- create.htt.model(shoulder=5000)
   tbl.model.strong.2 <- subset(x2$model, pcaMax >= 1)
   tbl.regions.strong.2 <- subset(x2$regions, geneSymbol %in% tbl.model.strong.2$gene)
   study.2 <- list(model=tbl.model.strong.2, regions=tbl.regions.strong.2)

   models <- list("1kb"=study.1, "5kb"=study.2)

   g <- buildMultiModelGraph(tv, targetGene="HTT", models)
   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   setGraph(tv, g.lo, names(models))
   setStyle(tv, "style.js")


#   tbl.geneModel.small <- tbl.geneModel[1:3,]
#
#   tbl.motifs.filtered <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel$gene)
#   tbl.motifs.filtered.small <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel.small$gene)
#
#   model.1 <- list(model=tbl.geneModel.small, regions=tbl.motifs.filtered.small)
#   model.2 <- list(model=tbl.geneModel, regions=tbl.motifs.filtered)
#
#   models <- list(subject.1009=model.1, big=model.2)
#
#   g <- buildMultiModelGraph(tv, targetGene="HTT", models)
#   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
#   setGraph(tv, g.lo, names(models))
#   setStyle(tv, "style.js")
#   fit(tv)


} # run
#----------------------------------------------------------------------------------------------------
if(!exists("all.data")){
  printf("loading all allen brain atlas data into 'all.data', a list of six matrix pairs, expression & metadata")
  all.data <- read.expression.data()
  }
