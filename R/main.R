#' A sequence-based method to evaluate the effectiveness of siRNA in human and other models
#'
#' A sequence-based method to evaluate the effectiveness of siRNA in human and other models
#'
#' @import parallel
#' @import data.table
#' @import Rsamtools
#' @import GenomicRanges
#' @import rtracklayer
#' @import Rbowtie
#' @import GenomicAlignments
#' @import stringr
#' @import dplyr
#'
#' @param sirna A data.frame object, with at least two columns of "id","antisense"
#' @param bam.files One or multiple bam file. sorted and indexed bam are recommended
#' @param gtf_file A gtf file for gene and transcript annotation. "gtf_file" is not need if  "genomic_location" is not null
#' @param gene Symbol of interest gene, e.g. E2F1. It is needed if "gtf_file" is not null
#' @param genomic_location A string with format like "chr1:1000:2000". "genomic_location" is not needed if "gtf_file" and "gene" is not null
#' @param tmp_dir A directory store the temporary outputs
#' @param remove_header A numeric value. For siRNA the 1st 5'nucleic acid is usually not needed during alignment
#' @param remove_overhang A numeric value. siRNA usually have a 2-nt overhang and they will be removed during alignment
#' @param paired Loigic value. The bam file is derived from paired-end RNA-seq
#' @param num_threads The number of threads
#'
#' @author Guofeng Meng
#'
#'
#' @details This function can map the input siRNA sequence to RNA-seq reads and evaluate if (1) siRNA match
#'
#' @return data.frame
#'
#' @examples
#' #siRNASeqEva(sirna, bam.files, gtf_file="/path/gtf_file.gtf", gene="E2F1")
#'
#' @export




siRNASeqEva<-function(sirna, bam.files=NULL, gtf_file=NULL, gene_symbol=NULL, genomic_location=NULL, tmp_dir="temp", remove_header=1, remove_overhang=2, paired=TRUE, num_threads=1, add_chr=FALSE){
    gr_target=NULL
    final_gtf=NULL
    if(!is.null(gtf_file) ){
        if(!is.null(gene_symbol)){
             gtf <- as.data.table(rtracklayer::import(gtf_file))
             if(!"gene_name" %in% names(gtf) & "gene" %in% names(gtf))
               gtf$gene_name=as.vector(gtf$gene)
             gtf = gtf[ gene_name ==gene_symbol & (type=="exon" | type=="five_prime_utr"| type=="three_prime_utr")]
             final_gtf=gtf
             chrs=unique(gtf$seqnames)
             if(add_chr)
				        chrs=paste0("chr", chrs)
            gr_target=makeGRangesFromDataFrame(data.frame(chr=chrs, start=min(gtf$start-1000), end=max(gtf$end+1000),strand="+"))
        }
    }else{
        if(!is.null(genomic_location)){
            genomic_location=sub("-",":", genomic_location)
            locs=strsplit(genomic_location,":")[[1]]
            if(length(locs)>=3){
                gr_target=makeGRangesFromDataFrame(data.frame(chr=locs[1], start=max(locs[2]-1000,0), end=locs[3] +1000,strand="+"))
            }else{
                gr_target=makeGRangesFromDataFrame(data.frame(chr=locs[1], start=1, end=100000000,strand="+"))
            }
        }else{
            print("Error: No target information is provided!")
            return(NULL)
        }
    }
    if(is.null(bam.files)){
        print("Error: no input bams!")
        return(NULL)
    }
    which <- gr_target
    what <- c("qname", "rname","strand", "pos", "qwidth", "cigar","seq")
    
    tmp_dir=sub("/$","",tmp_dir)
    if(!dir.exists(tmp_dir))
        dir.create(tmp_dir)
    used_bam_file=NULL
    if(length(bam.files)==1){
        param <- ScanBamParam(which=which, what=what, flag=scanBamFlag(isUnmappedQuery=FALSE))
        if(!file.exists(paste0(bam.files,".bai"))){
                out_ff=paste0(sub(".bam$","",bam.files), "_sorted")
                out_ff2=paste0(sub(".bam$","",bam.files), "_sorted.bam")
                sortBam(bam.files, out_ff, nThreads=num_threads)
                indexBam(out_ff2)
                in_bam=scanBam(out_ff2, param=param)[[1]]
                used_bam_file=out_ff2
                if(paired){
                  gal<-readGAlignmentPairs(out_ff2, param=param)
                }else{
                  gal<-readGAlignments(out_ff2, param=param)
                }

            }else{
                in_bam=scanBam(bam.files, param=param)[[1]]
                used_bam_file=files
                if(paired){
                  gal<-readGAlignmentPairs(bam.files, param=param)
                }else{
                  gal<-readGAlignments(bam.files, param=param)
                }
            }

    }else{
        new.bam.files=sapply(bam.files, function(ff){
             if(!file.exists(paste0(ff,".bai"))){
                out_ff=paste0(sub(".bam$","",ff), "_sorted")
                out_ff2=paste0(sub(".bam$","",ff), "_sorted.bam")
                sortBam(ff, out_ff, nThreads=num_threads)
                indexBam(out_ff2)
                return(out_ff2)
             }else{
                return(ff)
             }
        })
        param <- ScanBamParam( what=what, flag=scanBamFlag(isUnmappedQuery=FALSE))
        out_file=mergeBam(new.bam.files, destination=paste0(tmp_dir,"/merged_bam.bam"), region=gr_target, indexDestination=TRUE, overwrite=TRUE)
        in_bam=scanBam(out_file, param=param)[[1]]
        used_bam_file=out_file
        if(paired){
          gal<-readGAlignmentPairs(out_file,param=param)
        }else{
          gal<-readGAlignments(out_file,param=param)
        }
    }

    gv=GenomicAlignments::coverage(gal)[[ as.vector(seqnames(gr_target))]]
    bam_df=data.frame(chr=in_bam[["rname"]], star=in_bam[["pos"]],width=in_bam[["qwidth"]],strand=in_bam[["strand"]],seqs=in_bam[["seq"]], id=in_bam[["qname"]], cigar=in_bam[["cigar"]])
    
    if(nrow(bam_df)==0)
        return(NULL)
    bam_df$seq_name=paste0("R", 1:nrow(bam_df))
    read_pos = as.vector(bam_df$star)
    read_cigar = as.vector(bam_df$cigar)
    names(read_pos)<-as.vector(bam_df$seq_name)
    names(read_cigar)<-as.vector(bam_df$seq_name)
    read_cigar=str_replace_all(read_cigar,"([MSNIDHX=])","#\\1#")
    read_cigar=strsplit(read_cigar,"#")
    read_cigar=as.data.frame(t(sapply(read_cigar, function(x) append(x, rep(-1, 30-length(x))))))
    read_cigar$star=as.vector(bam_df$star)
    row.names(read_cigar)<-as.vector(bam_df$seq_name)

    write.table(paste(paste0(">",as.vector(bam_df$seq_name)), as.vector(bam_df$seqs), sep="\n"), paste0(tmp_dir,"/reads_fasta.fa"),row.names=F,col.name=F,quote=F)
    reference_build_report2 <- bowtie_build(references=paste0(tmp_dir,"/reads_fasta.fa"), outdir=tmp_dir, prefix="reads", force=TRUE, threads=num_threads)
    if(is.data.frame(sirna)){
        write.table(paste(paste0(">",as.vector(sirna$id)), gsub("U","T",toupper(as.vector(sirna$antisense))), sep="\n"), paste0(tmp_dir,"/sirna_fasta.fa"),row.names=F,col.name=F,quote=F)
    }else{
        print("error")
    }
    bowtie(sequences=paste0(tmp_dir,"/sirna_fasta.fa"), index=paste0(tmp_dir,"/reads"), type="single", outfile=paste0(tmp_dir,"/sirna.sam"), all=TRUE, f=TRUE, v=3, trim5=remove_header, trim3=remove_overhang, threads=num_threads, force=TRUE)
    res=fread(paste0(tmp_dir,"/sirna.sam"))[,c("V1","V2","V3","V4","V8")]
    if(nrow(res)==0)
        return(NULL)
    spos = as.vector(res$V4)
    vpos=rep(-1, length(spos))
    used.cigar=read_cigar[as.vector(res$V3),]
    done=rep(FALSE, length(spos))
    done[is.na(as.vector(used.cigar$star))]=TRUE

    while(!all(done)){
        #print("ok")
        names(used.cigar)<-c(paste0("S",1:(ncol(used.cigar)-1)),"star")
        done[used.cigar$S1=="-1"]=TRUE
        wh0= as.numeric(as.vector(used.cigar$S1)) >= spos
        wh1= !done & as.vector(used.cigar$S2)== "M"
        wh2= !done & as.vector(used.cigar$S2)== "N"
        wh3= !done & as.vector(used.cigar$S2)== "S"
        wh4= !done & as.vector(used.cigar$S2)== "I"
        wh5= !done & as.vector(used.cigar$S2)== "D"
        if(length(which(wh0 & wh1))!=0 ){
            vpos[wh0 & wh1]=as.vector(used.cigar$star)[wh0 & wh1] + spos[wh0 & wh1]
            done[wh0 & wh1]=TRUE
        }
        if(length(which(!wh0 & wh1))!=0 ){
            used.cigar$star[!wh0 & wh1] = used.cigar$star[!wh0 & wh1] + as.numeric(as.vector(used.cigar$S1[!wh0 & wh1]))
            spos[!wh0 & wh1]=spos[!wh0 & wh1] - as.numeric(as.vector(used.cigar$S1[!wh0 & wh1]))
        }
        if(length(which(wh0 & wh3))!=0 ){
            done[wh0 & wh3]=TRUE
        }
        if(length(which(!wh0 & wh3))!=0 ){
            spos[!wh0 & wh3]=spos[!wh0 & wh3] - as.numeric(as.vector(used.cigar$S1[!wh0 & wh3]))
        }
        if(length(which(wh2| wh5 ))!=0 ){
            used.cigar$star[wh2| wh5] = used.cigar$star[wh2| wh5] + as.numeric(as.vector(used.cigar$S1))[wh2| wh5]
        }
        if(length(which(wh4))!=0 ){
            used.cigar$star[wh4] = used.cigar$star[wh4] - as.numeric(as.vector(used.cigar$S1))[wh4]
        }
        used.cigar=used.cigar[,c(-1,-2)]
        if(ncol(used.cigar)==0)
            break()
    }
    res$map_loc=vpos
    rm(wh0, wh1, wh2, wh3, wh4, wh5, spos, vpos)

    mismatch2=unlist(mclapply(as.vector(res$V8), function(x) {
      if(x=="")
          return(0)
      return(length(strsplit(x,",")[[1]]))
    }, mc.cores=num_threads))

    res$mismatch=mismatch2
    
    mapping_info=mclapply(as.vector(sirna$id), function(sr){
        #print(sr)
        sub.res=subset(res, V1==sr)
        locs=unique(as.vector(sub.res$map_loc) )
        locs=locs[locs!=-1]
        if(length(locs)==0)
            return(NULL)
        cv=as.vector(gv[locs])
        if(length(locs) >1){
            out=t(sapply(locs, function(lc){
                table(subset(sub.res, map_loc==lc)$mismatch)[c("0","1","2","3")]
            }))
            out[is.na(out)]=0
            colnames(out)<-c("Mismatch0","Mismatch1","Mismatch2","Mismatch3")
            out[,2]=out[,2] +out[,1]
            out[,3]=out[,3] +out[,2]
            out[,4]=out[,4] +out[,3]
            cbind(data.frame(id=sr, chr=as.vector(seqnames(gr_target)) , map_loc=locs, coverage=cv), out)
        }else{
            out=as.vector(table(subset(sub.res, map_loc==locs)$mismatch)[c("0","1","2","3")])
            out[is.na(out)]=0
            out[2]=out[2] +out[1]
            out[3]=out[3] +out[2]
            out[4]=out[4] +out[3]
            data.frame(id=sr, chr=as.vector(seqnames(gr_target)) , map_loc=locs, coverage=cv, Mismatch0=out[1],Mismatch1=out[2],Mismatch2=out[3],Mismatch3=out[4])
        }
    }, mc.cores=num_threads)
    mapping_info=mapping_info[!sapply(mapping_info, is.null)]
    if(length(mapping_info)==0)
      return(NULL)
    mapping_info=bind_rows(mapping_info)
    mapping_info=mapping_info[as.vector(mapping_info$Mismatch3)/(1+ as.vector(mapping_info$coverage)) > 0.1 & as.vector(mapping_info$coverage) > 3,]
    outcome=list(mapping_info=mapping_info, sirna=sirna, gtf=final_gtf, ranges=gr_target, bam = used_bam_file)
    return(outcome)
}



filter_redundency <-function(obj, ordered=TRUE){
	allids=unique(as.vector(obj$sirna$id))
	xx=dplyr::bind_rows(lapply(allids, function(ii){
		sub.x=subset(obj$mapping_info, id==ii)
		if(nrow(sub.x)==0){
			return(data.frame(id=ii, chr=NA, map_loc=NA, coverage=0, Mismatch0=0, Mismatch1=0, Mismatch2=0, Mismatch3=0, mark="NF"))
		}
		match0=as.vector(sub.x$Mismatch0)
		if(any(match0 > 3)){
			sub.x=subset(sub.x, Mismatch0 > 3)
		}
		sub.x$mark="UN"
		if(nrow(sub.x)==1)
			return(sub.x)
		match1=as.vector(sub.x$Mismatch1)
		match0=as.vector(sub.x$Mismatch0)
		max.match1=max(match1)
		wh=which(match0/match1 > 0.8  &  match1 > max.match1 * 0.6)
		if(length(wh[wh])==1){
			return(sub.x[wh,])
		}
		if(length(wh[wh])>1){
			sub.x=sub.x[wh,]
		}
		has.pos=as.vector(sub.x$map_loc)
		wh=which(allids==ii)
		left.tag1=wh-30
		right.tag1=wh+10
		left.tag2=wh-10
		right.tag2=wh+30
		if(right.tag1 > length(allids))
			right.tag1 = length(allids)
		if(right.tag2 > length(allids))
			right.tag2 = length(allids)
		if(left.tag1 <0)
			left.tag1=1
		if(left.tag2 <0)
			left.tag2=1
		  
		rg1=left.tag1:right.tag1
		rg2=left.tag1:right.tag2 
		
		nearby.id1=allids[rg1]
		nearby.id2=allids[rg2]
		
		nearby.pos1=as.vector(subset(obj$mapping_info, id%in%nearby.id1)$map_loc)
		nearby.pos2=as.vector(subset(obj$mapping_info, id%in%nearby.id2)$map_loc)
		sel=which.min(abs(sapply(has.pos, function(pos) min(mean(nearby.pos1-pos), mean(nearby.pos2-pos)))))
		sel.x=sub.x[sel,]
		if(as.vector(sel.x$Mismatch0)/as.vector(sel.x$Mismatch3) < 0.6){
			sub.x$mark="DU"
			return(sub.x)
		}
		return(sel.x)
	}))
	
	obj$refined_mapping_info=xx
	return(obj)
}

plot_siRNA<-function(obj, transcript=c("ENST00000349496","ENST00000396185"), chromosome=NULL, from=NULL, to=NULL, ratio_mismatch=0.6){
	gtf=obj$gtf
	sub.gtf=gtf[gtf$transcript_id %in% transcript,]
	
	if(is.null(chromosome))
		chromosome=unique(as.vector(sub.gtf$seqnames))
	if(is.null(from))
		from=min(as.vector(sub.gtf$start))
	if(is.null(to))
		to=max(as.vector(sub.gtf$end))
	# mapping track
	sub.zz2=subset(obj$refined_mapping_info, !is.na(map_loc))
	sub.zz2$start=as.vector(sub.zz2$map_loc)
	sub.zz2$end=as.vector(sub.zz2$map_loc) +0
	tps=apply(sub.zz2,1, function(x) {
		if(x["Mismatch3"]==0)
			return("No match")
		if(as.numeric(x["Mismatch2"])/as.numeric(x["Mismatch3"]) <=ratio_mismatch)
			return("Mismatch3")
		if(as.numeric(x["Mismatch1"])/as.numeric(x["Mismatch3"]) <=ratio_mismatch)
			return("Mismatch2")
		if(as.numeric(x["Mismatch0"])/as.numeric(x["Mismatch3"]) <=ratio_mismatch)
			return("Mismatch1")
		return("Full")
	})
	sub.zz2$tps=tps
	input=makeGRangesFromDataFrame(sub.zz2, keep.extra.columns=TRUE)
	options(ucscChromosomeNames = FALSE)
	aTrack <- AnnotationTrack(input, group=tps, feature=tps, stacking="full", name="siRNA",groupAnnotation="group",just.group = "above")
  
	# gtf track
	
	gtf.gr=makeGRangesFromDataFrame(sub.gtf, keep.extra.columns=TRUE)
	GeneTrack=GeneRegionTrack(range=gtf.gr,  stacking="full", name="Gene", transcript=as.vector(gtf.gr$transcript_id), gene=as.vector(gtf.gr$transcript_id),showId=TRUE )
	
	# bam coverage track
	coverageTrack <- DataTrack(range=obj$bam,  name="Coverage", chromosome=chromosome, type = "histogram")
	
	#snp track
	snp=fread("/data/meng/database/snp/snp.txt",skip=1, header=F)
	names(snp)<-c("chr","start","end","snp")
	sub.snp=unique(subset(snp, chr==chromosome | chr==paste0("chr",chromosome) & start >=from & end <= to,1:3))
	if(nrow(sub.snp)!=0){
		sub.snp$chr=chromosome
		sub.snp$end=as.vector(sub.snp$start)
		gr.snp=makeGRangesFromDataFrame(sub.snp, keep.extra.columns=TRUE) 
		gr.snp=gr.snp[unique(queryHits(findOverlaps(gr.snp, gtf.gr)))]
		 
		snpTrack <- AnnotationTrack(gr.snp,  name="SNP")
		displayPars(snpTrack) <- list(fill="red",col="red")
	}else{
		gr.snp=makeGRangesFromDataFrame(sub.snp, keep.extra.columns=TRUE) 
		snpTrack <- AnnotationTrack(gr.snp,  name="SNP")
	}   
  
	htTrack <- HighlightTrack(trackList = list(aTrack,snpTrack), start = start(gr.snp), width = 1, chromosome = chromosome, fill="yellow")
	displayPars(htTrack) <- list(fill="lightgreen",col="lightgreen")                     
	plotTracks(list(htTrack,coverageTrack, GeneTrack), chromosome=chromosome, from=from, to=to, Mismatch1 = "darkred", Mismatch0 = "darkgreen", Mismatch2 = "brown", Full="cornflowerblue", just.group = "above",sizes=c(3,0.3,0.5,0.5),transcriptAnnotation = "transcript")
  
}
