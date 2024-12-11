torque.template <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/methods/kallistoQuantification.sh"

escape.forwardslash <- function(fp) {
    gsub("/", "\\/", fp, fixed = TRUE)
}

generateTorqueQuantificationScripts <- function(out.dir, raw.read.dirs, 
    torque.jobs.dir, reference.transcriptome.index) {
    for (reads.dir in raw.read.dirs) {
        sample.no <- sub("^.*/", "", reads.dir)
        quants.dir <- file.path(out.dir, sample.no)
        if (!dir.exists(quants.dir)) {
            dir.create(quants.dir, recursive = TRUE)
        }
        trimmed.reads <- system(paste0("find ", reads.dir, " -type f -name '*T.fq.gz'"), 
            intern = TRUE)
        if (!dir.exists(torque.jobs.dir)) {
            dir.create(torque.jobs.dir, recursive = TRUE)
        }
        torque.file <- file.path(torque.jobs.dir, paste0("quantify_", sample.no, 
            ".sh"))
        system(paste0("sed -e 's/OUT_DIR/", escape.forwardslash(quants.dir), 
            "/' -e 's/REFERENCE_TRANSCRIPTOME_INDEX/", escape.forwardslash(reference.transcriptome.index), 
            "/' -e 's/FRWRD/", escape.forwardslash(trimmed.reads[[1]]), 
            "/' -e 's/BCKWRD/", escape.forwardslash(trimmed.reads[[2]]), 
            "/' ", torque.template, " > ", torque.file))
    }
}


#' Leaves PCP1
out.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/results/quantifications/LeavesPCP1"
raw.read.dirs <- system("find /mnt/data/usadel/VANESSA/LeavesPCP1/raw_data/F19FTSEUHT0641_POTtzjE/results -type d -regextype posix-awk -regex '.*/[0-9]+'", 
    intern = TRUE)
torque.jobs.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/methods/torque_scripts/LeavesPCP1"
reference.transcriptome.index <- "/mnt/data/asis/Potato-Genome/RNA-Seq-Projects/PCP1_Leaves/DM_1-3_516_R44_potato.v6.1.transcriptome_curated.fa.index"
generateTorqueQuantificationScripts(out.dir, raw.read.dirs, torque.jobs.dir, reference.transcriptome.index)


#' Leaves T6P
out.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/results/quantifications/LeavesT6P"
raw.read.dirs <- system("find /mnt/data/usadel/VANESSA/LeavesT6P/raw_data/F19FTSEUHT1836_POTwpeE/results -type d -regextype posix-awk -regex '.*/Sample[0-9]+'", 
    intern = TRUE)
torque.jobs.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/methods/torque_scripts/LeavesT6P"
reference.transcriptome.index <- "/mnt/data/asis/Potato-Genome/RNA-Seq-Projects/T6P_Leaves/DM_1-3_516_R44_potato.v6.1.transcriptome_curated.fa.index"
generateTorqueQuantificationScripts(out.dir, raw.read.dirs, torque.jobs.dir, reference.transcriptome.index)


#' Apex
out.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/results/quantifications/Apex"
raw.read.dirs <- system("find /mnt/data/usadel/VANESSA/apex/F20FTSEUHT0704_POTojcE/results -type d -regextype posix-awk -regex '.*/[0-9]+'", 
    intern = TRUE)
torque.jobs.dir <- "/mnt/data/asis/Potato-PCP-T6P-Vanessa-Wahl/methods/torque_scripts/Apex"
reference.transcriptome.index <- "/mnt/data/asis/Potato-Genome/RNA-Seq-Projects/Apex/DM_1-3_516_R44_potato.v6.1.transcriptome_curated.fa.index"
generateTorqueQuantificationScripts(out.dir, raw.read.dirs, torque.jobs.dir, reference.transcriptome.index)
