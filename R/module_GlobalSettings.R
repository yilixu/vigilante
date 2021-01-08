#' Prepare data for vigilante analyses
#'
#' (Internal) In order to properly perform downstream analyses, vigilante needs to prepare all input data files in the working directory before continuing. This internal process includes several steps such as renaming data files, setting group-related parameters, setting ENSEMBL refernece, etc.
#'
#' @keywords internal

# v_prepareVdata function
v_prepareVdata = function(studyID, studyID_regex, studyID_altered, speciesID, fileNum, fileDNA_mafNum, clinicalFeature) {

  # set file name matching patterns
  name_pattern = paste0(studyID_regex, "_[[:digit:]]+_[[:digit:]]+")

  # capture original file name
  file_name = list()
  for (i in 1:length(name_pattern)) {
    file_name[[i]] = dir(path = ".", pattern = name_pattern[i])
  }
  rm(i)
  file_name = unlist(file_name)
  realID = gsub(glue::glue("^({name_pattern}).+$"), "\\1", file_name)
  fileID = realID

  # extract assayID, set realID, aliasID and sampleNum, generate default groupInfo
  if (!("groupInfo.csv" %in% dir())) {

    # groupInfo.csv file generation notice
    print(as.character(glue::glue("groupInfo.csv file not detected, a default one will be generated in the working directory {getwd()}")))

    # reform original file name, extract assayID and set realID and aliasID
    realID = unique(realID)
    assayID = sprintf("A%05d", 1:length(realID))
    realID = paste0(assayID, "_", realID)
    realID = sort(realID)
    realID = gsub("A[[:digit:]]{5}_", "", realID)
    aliasID = gsub(paste0(studyID_regex, "_"), "", realID)
    assayID = sort(assayID)

    # if condition for whether swap realID and aliasID
    if (studyID_altered == TRUE) {
      aliasID = realID
      realID = gsub("A0", paste0(studyID, "_"), assayID)
    }

    # set default groupInfo and output it to csv file
    if (clinicalFeature == TRUE) {
      groupInfo = data.frame(assayID, Group = "Undefined", MAF_Group = "All", realID, aliasID, CliFeaC_ = "Undefined", CliFeaN_ = "Undefined", stringsAsFactors = FALSE)
    } else {
      groupInfo = data.frame(assayID, Group = "Undefined", MAF_Group = "All", realID, aliasID, stringsAsFactors = FALSE)
    }
    write.csv(groupInfo, file = "./groupInfo.csv", row.names = FALSE, quote = FALSE)
    print("groupInfo.csv generated in the working directory, please follow the instruction to fill in the groupInfo.csv before continue to the next step")

    # remove default extracted aliasID
    rm(aliasID)

    # if groupInfo.csv already exists
  } else {
    groupInfo = read.csv("groupInfo.csv", stringsAsFactors = FALSE)
    realID = groupInfo["realID"]
    realID = unlist(realID, use.names = FALSE)
    assayID = groupInfo["assayID"]
    assayID = unlist(assayID, use.names = FALSE)
  }

  # set sample number
  sampleNum = length(assayID) # number of samples

  # remove default extracted groupInfo
  rm(groupInfo)

  # ask user to make sure groupInfo.csv is properly filled in
  status_groupInfo = menu(choices = c("Yes", "No"), title = "\nIs groupInfo.csv properly filled in?")
  if (status_groupInfo == 1) {
    print("groupInfo.csv is properly filled in, now continue to the next step")
  } else if (status_groupInfo == 2) {
    print("Please follow the instruction to fill in the groupInfo.csv before continue to the next step")
    stop()
  } else {
    print("Please choose a valid answer")
    stop()
  }
  rm(status_groupInfo)

  # set group-related parameters and plot_name
  print("Setting group-related parameters")

  # load group info
  groupInfo = read.csv("groupInfo.csv", stringsAsFactors = FALSE)

  # include CliFea columns based on clinicalFeature condition
  if (clinicalFeature == TRUE) {
    groupInfo = groupInfo
  } else {
    groupInfo = groupInfo[, 1:5]
  }

  # set group name
  grpName = unique(groupInfo[, 2])
  grpName_maf = unique(groupInfo[, 3])

  # group up
  grp.list = list()
  for (i in 1:length(grpName)) {
    grp.list[[i]] = subset(groupInfo, subset = Group == grpName[i])[, 1]
  }
  rm(i)
  names(grp.list) = grpName

  # group up maf data
  grp.list.maf = list()
  for (i in 1:length(grpName_maf)) {
    grp.list.maf[[i]] = subset(groupInfo, subset = MAF_Group == grpName_maf[i])[, 1]
  }
  rm(i)
  names(grp.list.maf) = grpName_maf

  # set plotting number to reflect actual plotting order
  aliasID = unlist(groupInfo$aliasID)
  plot_name = c(aliasID, paste0("Group_", 1:length(grpName), "_", grpName))
  rm(aliasID)

  # check if input data files are already in position (in separate folders)
  temp_folderName = dir(path = ".", pattern = paste0(studyID, "_A[[:digit:]]{5}dir"))
  if (length(temp_folderName) == 0) {

    # ask user to make sure original data files are backed up before renaming
    dataBackUp_status = menu(choices = c("YES, they are all backed up and ready to go through renaming (Continue now)", "Not yet, but I will do it (Come back later)", "No, I don't want them to be renamed (Stop and quit)"), title = "\nIn order to properly perform downstream analyses, vigilante needs to rename all input data files in the working directory before continuing (only file names will be modified). Are original input data files safely backed up in places other than the working directory?")
    if (dataBackUp_status == 1) {
      print("Original input data files are safely backed up, now continue to the next step")
    } else if (dataBackUp_status == 2) {
      print("It is recommended to back up original input data files before continue to the next step")
      stop()
    } else if (dataBackUp_status == 3) {
      print("Thank you for using vigilante")
      stop()
    } else {
      print("Please choose a valid answer")
      stop()
    }
    rm(dataBackUp_status)

    # rename files
    print("Start renaming input data files")

    # capture old file name and set corresponding new file name
    file_name_old = file_name
    file_name_new = data.frame(cbind(file_name, fileID), stringsAsFactors = FALSE)
    if (studyID_altered == TRUE) {
      aliasID = groupInfo["aliasID"]
      aliasID = unlist(aliasID, use.names = FALSE)
      temp_pair = data.frame(cbind(aliasID, assayID), stringsAsFactors = FALSE)
      rm(aliasID)
    } else {
      temp_pair = data.frame(cbind(realID, assayID), stringsAsFactors = FALSE)
    }
    colnames(temp_pair) = c("fileID", "assayID")
    file_name_new = plyr::join(file_name_new, temp_pair, by = "fileID", type = "left")
    for (i in 1:nrow(file_name_new)) {
      file_name_new[i, 4] = gsub(name_pattern, paste0(studyID, "_", file_name_new[i, 3]), file_name_new[i, 1])
    }
    rm(i)
    file_name_new = file_name_new[, 4]

    # rename files
    for (i in 1:length(file_name_old)) {
      file.rename(file_name_old[i], file_name_new[i])
    }
    rm(i)

    # remove temp name data
    rm(name_pattern, file_name, temp_pair, fileID, file_name_old, file_name_new)

    # create dummy files for incomplete samples, vks pkg required
    print("Checking if all the samples have the same number of input data files, if not, corresponding dummmy input data files (empty-content placeholders) will be generated to cover the missing parts")

    # capture existing file names
    file_baf = dir(path = ".", pattern = "baf.tsv")
    file_cna = dir(path = ".", pattern = "cna.tsv")
    file_maf = dir(path = ".", pattern = ".maf")
    file_ge = dir(path = ".", pattern = "cDNA_genes.sf")
    file_fuca = dir(path = ".", pattern = "-fusion-genes.txt")
    file_star = dir(path = ".", pattern = "star-fusion.*abridged")
    file_te = dir(path = ".", pattern = "cufflinks.isoforms.fpkm_tracking")

    # check and create dummy files, baf (1/7)
    if (length(file_baf) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_baf)) == 0) {
          file.copy(system.file("extdata", "template_dummy_exo.baf.tsv", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy_exo.baf.tsv"))
        }
      }
      rm(i)
    }
    rm(file_baf)

    # check and create dummy files, cna (2/7)
    if (length(file_cna) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_cna)) == 0) {
          file.copy(system.file("extdata", "template_dummy_exo.cna.tsv", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy_exo.cna.tsv"))
        }
      }
      rm(i)
    }
    rm(file_cna)

    # check and create dummy files, maf (3/7)
    if ((fileDNA_mafNum == 1 & length(file_maf) > 1)) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_maf)) == 0) {
          file.copy(system.file("extdata", "template_dummy.MuTect.seurat.strelka.maf", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.MuTect.seurat.strelka.maf"))
        }
      }
      rm(i)
    } else if (fileDNA_mafNum == 2 & length(file_maf) > 2) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_maf)) == 0) {
          file.copy(system.file("extdata", "template_dummy.MuTect.seurat.strelka.snvs.maf", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.MuTect.seurat.strelka.snvs.maf"))
          file.copy(system.file("extdata", "template_dummy.MuTect.seurat.strelka.indels.maf", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.MuTect.seurat.strelka.indels.maf"))
        } else if (length(grep(assayID[i], file_maf)) == 1) {
          file.copy(system.file("extdata", "template_dummy.MuTect.seurat.strelka.indels.maf", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.MuTect.seurat.strelka.indels.maf"))
        }
      }
      rm(i)
    }
    rm(file_maf)

    # check and create dummy files, ge (4/7)
    if (length(file_ge) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_ge)) == 0) {
          file.copy(system.file("extdata", "template_dummy_salmon_bc_cDNA_genes.sf", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy_salmon_bc_cDNA_genes.sf"))
        }
      }
      rm(i)
    }
    rm(file_ge)

    # check and create dummy files, fuca (5/7)
    if (length(file_fuca) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_fuca)) == 0) {
          file.copy(system.file("extdata", "template_dummy.final-list_candidate-fusion-genes.txt", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.final-list_candidate-fusion-genes.txt"))
        }
      }
      rm(i)
    }
    rm(file_fuca)

    # check and create dummy files, star (6/7)
    if (length(file_star) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_star)) == 0) {
          file.copy(system.file("extdata", "template_dummy.star-fusion.fusion_candidates.final.abridged.txt", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.star-fusion.fusion_candidates.final.abridged.txt"))
        }
      }
      rm(i)
    }
    rm(file_star)

    # check and create dummy files, te (7/7)
    if (length(file_te) > 1) {
      for (i in 1:length(assayID)) {
        if (length(grep(assayID[i], file_te)) == 0) {
          file.copy(system.file("extdata", "template_dummy.cufflinks.isoforms.fpkm_tracking", package = "vigilante.knights.sword", mustWork = TRUE), paste0(studyID, "_", assayID[i], "_dummy.cufflinks.isoforms.fpkm_tracking"))
        }
      }
      rm(i)
    }
    rm(file_te)

    # create separate folders and move files
    print("Checking completed, now moving all input data files into position")

    # capture old file path and set corresponding new file path
    file_path_old = dir(path = ".", pattern = paste0(studyID, "_A[[:digit:]]{5}[[:punct:]]{1}"))
    folder_name = paste0(studyID, "_", assayID, "dir")
    file_path_new = paste0("./", rep(folder_name, each = fileNum), "/")

    # create separate folders
    for (i in assayID) {
      dir.create(paste0(studyID, "_", i, "dir"))
    }
    rm(i)

    # move files
    for (i in 1:length(file_path_old)) {
      file.rename(file_path_old[i], paste0(file_path_new[i], file_path_old[i]))
    }
    rm(i)

    # remove temp path data
    rm(file_path_old, file_path_new)
  } else if (length(temp_folderName) > 0) {

    # skip renaming, dummy files, creating folders and moving files (first-time run preparation steps)
    print("Input data files are already in position, skipping first-time run preparation steps")

    # extract and reset folder_name
    folder_name = dir(path = ".", pattern = paste0(studyID, "_A[[:digit:]]{5}dir"))
  }
  rm(temp_folderName)

  # set ENSEMBL reference based on speciesID, vks pkg required
  print("Setting ENSEMBL reference for downstream analyses")
  GRCh38_ref = NULL
  GRCh37_ref = NULL
  GRCm38_ref = NULL
  if (speciesID == "hg38") {
    GRCh38_ref = vigilante.knights.sword::GRCh38_ref
  } else if (speciesID == "hg19") {
    GRCh37_ref = vigilante.knights.sword::GRCh37_ref
  } else if (speciesID == "mm10") {
    GRCm38_ref = vigilante.knights.sword::GRCm38_ref
  }

  # add variables to temp_globalSettings_returnList
  temp_prepareVdata_returnList = list(
    "groupInfo" = groupInfo,
    "grp.list" = grp.list,
    "grp.list.maf" = grp.list.maf,
    "assayID" = assayID,
    "folder_name" = folder_name,
    "grpName" = grpName,
    "grpName_maf" = grpName_maf,
    "plot_name" = plot_name,
    "realID" = realID,
    "sampleNum" = sampleNum
  )

  # remove variables already added to temp_prepareVdata_returnList
  rm(list = names(temp_prepareVdata_returnList))

  # add ENSEMBL reference to temp_prepareVdata_returnList
  temp_prepareVdata_returnList = c(temp_prepareVdata_returnList, GRCh38_ref, GRCh37_ref, GRCm38_ref)

  # remove ENSEMBL reference already added to temp_prepareVdata_returnList
  if (speciesID == "hg38") {
    rm(GRCh38_ref)
  } else if (speciesID == "hg19") {
    rm(GRCh37_ref)
  } else if (speciesID == "mm10") {
    rm(GRCm38_ref)
  }

  # end of v_prepareVdata function
  print("Data preparation for vigilante completed")
  return(temp_prepareVdata_returnList)
}



#' Global settings for vigilante
#'
#' vigilante is designed to work in a highly self-contained pipeline style. To ensure that, there are several global settings need to be specified by the user. Also, in order to properly perform downstream analyses, vigilante needs to prepare all input data files in the working directory before continuing. This internal process includes several steps such as renaming data files, setting group-related parameters, setting ENSEMBL refernece, etc.
#'
#' @param studyID character, study or project name, will be used in multiple naming situations (such as on the plot, or in the output file names), it is recommended to be concise and meaningful (e.g. "ProstateCancer", or "PCa" for short; "TripleNegativeBreastCancer", or "TNBC" for short), also avoid using special characters.
#' @param studyID_regex regular expression, type "?regex" in the console to see more instruction about how to write a valid regex. Here studyID_regex asks for a regular expression that can capture all the input data files if there are multiple studies/projects involved, otherwise set it to be the same as "studyID".
#' @param studyID_altered logical, default FALSE, but if there are multiple studies/projects involved, should be set to TRUE so that all the involved studies/projects will be treated as a whole in an higher level (they will still remain separated in their individual level). See "Details" section for more details about "studyID", "studyID_regex" and "studyID_altered".
#' @param speciesID character, choose one from c("hg38", "hg19", "mm10"), default "hg38". This will affect downstream analysis methods based on different genome reference builds: "hg38" for GRCh38, "hg19" for GRCh37, "mm10" for GRCm38. In addition, each vigilante run should set only one "speciesID" (genome reference build). If input data files involve multiple genome reference build, separate them into multiple vigilante runs of which each uses a single "speciesID".
#' @param fileNum,fileDNA_mafNum integer, "fileNum" for the maximum total number of input data files per each sample (e.g. there are 50 samples, 35 of them have 5 files per sample, 12 of them have 7 files per sample, 3 of them have 6 files per sample, set "fileNum" to 7); similarly, "fileDNA_mafNum" for the maximum total number of .maf files per each sample (e.g. there are .maf files from Strelka came in a set of 2, while other .maf files from MuTect came in 1, set "fileDNA_mafNum" to 2).
#' @param clinicalFeature logical, whether input data have associated clinical information (e.g. Gleason score, Race, Ethnicity). If TRUE, addtional columns for specifying clinical information will be added to the vigilante-generated groupInfo.csv file and user needs to fill them before downstream analyses can properly take these clinical information into consideration.
#' @param createOutputFolders logical, whether to allow vigilante create specific output folders as the default place for storing downstream analyses output files (e.g. plots, results tables). If TRUE, a set of output folders will be created in the working directory under ./_VK/. It is recommended to set it to TRUE so that downstream analyses output files can be better organized inside "_VK" folder; if FALSE, user needs to specify a output path each time user chooses to generate a output file.
#' @param prepareVdata logical, whether to allow vigilante prepare all input data files in the working directory. This internal process includes several steps such as renaming data files, setting group-related parameters, setting ENSEMBL refernece, etc. If TRUE, user will be asked to backup input data files before continuing; if FALSE, vigilante will stop the run and no input data files will be affected.
#'
#' @details
#' The workflow of vigilante is highly module-based. Modules are connected and there are certain settings that are shared across all modules. To ensure a successful and smooth run, these settings need to be properly specified by the user.
#'
#' Take "studyID", "studyID_regex" and "studyID_altered" for example. Oftentimes input data files generated by upstream tools came with diverse naming conventions. It might be easy for the user to recognize those files, but not for vigilante if there is no recognizable patterns.
#'
#' To make input data files clear to vigilante, it would be nice to have them named something like "studyID_sampleID_(other descriptions).file extension". Here "studyID" is the name of the study or project, and it will be used in multiple naming situations (such as on the plot, or in the output file names), so it is recommended to be concise and meaningful.
#'
#' In the first demo example below, "studyID" is set to "KSCWUSCRF" (K - certain prefix, SCW - schwannoma, USC - University of Southern California, RF - certain suffix). Because it is a single study focusing on schwannoma and doesn't contain input data files from other studies (i.e. all input data files are named after "KSCWUSCRF" by upstream tools), "studyID_regex" is set to "[[:upper:]]{9}" to properly capture the study name pattern in the file names, or to make it simpler, here "studyID_regex" can be set to "KSCWUSCRF" as only one study is involved. Moreover, "studyID_altered" is set to FALSE because for single study, it is not necessary to change the study name inherited from upstream tools; however, if the user chooses to alter the study name, set "studyID_altered" to TRUE, and then "studyID" will be used as the new study name.
#'
#' Sometimes the working study isn't a single study but a combined study, and contains input data files from different sub-studies. In that case, "studyID_altered" should be set to TRUE so that all the involved sub-studies will be treated as a whole in an higher level (they will still remain separated in their individual level), and here "studyID_regex" should be a regular expression that can properly capture the study name patterns in the file names across all involved sub-studies.
#'
#' In the second demo example below, "studyID" is set to "KPROCOMBI" (K - certain prefix, PRO - prostate cancer, COMBI - combined study). Because sub-studies of "KPROCOMBI" are named such as "KGARUSCAG", "KPINSKIJC", "KRPLUSCJC" etc., here "studyID_regex" is set to "[[:upper:]]{9}" to properly capture them all.
#'
#' Another very important thing about vigilante is the "groupInfo.csv" file. This file contains the meta-info required by downstream analyses. By default, "groupInfo.csv" has five columns: "assayID", "Group", "MAF_group", "realID" and "aliasID." If "clinicalFeature" is set to TRUE, there can be more columns. Usually, user should leave "assayID", "MAF_group" and "realID" unchanged as they are auto-populated by vigilante and already in the right format. The "aliasID" column can be changed if user wants to set specific names for their samples, otherwise can be left unchanged as well. The "Group" column (and possible additional "CliFea" columns) is where user should properly fill in.
#'
#' For example, if samples are divided into training, validation and testing groups, the "Group" column should be filled with "Training", "Validation" and "Testing" accordingly. Similarly, if additional clinical information are available, user can specify them in the "CliFea" columns. Here "CliFea" is only a placeholder name, and user should change the column name to reflect that clinical feature (e.g. change it to "Race" column and fill in race information like "White", "Black or African American", "Asian" etc.; "Gleason Score" and fill in Gleason score values). Also, "CliFea" columns are not limited to two. User can add more columns to the right following the above instruction.
#'
#' @return list, because R CMD check discourages assignments to the global environment within functions, user needs to run the function with explicitly assigning the return value to a global variable named "globalSettings_returnList", which will be a list containing the required variables for downstream analyses.
#'
#' @examples
#' \dontrun{
#' # single study of schwannoma
#' globalSettings_returnList = v_globalSettings(studyID = "KSCWUSCRF",
#' studyID_regex = "[[:upper:]]{9}", studyID_altered = FALSE, speciesID =
#' "hg19", fileNum = 7, fileDNA_mafNum = 2, clinicalFeature = TRUE,
#' createOutputFolders = TRUE, prepareVdata = TRUE)
#'
#' # combined study of prostate cancer
#' globalSettings_returnList = v_globalSettings(studyID = "KPROCOMBI",
#' studyID_regex = "[[:upper:]]{9}", studyID_altered = TRUE, speciesID =
#' "hg19", fileNum = 7, fileDNA_mafNum = 2, clinicalFeature = TRUE,
#' createOutputFolders = TRUE, prepareVdata = TRUE)
#' }
#'
#' @import utils
#'
#' @export
#'
# v_globalSettings function
v_globalSettings = function(studyID, studyID_regex, studyID_altered = FALSE, speciesID = "hg38", fileNum, fileDNA_mafNum, clinicalFeature = FALSE, createOutputFolders = FALSE, prepareVdata = FALSE) {

  # ask user to make sure "vigilante.knights.sword" package is already installed
  status_sword = menu(choices = c("Yes", "No"), title = "\nDue to the requirement of CRAN that general packages should not exceed 5MB, the supplemental workbook or reference datasets (sword) required by vigilante & knights have been extracted and put in the standalone package vigilante.knights.sword. Please check https://github.com/yilixu/vigilante.knights.sword for more information. vigilante & knights will need 'sword' to perform downstream analysis. Is 'vigilante.knights.sword' package already installed?")
  if (status_sword == 1) {
    print("Checking 'vigilante.knights.sword' package status")
    status_sword_check = "vigilante.knights.sword" %in% utils::installed.packages()
    if (status_sword_check == TRUE) {
      print("Checking completed, now continue to the next step")
      rm(status_sword_check)
    } else {
      temp_install_choice1 = menu(choices = c("Yes (install from GitHub)", "No (stop and quit)"), title = "\n'vigilante.knights.sword' package not detected, do you want to re-install/update it now?")
      if (temp_install_choice1 == 1) {
        print("Start installing 'vigilante.knights.sword' from GitHub")
        devtools::install_github("yilixu/vigilante.knights.sword", ref = "main")
        print("Installation completed, now continue to the next step")
        rm(temp_install_choice1, status_sword_check)
      } else {
        print("Thank you for using vigilante & knights")
        stop()
      }
    }
  } else if (status_sword == 2) {
    temp_install_choice2 = menu(choices = c("Yes (install from GitHub)", "No (stop and quit)"), title = "\nDo you want to install 'vigilante.knights.sword' package now?")
    if (temp_install_choice2 == 1) {
      print("Start installing 'vigilante.knights.sword' from GitHub")
      devtools::install_github("yilixu/vigilante.knights.sword", ref = "main")
      print("Installation completed, now continue to the next step")
      rm(temp_install_choice2)
    } else {
      print("Thank you for using vigilante & knights")
      stop()
    }
  } else {
    print("Please choose a valid answer")
    stop()
  }
  rm(status_sword)

  # ask user to make sure v_globalSettings function return value is assigned to the global variable named "globalSettings_returnList"
  status_returnValue = menu(choices = c("Yes", "No"), title = "\nIs v_globalSettings function return value assigned to the global variable named 'globalSettings_returnList'?")
  if (status_returnValue == 1) {
    print("Function return value is properly assigned, now continue to the next step")
  } else if (status_returnValue == 2) {
    print("Please follow the instruction to assign the function return value before continue to the next step")
    stop()
  } else {
    print("Please choose a valid answer")
    stop()
  }
  rm(status_returnValue)

  # suppress experimental messages from dplyr.summarise
  options(dplyr.summarise.inform = FALSE)

  # set global variables required for downstream process
  temp_globalSettings_returnList = list(
    "studyID" = studyID,
    "studyID_regex" = studyID_regex,
    "studyID_altered" = studyID_altered,
    "speciesID" = speciesID,
    "fileNum" = fileNum,
    "fileDNA_mafNum" = fileDNA_mafNum,
    "clinicalFeature" = clinicalFeature
  )

  # create output folders
  if (createOutputFolders == TRUE) {
    if (!("_VK" %in% dir())) {
      print("Based on user's choice, a set of output folders will be created in the working directory under ./_VK/")
      dir.create("_data")
      dir.create("_data/_MuSiCa")
      dir.create("_VK")
      dir.create("_VK/_Circos")
      dir.create("_VK/_CHM")
      dir.create("_VK/_xCell")
      dir.create("_VK/_LLP")
      dir.create("_VK/_IGV")
      dir.create("_VK/_CSE")
      dir.create("_VK/_CSE/_Fusion")
      dir.create("_VK/_CSE/_Mutation")
      dir.create("_VK/_SSE")
      dir.create("_VK/_MuSiCa")
      dir.create("_VK/_DESeq2")
      dir.create("_VK/_GWAS")
      dir.create("_VK/_PWC")
      dir.create("_VK/_TCGA")
    } else {
      print("Output folders already exist, skipping to the next step")
    }
  } else {
    print("Based on user's choice, no output folders will be created in the file system")
  }

  # prepare data files for VIGILANTE format
  if (prepareVdata == TRUE) {
    prepareVdata_returnList = v_prepareVdata(studyID, studyID_regex, studyID_altered, speciesID, fileNum, fileDNA_mafNum, clinicalFeature)
  } else {
    print("In order to properly perform downstream analyses, vigilante needs to prepare all input data files in the working directory before continuing. Please set prepareVdata to TRUE and run the function again")
    stop()
  }

  # add variables from prepareVdata to global settings
  temp_globalSettings_returnList = c(temp_globalSettings_returnList, prepareVdata_returnList)

  # remove variables already added to temp_globalSettings_returnList
  suppressWarnings(rm(list = names(temp_globalSettings_returnList)))

  # assign globalSettings_returnList in the global environment
  # assign(x = "globalSettings_returnList", value = temp_globalSettings_returnList, envir = .GlobalEnv)

  # end of v_globalSettings function
  print("v_globalSettings run completed, return value saved to the global environment in **globalSettings_returnList**")
  return(temp_globalSettings_returnList)
}

# preset globalVariables for R CMD check
utils::globalVariables(c("ENSG", "Group", "MAF_Group"))
