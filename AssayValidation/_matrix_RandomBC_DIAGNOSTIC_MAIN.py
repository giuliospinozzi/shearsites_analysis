#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Mon Jan 18 09:44:56 2016

@author: stefano
"""


#++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import matrix_RandomBC_globModule
import matrix_RandomBC_assoModule
import matrix_RandomBC_dataModule
import matrix_RandomBC_processingModule
import matrix_RandomBC_filterModule
import matrix_RandomBC_outputModule

import matrix_RandomBC_CEMmodule
import matrix_RandomBC_barcodeDiagnosis


#++++++++++++++++++++++ Global Vars from globModule +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Import Functions
humanSorted = matrix_RandomBC_globModule.humanSorted

### Screen print configs 
verbose = matrix_RandomBC_globModule.verbose
verbosePrint = matrix_RandomBC_globModule.verbosePrint

### Association File Loading configs - assoModule
asso_folder = matrix_RandomBC_globModule.asso_folder
asso_file_name = matrix_RandomBC_globModule.asso_file_name
asso_delimiter = matrix_RandomBC_globModule.asso_delimiter

### Data Loading configs - dataModule
ground_dir = matrix_RandomBC_globModule.ground_dir
DISEASE = matrix_RandomBC_globModule.DISEASE
PATIENT = matrix_RandomBC_globModule.PATIENT
POOL = matrix_RandomBC_globModule.POOL
data_files_delimiter = matrix_RandomBC_globModule.data_files_delimiter
data_files_name_filter = matrix_RandomBC_globModule.data_files_name_filter

### Filter Data configs - filterModule
filter_data = matrix_RandomBC_globModule.filter_data
byHeaders = matrix_RandomBC_globModule.byHeaders
bySC = matrix_RandomBC_globModule.bySC

### COMMON OUTPUT GROUND DIR
common_output_ground_dir = matrix_RandomBC_globModule.common_output_ground_dir

### Matrix output configs - outputModule
matrix_outfolder = matrix_RandomBC_globModule.matrix_outfolder
matrix_files_delimiter = matrix_RandomBC_globModule.matrix_files_delimiter
relabelling = matrix_RandomBC_globModule.relabelling

### Export CEM Data ?
export_cem_data = matrix_RandomBC_globModule.export_cem_data
### Export CEM matrixes ?
export_matrixes = matrix_RandomBC_globModule.export_matrixes
### Export CEM Diagnostics ?
export_diagnostics = matrix_RandomBC_globModule.export_diagnostics
diagnostic_outfolder = matrix_RandomBC_globModule.diagnostic_outfolder
# Data selection configs
specific_samples = matrix_RandomBC_globModule.specific_samples
dilution_to_process = matrix_RandomBC_globModule.dilution_to_process
condition_to_process = matrix_RandomBC_globModule.condition_to_process
# Task to perform configs
checkNucleotidesBalancing = matrix_RandomBC_globModule.checkNucleotidesBalancing
FragmentLengthDistribution = matrix_RandomBC_globModule.FragmentLengthDistribution
checkShearSitesOccurrency = matrix_RandomBC_globModule.checkShearSitesOccurrency
checkRandomBCoccurrency = matrix_RandomBC_globModule.checkRandomBCoccurrency
checkEditDistance_diagonal = matrix_RandomBC_globModule.checkEditDistance_diagonal
checkEditDistance_extensive = matrix_RandomBC_globModule.checkEditDistance_extensive
plot_heatmap = matrix_RandomBC_globModule.plot_heatmap
limit_heatmap_plot = matrix_RandomBC_globModule.limit_heatmap_plot
plot_heatmap_byChunks = matrix_RandomBC_globModule.plot_heatmap_byChunks
ShS_chunk_size = matrix_RandomBC_globModule.ShS_chunk_size
checkBCcountRatio = matrix_RandomBC_globModule.checkBCcountRatio

#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


### Load Association File Data ###################################################################
asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
##################################################################################################

### Load Data ################################################################################################################################################
POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
##############################################################################################################################################################

### Build Data ##################################################################
verbosePrint("\n>>> Shape data as DataFrame ...")
#df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
df = matrix_RandomBC_processingModule.buildExhaustiveDataFrame(POOL_alldata_dict)
verbosePrint(">>> Dataframe built!")
#################################################################################

### Filter Data #############################################################################################
if filter_data:
    verbosePrint("\n>>> Filtering DataFrame ...")
    if byHeaders:
        headers_file_dir = matrix_RandomBC_filterModule.buildInputPath(ground_dir, DISEASE, PATIENT)
        lane_ID = "{POOL}".format(POOL=str(POOL)).lower().replace("_", "")
        headers_file_name = "{lane_ID}.quality.r1r2-100rbc.toRemove.list.gz".format(lane_ID=str(lane_ID))
        headers_file_path = os.path.normpath(os.path.join(headers_file_dir, headers_file_name))
        headers_to_remove = matrix_RandomBC_filterModule.loadGzFile(headers_file_path)
        df = matrix_RandomBC_filterModule.filterDF_byHeaders(df, headers_to_remove)
    if bySC:
        pass
    verbosePrint(">>> Done!")
#############################################################################################################

### Process Data as Matrixes #####################################################################################################
seqCount_matrix, ShsCount_matrix, barcodeCount_matrix, cellCount_matrix, fragmentEstimate_matrix = None, None, None, None, None
if (export_cem_data or export_matrixes):
    verbosePrint("\n>>> Computing matrixes:")
    verbosePrint("> seqCount matrix ...")
    seqCount_matrix = matrix_RandomBC_processingModule.buildSeqCountMatrix(df)
    verbosePrint("> ShsCount matrix ...")
    ShsCount_matrix = matrix_RandomBC_processingModule.buildShsCountMatrix(df)
    verbosePrint("> barcodeCount matrix ...")
    barcodeCount_matrix = matrix_RandomBC_processingModule.buildBarcodeCountMatrix(df)
    verbosePrint("> cellCount matrix ...")
    cellCount_matrix = matrix_RandomBC_processingModule.buildCellCountMatrix(df)
    verbosePrint("> fragmentEstimate matrix ...")
    fragmentEstimate_matrix = matrix_RandomBC_processingModule.buildFragmentEstimateMatrix(df)
    verbosePrint(">>> Done!")
#################################################################################################################################


### Export ###################################################################################################################################################

### EXPORT MATRIXES
if export_matrixes:
    # Set column relabelling
    metadata = None
    if relabelling:
        metadata = asso_dict
    # Export: outdir -> outPath -> writeMatrix
    verbosePrint("\n>>> Export matrixes ...")
    OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(common_output_ground_dir, matrix_outfolder)
    seqCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "seqCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(seqCount_matrix, seqCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
    ShsCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "ShsCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(ShsCount_matrix, ShsCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
    barcodeCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "barcodeCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(barcodeCount_matrix, barcodeCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
    cellCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "cellCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(cellCount_matrix, cellCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
    fragmentEstimate_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "fragmentEstimate_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(fragmentEstimate_matrix, fragmentEstimate_matrix_outPath, matrix_files_delimiter, metadata=metadata)
    verbosePrint(">>> Matrix Files Created!")

### EXPORT CEM DATA
if export_cem_data:
    verbosePrint("\n>>> Export CEM data ...")
    # Normalize matrixes
    seqCount_matrix_norm = seqCount_matrix.apply(lambda x: x/x.sum())
    ShsCount_matrix_norm = ShsCount_matrix.apply(lambda x: x/x.sum())
    barcodeCount_matrix_norm = barcodeCount_matrix.apply(lambda x: x/x.sum())
    cellCount_matrix_norm = cellCount_matrix.apply(lambda x: x/x.sum())
    fragmentEstimate_matrix_norm = fragmentEstimate_matrix.apply(lambda x: x/x.sum())
    # Prepare: more matrixes can be added!
    df_dict = {'sequence_count': seqCount_matrix, 'shearsite_count': ShsCount_matrix, 'TAG_count': barcodeCount_matrix, 'shearsite&TAG_count': cellCount_matrix, 'fragmentEstimate_count': fragmentEstimate_matrix,
               'sequence_count_quantification': seqCount_matrix_norm, 'shearsite_count_quantification': ShsCount_matrix_norm, 'TAG_count_quantification': barcodeCount_matrix_norm, 'shearsite&TAG_count_quantification': cellCount_matrix_norm, 'fragmentEstimate_count_quantification': fragmentEstimate_matrix_norm}
    # Export
    outpath = matrix_RandomBC_CEMmodule.exportCEM(df_dict, asso_dict)
    verbosePrint(">>> CEM data exported!")
    
### EXPORT DIAGNOSTICS

## Rebuild df for deeper diagnostics! NOT REQUIRED
#verbosePrint("\n>>> Rebuild and extend DataFrame for deeper diagnostics ...")
#df = matrix_RandomBC_processingModule.buildExhaustiveDataFrame(POOL_alldata_dict)
#verbosePrint(">>> Done! Headers loaded.")

if export_diagnostics:
    verbosePrint("\n\n *** START DIAGNOSTICS ***")
    # Prelimiary operations
    OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(common_output_ground_dir, diagnostic_outfolder)
    cem_IDs_extended = []  # it implements 'position tolerance' through matrix_RandomBC_CEMmodule.cem_locations
    for t, s in zip(matrix_RandomBC_CEMmodule.cem_locations, matrix_RandomBC_CEMmodule.cem_strands):
        cem_IDs_extended += matrix_RandomBC_CEMmodule.buildIndexes(t, s)
    # Get barcodes (and human Labels) and loop into!
    sampleLabelsDict = matrix_RandomBC_CEMmodule.buildRelabellingDict(asso_dict, "_")  # key=barcode, item=humanLabel _-separated
    sample_barcodes = humanSorted(df['barcode'].unique())
    for barcode in sample_barcodes:
        sample_label = sampleLabelsDict[barcode]
        verbosePrint("\n* Processing barcode {barcode} - {label}:".format(barcode=barcode, label=sample_label))
        if specific_samples:
            # only condition_to_process
            if asso_dict[barcode][2] not in condition_to_process:
                verbosePrint("  ...SKIP! (not in condition_to_process={condition_to_process})".format(condition_to_process=str(condition_to_process)))
                continue
            # only dilution_to_process
            if asso_dict[barcode][4] not in dilution_to_process:
                verbosePrint("  ...SKIP! (not in dilution_to_process={dilution_to_process})".format(dilution_to_process=str(dilution_to_process)))
                continue
        # Slice df by barcode -> barcode_df
        barcode_df = df[df['barcode']==barcode]
        # Prepare sample_OUTDIR
        sample_OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(OUTDIR, barcode+"_"+sample_label)
        # checkNucleotidesBalancing
        if checkNucleotidesBalancing:
            # By SC
            verbosePrint("  > Check Nucleotides Balancing by SC ...")
            barcode_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(barcode_df, how='by_seq_count')
            verbosePrint("  > Plot Nucleotides Balancing by SC ...")
            title = "PILED-UP RANDOM-BARCODES by SC - {barcode} {label}".format(barcode=barcode, label=sample_label)
            filename = "{barcode}_{label}_nucleotideBalancing_bySC.pdf".format(barcode=barcode, label=sample_label)
            export = os.path.normpath(os.path.join(sample_OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(barcode_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
            # By distinct-over-ShS
            verbosePrint("  > Check Nucleotides Balancing by distinct-over-ShS ...")
            barcode_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(barcode_df, how='simple')
            verbosePrint("  > Plot Nucleotides Balancing  by distinct-over-ShS ...")
            title = "PILED-UP RANDOM-BARCODES by distinct-over-ShS - {barcode} {label}".format(barcode=barcode, label=sample_label)
            filename = "{barcode}_{label}_nucleotideBalancing_by-distinct-over-ShS.pdf".format(barcode=barcode, label=sample_label)
            export = os.path.normpath(os.path.join(sample_OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(barcode_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
            # By distinct-overall
            verbosePrint("  > Check Nucleotides Balancing by distinct-overall ...")
            barcode_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(barcode_df, how='distinct')
            verbosePrint("  > Plot Nucleotides Balancing by distinct-overall ...")
            title = "PILED-UP RANDOM-BARCODES by distinct-overall - {barcode} {label}".format(barcode=barcode, label=sample_label)
            filename = "{barcode}_{label}_nucleotideBalancing_by-distinct-overall.pdf".format(barcode=barcode, label=sample_label)
            export = os.path.normpath(os.path.join(sample_OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(barcode_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
        # FragmentLengthDistribution
        if FragmentLengthDistribution:
            verbosePrint("  > Check Fragment Length Distribution ...")
            barcode_FragmentLengthDistribution_DF = matrix_RandomBC_barcodeDiagnosis.checkFragmentLengthDistribution(barcode_df)
            verbosePrint("  > Plot Fragment Length Distribution ...")
            title = "FRAGMENT LENGTH DISTRIBUTION - {barcode} {label}".format(barcode=barcode, label=sample_label)
            filename = "{barcode}_{label}_fragmentLengthDistribution.pdf".format(barcode=barcode, label=sample_label)
            export = os.path.normpath(os.path.join(sample_OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotFragmentLengthDistribution(barcode_FragmentLengthDistribution_DF, title=title, binning='Freedman-Diaconis', normalize=True, rug=False, kde=True, sonicLengthEstimation=True, show_live=False, export=export)
        # Get ISs in this sample to loop into!
        IS_IDs = humanSorted(barcode_df['genomic_coordinates'].unique())
        # Search for CEM data
        verbosePrint("  > Search for CEM data inside this sample ...")
        CEM_found = humanSorted(list(set(cem_IDs_extended).intersection(set(IS_IDs))))
        if CEM_found == list():
            verbosePrint("    Not Found, SKIP!")
            continue
        else:
            verbosePrint("    Found! {IS_found}.".format(IS_found=str(CEM_found)))
            # Restrict barcode_df/IS_IDs to CEM data: OVERRIDE
            verbosePrint("    Extracting Data...")
            N = len(barcode_df)
            barcode_df = matrix_RandomBC_CEMmodule.extractCEM_fromDF(barcode_df)
            verbosePrint("    Done! {n} rows taken over {N}.".format(n=str(len(barcode_df)), N=str(N)))
            IS_IDs = humanSorted(barcode_df['genomic_coordinates'].unique())
            verbosePrint("    Check CEM ISs found: {IS_IDs}.".format(IS_IDs=str(IS_IDs)))
        for IS in IS_IDs:
            # Get CEM info
            CEM_info_dict = matrix_RandomBC_CEMmodule.CEM_info(IS)
            CEM_real_coordinate = IS
            CEM_nominal_coordinate = CEM_info_dict['coordinate']
            CEM_name = CEM_info_dict['symbol']
            verbosePrint("    > Processing {CEM_name}:{CEM_nominal_coordinate} (coordinate found:{CEM_real_coordinate})".format(CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate)))
            # Slice barcode_df by IS -> IS_df
            IS_df = barcode_df[barcode_df['genomic_coordinates']==IS]
            # Prepare IS_OUTDIR
            IS_OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(sample_OUTDIR, "{CEM_name}:{CEM_nominal_coordinate}".format(CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate)))
            # checkShearSitesOccurrency
            if checkShearSitesOccurrency:
                verbosePrint("      > Check Shear Sites Occurrency ...")
                IS_ShearSitesOccurrency_DF = matrix_RandomBC_barcodeDiagnosis.checkShearSitesOccurrency(IS_df)
                verbosePrint("      > Plot Shear Sites Occurrency ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} SHEAR SITE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_shearSitesOccurrency.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotShearSitesOccurrency(IS_ShearSitesOccurrency_DF, title=title, normalize=100, show_live=False, export=export)
            # checkNucleotidesBalancing
            if checkNucleotidesBalancing:
                verbosePrint("      > Check Nucleotides Balancing ...")
                IS_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(IS_df, how='simple')
                verbosePrint("      > Plot Nucleotides Balancing ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} PILED-UP RANDOM-BARCODES by-distinct-ShS - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_nucleotideBalancing_by-distinct-ShS.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(IS_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
            # checkRandomBCoccurrency
            if checkRandomBCoccurrency:
                verbosePrint("      > Check Random Barcodes Occurrency ...")
                IS_distinctBC_DF = matrix_RandomBC_barcodeDiagnosis.checkRandomBCoccurrency(IS_df)
                verbosePrint("      > Plot Random Barcodes Occurrency ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_RandomBCoccurrency.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotRandomBCoccurrency(IS_distinctBC_DF, title=title, show_top_ranked=10, annot=True, show_live=False, export=export)
            # checkEditDistance - diagonal
            if checkEditDistance_diagonal:
                verbosePrint("      > Check Edit Distances whitin Shear Sites ...")
                IS_diagonal_editDistance_DF = matrix_RandomBC_barcodeDiagnosis.checkEditDistance(IS_df, all_combinations=False)
                verbosePrint("        Done! Square matrix {dim}X{dim}.".format(dim=str(len(IS_diagonal_editDistance_DF))))
                if plot_heatmap:
                    if limit_heatmap_plot[0]:
                        if len(IS_diagonal_editDistance_DF) <= limit_heatmap_plot[1]:
                            verbosePrint("      > Plot Edit Distance Heatmap whitin Shear Sites ...")
                            title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                            export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                            matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_diagonal_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
                        else:
                            verbosePrint("      [SKIP] Edit Distance Heatmap Plot! Limit {lim}X{lim}.".format(lim=str(limit_heatmap_plot[1])))
                    else:
                        verbosePrint("      > Plot Edit Distance Heatmap whitin Shear Sites ...")
                        title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                        matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_diagonal_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
                verbosePrint("      > Plot Edit Distance Distribution ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} EDIT DISTANCE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceDistribution.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotEditDistanceOccurrency(IS_diagonal_editDistance_DF, title=title, vmin=0, vmax=12, annot=True, percentile_colors=False, show_live=False, export=export)
                # chunking
                if plot_heatmap_byChunks:
                    verbosePrint("      > Trying to chunck Edit Distance Matrix (chunk size={ShS_chunk_size})...".format(ShS_chunk_size=str(ShS_chunk_size)))
                    IS_diagonal_EditDistance_DFchunks = matrix_RandomBC_barcodeDiagnosis.chunkEditDistance_DF(IS_diagonal_editDistance_DF, ShS_chunk_size=ShS_chunk_size)
                    verbosePrint("        Done! {N} chunk(s).".format(N=str(len(IS_diagonal_EditDistance_DFchunks))))
                    verbosePrint("      > Plot Edit Distance Heatmap chunk-by-chunk...")
                    EDheatmap_byChunk_OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(IS_OUTDIR, "EDheatmap_byChunks")
                    for n, ED_DF in list(enumerate(IS_diagonal_EditDistance_DFchunks, start=1)):
                        ID = "CHUNK{n}".format(n=str(n))
                        verbosePrint("        > {ID} ...".format(ID=ID))
                        title = "{ID} - {CEM_name}:{CEM_nominal_coordinate} - {barcode} {label} {CEM_real_coordinate}".format(ID=ID, barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        filename = "{ID}_{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap.pdf".format(ID=ID, barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        export = os.path.normpath(os.path.join(EDheatmap_byChunk_OUTDIR, filename))
                        matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(ED_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
            # checkEditDistance - extensive
            if checkEditDistance_extensive:
                verbosePrint("      > Check Edit Distances all-VS-all ...")
                IS_extensive_editDistance_DF = matrix_RandomBC_barcodeDiagnosis.checkEditDistance(IS_df, all_combinations=True)
                verbosePrint("        Done! Square matrix {dim}X{dim}.".format(dim=str(len(IS_extensive_editDistance_DF))))
                if plot_heatmap:
                    if limit_heatmap_plot[0]:
                        if len(IS_extensive_editDistance_DF) <= limit_heatmap_plot[1]:
                            verbosePrint("      > Plot Edit Distance Heatmap all-VS-all...")
                            title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX all-VS-all - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap_ALL.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                            export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                            matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_extensive_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
                        else:
                            verbosePrint("      [SKIP] Edit Distance Heatmap Plot! Limit {lim}X{lim}.".format(lim=str(limit_heatmap_plot[1])))
                    else:
                        verbosePrint("      > Plot Edit Distance Heatmap all-VS-all...")
                        title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX all-VS-all - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap_ALL.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                        matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_extensive_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
                verbosePrint("      > Plot Edit Distance Distribution all-VS-all ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} EDIT DISTANCE OCCURRENCIES all-VS-all - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceDistribution_ALL.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotEditDistanceOccurrency(IS_extensive_editDistance_DF, title=title, vmin=0, vmax=12, annot=True, percentile_colors=False, show_live=False, export=export)
                # chunking
                if plot_heatmap_byChunks:
                    verbosePrint("      > Trying to chunck Edit Distance Matrix (chunk size={ShS_chunk_size})...".format(ShS_chunk_size=str(ShS_chunk_size)))
                    IS_extensive_EditDistance_DFchunks = matrix_RandomBC_barcodeDiagnosis.chunkEditDistance_DF(IS_extensive_editDistance_DF, ShS_chunk_size=ShS_chunk_size)
                    verbosePrint("        Done! {N} chunk(s).".format(N=str(len(IS_extensive_EditDistance_DFchunks))))
                    verbosePrint("      > Plot Edit Distance Heatmap chunk-by-chunk...")
                    EDheatmap_ALL_byChunk_OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(IS_OUTDIR, "EDheatmap_ALL_byChunks")
                    for n, ED_DF in list(enumerate(IS_extensive_EditDistance_DFchunks, start=1)):
                        ID = "CHUNK{n}".format(n=str(n))
                        verbosePrint("        > {ID} ...".format(ID=ID))
                        title = "{ID} - {CEM_name}:{CEM_nominal_coordinate} - {barcode} {label} {CEM_real_coordinate}".format(ID=ID, barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        filename = "{ID}_{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap_ALL.pdf".format(ID=ID, barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                        export = os.path.normpath(os.path.join(EDheatmap_ALL_byChunk_OUTDIR, filename))
                        matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(ED_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
            # checkBCcountRatio - 'diagonal'
            if checkBCcountRatio:
                verbosePrint("      > Check Random Barcodes seqCount ratios ...")
                IS_BCcountRatio_DF = matrix_RandomBC_barcodeDiagnosis.checkBCcountRatio(IS_df)  # all_combinations=False by default
                verbosePrint("      > Plot Random Barcodes seqCount ratios ...")
                title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODE SEQ-COUNT RATIOS - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_RandomBCseqCountRatios.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
                export = os.path.normpath(os.path.join(IS_OUTDIR, filename))
                matrix_RandomBC_barcodeDiagnosis.plotBCcountRatio(IS_BCcountRatio_DF, title=title, show_live=False, export=export)

            ###################################
            # HERE FURTHER DIAGNS AT IS-LEVEL #
            ###################################
            
            #shearsites = humanSorted(IS_df['shearsite'].unique())
            #for ShS in shearsites:
                #ShS_df = IS_df[IS_df['shearsite']==ShS]
                #randomBC_list = list(ShS_df['randomBC'])
                #sc_list = list(ShS_df['seq_count'])
                
                ########################################################
                # HERE FURTHER DIAGNS AT SAMPLE-CEM_IS-SHEARSITE LEVEL #
                ########################################################
    
    verbosePrint("\n *** END DIAGNOSTICS ***\n")


#############################################################################################################################################################


