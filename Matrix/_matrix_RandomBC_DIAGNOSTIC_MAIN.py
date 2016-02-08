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
import matrix_RandomBC_outputModule

import matrix_RandomBC_CEMmodule
import matrix_RandomBC_barcodeDiagnosis


#++++++++++++++++++++++ Global Vars from globModule +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Import Functions
humanSorted = matrix_RandomBC_globModule.humanSorted

# Screen print
verbose = matrix_RandomBC_globModule.verbose
verbosePrint = matrix_RandomBC_globModule.verbosePrint

# Association File - assoModule
asso_folder = matrix_RandomBC_globModule.asso_folder
asso_file_name = matrix_RandomBC_globModule.asso_file_name
asso_delimiter = matrix_RandomBC_globModule.asso_delimiter

# Data - dataModule
ground_dir = matrix_RandomBC_globModule.ground_dir
DISEASE = matrix_RandomBC_globModule.DISEASE
PATIENT = matrix_RandomBC_globModule.PATIENT
POOL = matrix_RandomBC_globModule.POOL
data_files_delimiter = matrix_RandomBC_globModule.data_files_delimiter
data_files_name_filter = matrix_RandomBC_globModule.data_files_name_filter

# Output - outputModule
#ground_dir, DISEASE, PATIENT, POOL as for Data
outfolder = matrix_RandomBC_globModule.outfolder
out_files_delimiter = matrix_RandomBC_globModule.out_files_delimiter
relabelling = matrix_RandomBC_globModule.relabelling

# Export CEM Data
export_cem_data = matrix_RandomBC_globModule.export_cem_data
# Export matrixes
export_matrixes = matrix_RandomBC_globModule.export_matrixes
# Export Diagnostics
export_diagnostics = matrix_RandomBC_globModule.export_diagnostics


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


### Load Association File Data ###################################################################
asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
##################################################################################################

### Load Data ################################################################################################################################################
POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
##############################################################################################################################################################

### Process Data ################################################################################################################
df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
seqCount_matrix, ShsCount_matrix, barcodeCount_matrix, cellCount_matrix, fragmentEstimate_matrix = None, None, None, None, None
if (export_cem_data or export_matrixes):
    seqCount_matrix = matrix_RandomBC_processingModule.buildSeqCountMatrix(df)
    ShsCount_matrix = matrix_RandomBC_processingModule.buildShsCountMatrix(df)
    barcodeCount_matrix = matrix_RandomBC_processingModule.buildBarcodeCountMatrix(df)
    cellCount_matrix = matrix_RandomBC_processingModule.buildCellCountMatrix(df)
    fragmentEstimate_matrix = matrix_RandomBC_processingModule.buildFragmentEstimateMatrix(df)
#################################################################################################################################


### Export ###################################################################################################################################################

### EXPORT MATRIXES
if export_matrixes:
    # Set column relabelling
    metadata = None
    if relabelling:
        metadata = asso_dict
    # Export: outdir -> outPath -> writeMatrix
    OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(ground_dir, DISEASE, PATIENT, POOL, outfolder)
    seqCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "seqCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(seqCount_matrix, seqCount_matrix_outPath, out_files_delimiter, metadata=metadata)
    ShsCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "ShsCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(ShsCount_matrix, ShsCount_matrix_outPath, out_files_delimiter, metadata=metadata)
    barcodeCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "barcodeCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(barcodeCount_matrix, barcodeCount_matrix_outPath, out_files_delimiter, metadata=metadata)
    cellCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "cellCount_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(cellCount_matrix, cellCount_matrix_outPath, out_files_delimiter, metadata=metadata)
    fragmentEstimate_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "fragmentEstimate_matrix.tsv"))
    matrix_RandomBC_outputModule.writeMatrix(fragmentEstimate_matrix, fragmentEstimate_matrix_outPath, out_files_delimiter, metadata=metadata)

### EXPORT CEM DATA
if export_cem_data:
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
    
### EXPORT DIAGNOSTICS

if export_diagnostics:
    # Prelimiary operations
    dilution_to_process = ['L', 'M', 'N']
    OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(ground_dir, DISEASE, PATIENT, POOL, "Diagnostics")
    cem_IDs_extended = []
    for t, s in zip(matrix_RandomBC_CEMmodule.cem_locations, matrix_RandomBC_CEMmodule.cem_strands):
        cem_IDs_extended += matrix_RandomBC_CEMmodule.buildIndexes(t, s)
    # Get barcodes (and human Labels) and loop into!
    sampleLabelsDict = matrix_RandomBC_CEMmodule.buildRelabellingDict(asso_dict, "_")  # key=barcode, item=humanLabel _-separated
    sample_barcodes = humanSorted(df['barcode'].unique())
    for barcode in sample_barcodes:
        sample_label = sampleLabelsDict[barcode]
        verbosePrint("\n* Processing barcode {barcode} - {label}:".format(barcode=barcode, label=sample_label))
        # L, M, N dilutions only
        if asso_dict[barcode][4] not in dilution_to_process:
            verbosePrint("  ...SKIP! (not in dilution_to_process={dilution_to_process})".format(dilution_to_process=str(dilution_to_process)))
            continue
        ###################################################
        # HERE SELECT ALSO BY OTHER FEATURES ##############
        #        ...TO DO!!!                 ##############
        ###################################################
        # Slice df by barcode -> barcode_df
        barcode_df = df[df['barcode']==barcode]
        # checkNucleotidesBalancing
        verbosePrint("  > Check Nucleotides Balancing ...")
        barcode_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(barcode_df)
        verbosePrint("  > Plot Nucleotides Balancing ...")
        title = "PILED-UP RANDOM-BARCODES - {barcode} {label}".format(barcode=barcode, label=sample_label)
        filename = "{barcode}_{label}_nucleotideBalancing.pdf".format(barcode=barcode, label=sample_label)
        export = os.path.normpath(os.path.join(OUTDIR, filename))
        matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(barcode_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
        # FragmentLengthDistribution
        verbosePrint("  > Check Fragment Length Distribution ...")
        barcode_FragmentLengthDistribution_DF = matrix_RandomBC_barcodeDiagnosis.checkFragmentLengthDistribution(barcode_df)
        verbosePrint("  > Plot Fragment Length Distribution ...")
        title = "FRAGMENT LENGTH DISTRIBUTION - {barcode} {label}".format(barcode=barcode, label=sample_label)
        filename = "{barcode}_{label}_fragmentLengthDistribution.pdf".format(barcode=barcode, label=sample_label)
        export = os.path.normpath(os.path.join(OUTDIR, filename))
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
            barcode_df = matrix_RandomBC_CEMmodule.extractCEM(barcode_df)
            IS_IDs = humanSorted(barcode_df['genomic_coordinates'].unique())
        for IS in IS_IDs:
            # Get CEM info
            CEM_info_dict = matrix_RandomBC_CEMmodule.CEM_info(IS)
            CEM_real_coordinate = IS
            CEM_nominal_coordinate = CEM_info_dict['coordinate']
            CEM_name = CEM_info_dict['symbol']
            verbosePrint("    > Processing {CEM_name}:{CEM_nominal_coordinate} (coordinate found:{CEM_real_coordinate})".format(CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate)))
            # Slice barcode_df by IS -> IS_df
            IS_df = barcode_df[barcode_df['genomic_coordinates']==IS]
            # checkShearSitesOccurrency
            verbosePrint("      > Check Shear Sites Occurrency ...")
            IS_ShearSitesOccurrency_DF = matrix_RandomBC_barcodeDiagnosis.checkShearSitesOccurrency(IS_df)
            verbosePrint("      > Plot Shear Sites Occurrency ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} SHEAR SITE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_shearSitesOccurrency.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotShearSitesOccurrency(IS_ShearSitesOccurrency_DF, title=title, normalize=100, show_live=False, export=export)
            # checkNucleotidesBalancing
            verbosePrint("      > Check Nucleotides Balancing ...")
            IS_nucleotidesCount_DF = matrix_RandomBC_barcodeDiagnosis.checkNucleotideBalancing(IS_df)
            verbosePrint("      > Plot Nucleotides Balancing ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} PILED-UP RANDOM-BARCODES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_nucleotideBalancing.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotNucleotideBalancing(IS_nucleotidesCount_DF, title=title, stacked_bar=True, show_live=False, export=export)
            # checkRandomBCoccurrency
            verbosePrint("      > Check Random Barcodes Occurrency ...")
            IS_distinctBC_DF = matrix_RandomBC_barcodeDiagnosis.checkRandomBCoccurrency(IS_df)
            verbosePrint("      > Plot Random Barcodes Occurrency ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_RandomBCoccurrency.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotRandomBCoccurrency(IS_distinctBC_DF, title=title, show_top_ranked=10, annot=True, show_live=False, export=export)
            # checkEditDistance - diagonal
            verbosePrint("      > Check Edit Distances whitin Shear Sites ...")
            IS_diagonal_editDistance_DF = matrix_RandomBC_barcodeDiagnosis.checkEditDistance(IS_df, all_combinations=False)
            verbosePrint("      > Plot Edit Distance Heatmap whitin Shear Sites ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_diagonal_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
            verbosePrint("      > Plot Edit Distance Distribution ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} EDIT DISTANCE OCCURRENCIES - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceDistribution.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotEditDistanceOccurrency(IS_diagonal_editDistance_DF, title=title, vmin=0, vmax=12, annot=True, percentile_colors=False, show_live=False, export=export)
            # checkEditDistance - extensive
            ###################################################
            # HERE PUT SOME CONTROLS             ##############
            #        ...TO DO!!!                 ##############
            ###################################################
            verbosePrint("      > Check Edit Distances all-VS-all ...")
            IS_extensive_editDistance_DF = matrix_RandomBC_barcodeDiagnosis.checkEditDistance(IS_df, all_combinations=True)
            verbosePrint("      > Plot Edit Distance Heatmap all-VS-all...")
            title = "{CEM_name}:{CEM_nominal_coordinate} RANDOM-BARCODES EDIT-DISTANCE MATRIX all-VS-all - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceHeatmap_ALL.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.editDistanceHeatmap(IS_extensive_editDistance_DF, title=title, cmap="RdYlBu", annot=True, show_live=False, export=export)
            verbosePrint("      > Plot Edit Distance Distribution all-VS-all ...")
            title = "{CEM_name}:{CEM_nominal_coordinate} EDIT DISTANCE OCCURRENCIES all-VS-all - {barcode} {label} {CEM_real_coordinate}".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            filename = "{barcode}_{label}_{CEM_name}:{CEM_nominal_coordinate}_{CEM_real_coordinate}_EditDistanceDistribution_ALL.pdf".format(barcode=barcode, label=sample_label, CEM_name=str(CEM_name), CEM_nominal_coordinate=str(CEM_nominal_coordinate), CEM_real_coordinate=str(CEM_real_coordinate))
            export = os.path.normpath(os.path.join(OUTDIR, filename))
            matrix_RandomBC_barcodeDiagnosis.plotEditDistanceOccurrency(IS_extensive_editDistance_DF, title=title, vmin=0, vmax=12, annot=True, percentile_colors=False, show_live=False, export=export)
            ###################################################
            # IF CONTROLS FAIL, try:                          #
            # - Only hist                                     #
            # - masked heatmap (show 0,1,2 only)              #
            #                                     ...TO DO!!! #
            ###################################################
            
            #shearsites = humanSorted(IS_df['shearsite'].unique())
            #for ShS in shearsites:
                #ShS_df = IS_df[IS_df['shearsite']==ShS]
                #randomBC_list = list(ShS_df['randomBC'])
                #sc_list = list(ShS_df['seq_count'])
                #### HERE SOMETHING AT SAMPLE-CEM_IS-SHEARSITE LEVEL
        


#############################################################################################################################################################


