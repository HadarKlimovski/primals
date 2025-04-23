import os
import sys
sys.path.insert(0, '/home/labs/hornsteinlab/hadarkl/Hadar/primals/code')


SELECTED_FEATURES_COX_OLD = ['LRG1',
 'GSN',
 #'ITIH4.1',
 'AGT',
 'IGF1',
 'REG1A',
 'CTSB',
 'IGFBP1',
 'VCAN',
 'FABP4',
 'GHR',
 'IGFBP2',
 'MRC1',
 'RNH1',
 'IGFALS',
 'NBL1',
 'NPC2',
 'FGL1',
 'PTPRJ',
 'ADAM15',
 'APOF',
 'CILP2',
 'SPON1',
 'FETUB',
 'VSIG4',
 'ALSFRS score (unit)',
 'ALSFRS_slope (-unit/month)']
#protein that chosen with replication ci=0.65 from the iteretaion algorithem
ITA_33_1 = ['NPC2', 'ADAM15', 'PTPRJ_CTNNB1', 'ALSFRS_slope (-unit/month)', 'LRG1', 'VSIG4_C3', 'CTSB', 'MT2A', 'FETUB_CETP', 'FABP4_VIM', 'ICOSLG', 'NAGLU', 'APOF_APOA4', 'ALSFRS score (unit)', 'AGT', 'LRG1_CYCS', 'CDH1', 'REG1A', 'PON1', 'NPC2_TTR']

SELECTED_FEATURES_COX_REP = ['ICOSLG',
 'AGT',
 'LRG1',
 'MT2A',
 'SHBG',
 'IGF1',
 'REG1A',
 'CTSB',
 'GHR',
 'CDH1',
 'FABP4',
 'IGFBP2',
 'MRC1',
 'PON1',
 'GRN',
 'CASP1',
 'IGFALS',
 'NBL1',
 'NAGLU',
 'NPC2',
 'FKBP1A',
 'FGL1',
 'PTPRJ',
 'ADAM15',
 'ART3',
 'APOF',
 'CSPG4',
 'CILP2',
 'IGDCC4',
 'SPON1',
 'APMAP',
 'TNFRSF12A',
 'CNTN3',
 'EXTL2',
 'FETUB',
 'VSIG4',
 'LGALS3',
 'Neutrophil degranulation',
 'ALSFRS_slope (-unit/month)',
 'STAT3 Pathway',
 'ALSFRS score (unit)',
 'smok p*y',
 '1 = Age Onset ≥ 60']
 
#-------------------------------------------------------------------------------------------------------------------------------------------
#progressive degeneration of neurons[increase HR]
NEURON_DEGENERATION = ['IGF1','NAGLU','CLU']
#Growth hormone synthesis[increase HR]
GROWTH_HORMONE_SYN = ['GHR', 'IGF1', 'IGFALS'] #'TGFB'
#Alzhimer [increase HR]
ALZHIMER_1 = ['APP','LRG1','APOA4', 'CLU', 'CTSB', 'IGFBP2','REG1A','ACE', 'NBL1', 'AGT']
ALZHIMER_2 = ['LRG1','CTSB', 'IGFBP2','REG1A', 'NBL1', 'AGT']
ALS =  ['APP', 'CTSB','MT2A', 'GRN']
#Lysosme
LYSOSOME = ['NPC2', 'CTSB', 'NAGLU']
#Progressive neurological disorder 
NEUROLOGICAL_DISORDER = [ 'GRN','MRC1','CASP1','VIM','AGT', 'LRG1','IGFBP2','CTSB','REG1A','ACE','CLU']
#Movment disorders[ increase HR]
MOVMENT_DISORDER = ['APP','CTSB', 'NPC2', 'GRN', 'CLU', 'CASP1','MT2A','CYCS']
#NEUROMUSCULAR proteins 
NEUROMUSCULAR = ["APP",
    "CTSB",
    "IGF1",
    "CTNNB1",
    "LRG1",
    "AGT",
    "NAGLU",
    "PON1",
    "VIM",
    "CASP1",
    "MT2A",
    "IGFALS",
    "CYCS",
    "FKBP1A",
    "APOA4",
    "CLU",
    "PAM",
    "PYCARD"
]

#LRG1 pathway from paper[LRG1: an emerging player in disease pathogenesis]
LRG1_PATHWAY = ['LRG1', 'AKT', 'TGFB', 'EGF', 'LRP', 'ALK5', 'IL6', 'ALK1']


#------------------------------------------------------------------------------------------------------------------
#PROTEIN INTERACTION FROM STRING 
protein_lists = {
    'LRG1': ['LRG1', 'APOA4', 'C9', 'TTR', 'CP', 'ENG', 'CYCS', 'KLK10'],
    'GSN': ['ACTA1', 'GSN', 'PFN1', 'TLN1', 'TLN2', 'CTTN', 'WASL', 'TWF2'],
    'ITIH4.1': ['CD207', 'PPP2R3C', 'DUSP5', 'CD4', 'CCR5','CCR1', 'CD209', 'SCD5', 'CLEC4M', 'ITIH4.1','C3', 'SERPING1'],
    'AGT': ['AGT', 'TAC1', 'REN', 'KNG1', 'ACE', 'GRP', 'AGTR2', 'EDNRA', 'DBKRBD2', 'ACE2','FYN', 'PFKM', 'NDUFA3', 'HSPA5', 'RIPK4','CANX', 'FAS', 'PAK3', 'TUBBB2B', 'NOMO1', 'HNRNPCL1'],
    'IGF1': ['IGF1', 'IGF1R', 'INS', 'IGFBP1', 'IGFBP3', 'IGFBP4', 'IGFBP6', 'IGFBP5', 'INSR', 'IGFALS', 'SPCS1'],
    'REG1A': ['REG1A', 'TFF1', 'CLPS', 'EXTL3', 'REG1B', 'DEFA5', 'DEFA6', 'REG3A','NEK4', 'REG1B', 'NAA20','NAA25'],
    'CTSB': ['CTSS', 'CTSB', 'CTSD', 'ANAXA2', 'CSTA', 'BCL2', 'NLRP3', 'CD74','BAG6','PROCR','MSI2','GSTP1', 'FXN', 'MYO1C','CLIC4' ],
    'IGFBP1': ['IGFBP1', 'IGF1', 'AHSG', 'AFP', 'INS', 'IGFBP4', 'IGFBP6', 'IGF1R', 'FOXO1', 'ATF4'],
    'VCAN': ['VCAN', 'BCAN', 'SELL', 'TLR2', 'PTPRZ1', 'NACT', 'BGN', 'FBN1', 'CD44', 'ACAN'],
    'FABP4': ['FABP4', 'PLIN1', 'CEBPA', 'CD36', 'ADIPOQ', 'LIPE', 'GOT2', 'FABP1', 'PPARG', 'SCARB2', 'VIM', 'ZBED1', 'CHD3','OSTF1'],
    'GHR': ['GHR', 'GH1', 'STAT5B', 'CSH2', 'CSH1', 'IGF1', 'SOCS2', 'STAT5A', 'PRL', 'JAK2'],
    'IGFBP2': ['IGFBP2', 'IGFR1', 'IGF2', 'INS', 'IGFBP3', 'MIIP', 'NIPA2', 'CYFIP1', 'TUBGCP5'],
    'MRC1': ['MRC1', 'ARG1', 'CD163', 'ITGAM', 'CD68', 'CDC45', 'TIPIN', 'WDHD1', 'PLA2R1', 'MCM10', 'MRC2'],
    'RNH1': ['RNASE4', 'ANG', 'EDDM3A', 'RNH1','KIF5B', ],
    'IGFALS': ['IGFBP4', 'IGF2', 'IGF1', 'F2', 'IGFBP5', 'IGFBP6', 'PAPPA2', 'IGFBP2', 'IGFBP1', 'IGFALS'],
    'NBL1': ['NBL1', 'BMP4', 'BMP5', 'BMP6', 'BMP8A', 'BMP7', 'GDF5', 'AMH', 'GDF7', 'BMP2', 'NCS1'],
    'EPHB4': ['EPHB4', 'ABL1', 'RASA1', 'EFNA1', 'NGEF', 'EFNA4', 'EFNB1', 'EFNA2', 'EFNB3', 'EFNA5', 'EFNB4'],
    'NPC2': ['SLC38A9', 'TMEM97', 'SREBF2', 'SREBF1', 'NPC2', 'NPC1', 'STARD3', 'GP2', 'GTPBP1', 'ABCA1'],
    'FGL1': ['FGL1', 'FGG', 'APOH', 'FGB', 'AERPINF2', 'LAG3', 'CLEC4G', 'THBS1', 'FGL2', 'FN1'],
    'PTPRJ': ['PTS', 'LCK', 'MET', 'PTPRJ', 'FLT3', 'FLT3LG', 'SDC2', 'THBS1', 'CDH5', 'CTNND1', 'CTNNB1'],
    'ADAM15': ['ADAM15', 'SRC', 'SH3GL', 'SNX9', 'MADL2L', 'GRB2', 'HCK', 'SH3PXD2A', 'NCK1', 'SH3PXD2B', 'SNX33', 'LCK', 'SH3PXD2A', 'LYN', 'FYN', 'YES1','SH3RF1','NPHP1', 'PACSIN3','ARHGEF6', 'ARHGEF7'  ],
    'APOF': ['APOF', 'APOC3', 'APOM', 'APOA4', 'CETP', 'APOA1', 'APOB', 'CLU', 'APOL1', 'APOA2', 'APOH'],
    'CILP2': ['CILP2', 'COMP', 'PSRC1', 'FURIN', 'PBX4', 'NCAN', 'TM6SF2', 'CELSR2', 'GALNT2', 'TRIB', 'ANGPLT3'],
    'SECTM1': ['SECTM1', 'BATF2', 'CD7', 'TEX19', 'ZMYND15', 'ALG6','PTPRD', 'OSMR'],
    'CHRDL1': ['CHRDL1', 'ATP1A2', 'MFAP4', 'GREM2', 'NOG', 'BMP2', 'BMP4', 'WFS1', 'SPARCL1', 'FAM20A', 'FAM20C'],
    'SPON1': ['APOE', 'SPON1', 'LRP2', 'VLDLR', 'LRP8', 'APBB1', 'APP', 'APLP2', 'BACE1', 'POFUT2', 'B3GLCT'],
    'APMAP': ['APMAP', 'EI5B', 'MTIF2', 'RPL23A', 'RPS9', 'EXOSC10', 'FDXACB1', 'IMP3', 'MRPL20', 'ECH1'],
    'FETUB': ['EEA1', 'ASTL', 'SERPINA7', 'AHSG', 'ORM2', 'FETUB', 'ORM1', 'CETP', 'MSR1', 'MEP1A', 'AGA'],
    'VSIG4': ['CSF1R', 'VSIG4', 'MS4A4A', 'CD163', 'C3AR1', 'C1QB', 'C1QC', 'C1QA', 'TMIGD3', 'C3', 'FOLR2']
}

#------------------------------------------------------------------------------------------------------------------
#KM dicotomized
KM_LIST = ['APOL1',
 'SUSD5',
 'SERPINC1',
 'MB',
 'C9',
 'LRG1',
 'F13B',
 'CELA3B',
 'LGALS3',
 'VCAM1',
 'AZGP1',
 'GRN',
 'NBL1',
 'PSMB3',
 'DSC2',
 'ITIH3',
 'APOF',
 'PAMR1',
 'IGDCC4',
 'CENPI',
 'ABHD14B',
 'SPON1',
 'TNFRSF12A']


###-----------------CORRELATION LIST---------
#proteins list from replication and discovery log before combat with correlation above 0.7
HIGHT_CORR_list_REP =['KLK13', 'PF4', 'MMRN1', 'GP6', 'C4BPA', 'LTBP1', 'CDC5L', 'TALDO1', 'TLN1', 'C18orf63', 'MPO', 'RAB11B', 'ABCA1', 'RNASE1', 'SAA1', 'SPARC', 'CDH13', 'C8B', 'PARK7', 'CD9', 'GAPDH', 'ADGRL2', 'FUCA2', 'TMSB4X', 'RPS6KC1', 'EHD1', 'ARPC3', 'S100A9', 'BLMH', 'DSG1', 'BLVRB', 'TAGLN2', 'TPM4', 'TUBA4A', 'TUBB1', 'C8A', 'JCHAIN', 'HSPA8', 'PTGDS', 'PPIA', 'PFN1', 'YWHAZ', 'DAG1', 'ARG1', 'ACTB', 'S100A4', 'PXK', 'NRGN', 'ZNF518A', 'ALAD', 'TPI1', 'LY6G6F', 'LTF', 'IGFBP4', 'PRDX6', 'ARPC2', 'PDHB', 'TREML1', 'ARPC1B', 'OXSR1', 'SLC2A1', 'IGFBP5', 'GSDMA', 'JUP', 'ROBO4', 'PARVB', 'PGLYRP1', 'CHMP4B', 'RUNDC3A', 'CTSG', 'WDR1', 'DSP', 'EIF4G1', 'CAT', 'S100A8', 'FLG2', 'CD5L', 'MMP8', 'HAGH', 'PDLIM1', 'CLEC3B', 'DEFA1', 'OSBP', 'C8G', 'CDK5RAP2', 'AGT', 'NPNT', 'ABHD14B', 'AZU1', 'ITGA6', 'CPN2', 'PNP', 'CA1', 'EIF5A', 'TPM3', 'C4A', 'CDSN', 'MYH9', 'ACTBL2', 'FERMT3', 'CASP14', 'SLC4A1', 'TUBB', 'SH3TC2', 'VASP', 'GP5', 'LIMS1', 'RSU1', 'AHSG', 'SOWAHB', 'ACTA1', 'NBL1', 'CFD', 'APLP2', 'PPBP', 'CCSER1', 'HBD', 'SERPINE1', 'YOD1', 'GCC2', 'PSMA6', 'STIP1', 'FLG', 'H3', 'SPDL1', 'ANXA2', 'ARPC5', 'BCAP31', 'CA2', 'TTN', 'WDHD1', 'HRG', 'ANKRD30B', 'ENO3', 'FLNA', 'MAP1A', 'ILK', 'HAL', 'PRDX2', 'TGM3', 'MMP9', 'CSTA', 'TJP1', 'CDH5', 'MAP3K20', 'IGFBP6', 'CPN1', 'CAP1', 'KPRP', 'SH3BGRL', 'CRP', 'WDR19', 'C4BPB', 'ALOX12', 'UBA52', 'FGB', 'F11', 'FABP5', 'CFL1', 'NEK1', 'HPSE', 'ZYX', 'PIBF1', 'BIN2', 'LGALSL', 'PGK1', 'IZUMO3', 'ACTR2', 'SAE1', 'GP9', 'PCLO', 'PROS1', 'LORICRIN', 'PSMB3', 'ENO1', 'SERPINB12', 'BPGM', 'ACTN1', 'H4C1', 'PKM', 'THBS1', 'SAA2', 'MB', 'CFHR2', 'CD59', 'FEN1', 'NME1', 'ST13', 'HSD17B4', 'RELN', 'PTK2', 'RAP1B', 'HNRNPK', 'MTPN', 'IFIT1', 'B2M', 'MSLNL', 'SPTBN5', 'TTR', 'HPRT1', 'NME2', 'RNASE7', 'CST3', 'VWF', 'PRDX1', 'TUBA1C', 'PSG5', 'LANCL1', 'SRGN', 'HBA1', 'F12', 'PDS5B', 'CFHR1', 'PAFAH1B3', 'HBB', 'FGG', 'SVEP1', 'HRNR', 'APOE', 'GP1BB', 'CTNNB1', 'H2AC18', 'TMPRSS6', 'DSC1', 'FUCA1', 'CA3', 'ITGA2B', 'SELENBP1', 'PLEK', 'LCN2', 'MCAM','APOC2','IGF2']
HIGHT_CORR_list_DIS = ['KLK13', 'GP6', 'CDC5L', 'TLN1', 'RAB11B', 'CAVIN2', 'ABCA1', 'SAA1', 'PYCARD', 'NOTCH4', 'CD9', 'FUCA2', 'TMSB4X', 'PAFAH1B2', 'RPS6KC1', 'S100A9', 'BLVRB', 'TAGLN2', 'S100A14', 'TUBA4A', 'TUBB1', 'JCHAIN', 'ABCE1', 'ARG1', 'ACTB', 'RAB27B', 'NRGN', 'ALAD', 'LY6G6F', 'COL15A1', 'UPF3B', 'TREML1', 'SLC2A1', 'OXSR1', 'GSDMA', 'BNIPL', 'PARVB', 'JUP', 'DSP', 'EIF4G1', 'S100A8', 'XP32', 'CD5L', 'PDLIM1', 'OSBP', 'CD58', 'CYCS', 'CUTA', 'ITGA6', 'KCTD12', 'CA1', 'TPSAB1', 'SERPINB8', 'EIF5A', 'ACTBL2', 'ARR3', 'FERMT3', 'CASP14', 'AGA', 'SLC4A1', 'TUBB', 'GP5', 'LIMS1', 'RSU1', 'ACTA1', 'COL1A1', 'HBD', 'CCSER1', 'GCC2', 'PSMA6', 'STIP1', 'H3', 'ANXA2', 'CAPG', 'CCDC126', 'TXN', 'CA2', 'YWHAE', 'CORO1C', 'DPM3', 'YWHAH', 'H2BC5', 'FLNA', 'PLXNB2', 'ILK', 'PKP1', 'PRDX2', 'TGM3', 'CPOX', 'CAP1', 'KPRP', 'WDR19', 'FGB', 'WDR43', 'ZYX', 'SERPINB7', 'LGALSL', 'IZUMO3', 'MAN1C1', 'GP9', 'LORICRIN', 'HSPA6', 'SPP1', 'PSMB3', 'ITGB3', 'EPB42', 'PEAR1', 'ACTN1', 'H4C1', 'TULP1', 'SAA2', 'CFHR2', 'GGCT', 'RAP1B', 'IFIT1', 'PCDH1', 'MSLNL', 'SPTBN5', 'TTR', 'RNASE7', 'LAMP1', 'TUBA1C', 'KLK7', 'HBA1', 'PDS5B', 'CTSV', 'CFHR1', 'PAFAH1B3', 'HBB', 'FGG', 'GP1BB', 'CTNNB1', 'H2AC18', 'DSC1', 'FUCA1', 'ITGA2B', 'GSS', 'TGM1', 'YWHAG', 'S100A6', 'PLEK']
SIGNIFICANT_PROTEIN_COX_UNI_LOW_CORR_DATA = ['CILP2',
 'REG1A',
 'GRN',
 'MRC1',
 'IGFBP2',
 'ADAM15',
 'EPHB4',
 'IGFALS',
 'VSIG4',
 'FGL1',
 'SPON1',
 'APMAP',
 'PON1',
 'PTPRJ',
 'IGDCC4',
 'EXTL2',
 'SHBG',
 'IGF1',
 'APOF',
 'CTSB',
 'LRG1',
 'FKBP1A',
 'CNTN3',
 'GHR',
 'SECTM1',
 'CDH1',
 'NPC2',
 'FETUB',
 'ICOSLG',
 'FRMPD3',
 'ALSFRS_slope (-unit/month)',
 'ALSFRS score (unit)']




#NO CROSS VALIDATION
SIG_PROTEIN_DIS_REP_COHORT = ['KIF5B',
 'GSN',
 'APOL1',
 'ICOSLG',
 'C1QA',
 'C9',
 'LRG1',
 'VIM',
 'FBP1',
 'CSF1',
 'LGALS3',
 'IGFBP2',
 'AZGP1',
 'ITIH3',
 'ACACA',
 'DYNC1H1',
 'NECTIN1',
 'ADAMTS13',
 'IGDCC4',
 'CENPI',
 'NRCAM',
 'ATF6B',
 'ALSFRS_slope (-unit/month)',
 '1 = Age Onset ≥ 60',
 'ALSFRS score (unit)']