import os

#MAIN_DIR = os.path.join('/home', 'labs','hornsteinlab','Collaboration','PRIMALS', 'Serum')
MAIN_DIR = os.path.join('/home', 'labs','hornsteinlab','hadarkl')
PROTEIN_DATA_DIR = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','proteins')

OUTPUTS_DIR = os.path.join('/home', 'labs', 'hornsteinlab', 'hadarkl', 'Hadar', 'primals', 'code', 'outputs', 'csv')
FIGURES_DIR = os.path.join('/home', 'labs', 'hornsteinlab', 'hadarkl', 'Hadar', 'primals', 'code', 'outputs', 'figures')
MODEL_DIR =  os.path.join('/home', 'labs', 'hornsteinlab', 'hadarkl', 'Hadar', 'primals', 'code', 'outputs', 'models')
PROTEOMIC_DISCOVERY_PROCESSED = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','proteins', 'discovery_data_processed.csv')
PROTEOMIC_DATA_REPLICATION = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','proteins', 'replication_data_processed.csv')
CLINICAL_DATA_REP = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','clinical_data', 'clinical_data_replication.csv')
CLINICAL_DATA_DISCOVERY = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','clinical_data', 'clinical_data_baseline.csv')
DELTA_ALSFRS_REP = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','clinical_data', 'delta_ALSFRS_rep.csv')
DELTA_ALSFRS_DIS = os.path.join(MAIN_DIR,'Hadar','primals', 'code', 'data', 'raw_data','clinical_data', 'delta_ALSFRS_discovery.csv')

seed = 42

# RAW_PROTEOMICS_DATA = os.path.join(MAIN_DIR, 'inputs', 'processed_data', 'proteomics', 'report_DIANN_perseus_renamed.xlsx')

# RAW_PROTEOMICS_REPLICATION_DATA = os.path.join('/home', 'labs', 'hornsteinlab', 'hadarkl', 'Hadar', 'report_ALS replication_perseus.xlsx')

PROCESSED_PROTEOMICS_DATA = os.path.join(MAIN_DIR, 'Hadar', 'als_protein_df_final.xlsx')

# CLINICAL_DATA = os.path.join(MAIN_DIR, 'Hadar', 'clinical_data_hadar.csv')
MAIN_FOLDER = os.path.join(MAIN_DIR, 'Hadar')


IMMUNE_PROTEINS = ["ORM1", "ORM2", "SERPINA1", "A2ML1", "A2M","APOA1",'APOA2', "HP", "HPR", "TFRC","TF", "ALB", "IGHA2", "IGHA1", "IGHD", "IGHG1","IGHG2",  "IGHG3", "IGHG4","IGHM", "ALB", "SERPINA3","IGKV3-7", "IGLV4-69","IGLV8-61",
                   "IGLV10-54","IGLV7-46",'IGLV4-69', 'IGLV8-61', 'IGLV4-60', 'IGLV10-54', 'IGLV5-48', 'IGLV7-46', 'IGLV5-37', 'IGLV2-18', 'IGLV3-12', 'IGLV3-10',
                   'IGLV3-9', 'IGLV4-3', 'IGLV5-45', 'IGLV5-52', 'IGLV1-36', 'IGLV9-49', 'IGLV5-39', 'IGLV1-47', 'IGLV1-51', 'IGLV1-40', 'IGLV2-14', 'IGLV2-23', 'IGLV2-11', 'IGLV2-8', 'IGLV3-19', 'IGLV3-1', 'IGLV3-25', 'IGLV3-27',
                   'IGLV6-57', 'IGLV3-21','IGHV3-64', 'IGHV4-4', 'IGHV3OR16-12', 'IGHV1OR15-1', 'IGHV3OR15-7', 'IGHV4-30-2', 'IGHV1OR21-1', 'IGHV1-45', 'IGHV3-49', 'IGHV6-1', 'IGHV3-15', 'IGHV2-26', 'IGHV3-73', 'IGHV3-74', 'IGHV3-43',
                   'IGHV3-72', 'IGHV1-3', 'IGHV1-18', 'IGHV3-20', 'IGHV1-24', 'IGHV4-28', 'IGHV3-35', 'IGHV3-38', 'IGHV5-51', 'IGHV2-70D',
                   'IGHA2', 'IGHV3-64D', 'IGHV5-10-1', 'IGHV1-69', 'IGHV1-46', 'IGHV3-11', 'IGHV3-48', 'IGHV3-23', 'IGHV3-13',
                   'IGHV3-53', 'IGHV3-30;IGHV3-30-5', 'IGHV3-7', 'IGHV3-9', 'IGHV2-70', 'IGHV2-5', 'IGHE', 'IGHG1',
                   'IGHG2', 'IGHG3', 'IGHG4', 'IGHM', 'IGHA1', 'IGHD', 'IGHV4-34', 'IGHV1-8', 'IGHV3-38-3', 'IGHV1-2','IGLV4-69',
                   'IGLV8-61', 'IGLV4-60', 'IGLV10-54', 'IGLV5-48', 'IGLV7-46', 'IGLV5-37', 'IGLV2-18', 'IGLV3-12', 'IGLV3-10', 'IGLV3-9', 'IGLV4-3', 'IGLV5-45', 'IGLV5-52', 'IGLV1-36', 'IGLV9-49', 'IGLV5-39', 'IGLV1-47', 'IGLV1-51', 'IGLV1-40',
                   'IGLV2-14', 'IGLV2-23', 'IGLV2-11', 'IGLV2-8', 'IGLV3-19', 'IGLV3-1', 'IGLV3-25', 'IGLV3-27', 'IGLV6-57', 'IGLV3-21','IGHV3-72.1','IGHA2.1',
                   'IGKV3-7', 'IGKV1-27', 'IGKV2D-30', 'IGKV3D-15', 'IGKV1D-8', 'IGKV2-40;IGKV2D-40', 'IGKV6D-21', 'IGKV1-13;IGKV1D-13', 'IGKV6-21', 'IGKV3D-20', 'IGKV1-8', 'IGKV2-24', 'IGKV1-12;IGKV1D-12',
                   'IGKJ4', 'IGKJ1', 'IGKJ3', 'IGLC7', 'IGKV2-29', 'IGLL5', 'IGKV1-33;IGKV1D-33', 'IGKV1-39;IGKV1D-39', 'IGKV1-17', 'IGKV1D-16', 'IGKV1-5', 'IGKV3-20', 'IGKC', 'IGKV1-16', 'IGKV3-11', 'IGKV2-30', 'IGKV4-1', 'IGLC3', 'IGLL1', "IGLV3-16"]


KERATIN = ['KRTAP4-9', 'KRT10', 'KRT87P', 'KRTAP16-1', 'KRT19', 'KRT86', 'KRT33A',  'KRT34', 'KRT36', 'KRT38', 'KRT75', 'KRT14', 'KRT6A', 'KRT6B', 'KRT1', 'KRT18', 'KRT8', 'KRT19', 'KRT7', 'KRT16', 'KRTAP2-3;KRTAP2-4', 'KRT3',  'KRT10', 'KRT13', 'KRT5', 'KRT15', 'KRT4', 'KRT9', 'KRT2', 'KRT6C', 'KRTAP10-10', 'KRTAP10-3', 'KRTAP10-9', 'KRTAP10-11', 'KRTAP10-12',
       'KRTDAP', 'KRT83', 'KRT85', 'KRT76', 'KRT17', 'KRT33B', 'KRT32',
       'KRT81', 'KRT72', 'KRT31', 'KRTAP20-2', 'KRTAP20-1', 'KRTAP6-1',
       'KRTAP19-5', 'KRTAP13-4', 'KRTAP24-1', 'KRT71', 'KRTAP19-7',
       'KRTAP13-2', 'KRT79', 'KRT40', 'KRT39', 'KRT80', 'KRT28',
       'KRT27', 'KRT25', 'KRT77', 'KRTAP13-1', 'KRTAP11-1', 'KRTAP8-1',
       'KRTAP7-1', 'KRTAP1-3', 'KRT78', 'KRT35', 'KRT12', 'KRTAP9-8',
       'KRTAP9-3', 'KRTAP4-6', 'KRTAP4-11', 'KRTAP4-1', 'KRTAP4-9', 'KRTAP4-8',
       'KRTAP4-5', 'KRTAP4-4', 'KRTAP4-3', 'KRTAP4-2', 'KRTAP3-3', 'KRTAP3-2',
       'KRTAP3-1', 'KRTAP1-5', 'KRTAP2-1;KRTAP2-2', 'KRT23', 'KRT84', 'KRT82']



REPLICATION_INDEX = ['705', '712', '717', '735', '746', '766', '789', '797', '815', '816', '828',
                   '844', '850', '852', '853', '859', '866', '869', '870', '873', '874', '879', '885',
                   '889', '897', '901', '904', '905', '915', '928', '931', '934', '935', '942', '951', '953',
                   '956', '959', '963', '985', '1011', '1021', '1024', '1025', '1033', '1036', '1043',
                   '1053', '1058', '1066', '1072', '1074', '1078', '1079', '1080', '1084', '1086', '1091',
                   '1094', '1096', '1101', '1112', '1113', '1123', '1131', '1144', '1145', '1146', '1163', '1166',
                   '1170', '1211', '1217', '1228', '1230', '1235', '1236', '1238', '1239', '1240', '1242', '1246',
                   '1249', '1250', '1254', '1257', '1271', '1278', '1280', '1283', '1284', '1304', '1314', '1317',
                   '1319', '1320', '1321', '1323', '1324', '1328', '1336', '1354', '1356', '1359', '1360', '1374',
                   '1387', '1392', '1394', '1408', '1421', '1427', '1431', '1435', '1440', '1441', '1444', '1463',
                   '1471', '1475', '1476', '1477', '1482', '1493', '1494', '1500', '1501', '1505', '1510', '1519',
                   '1523', '1527', '1535', '1551', '1554', '1556', '1561', '1573', '1576', '1581', '1587', '1590',
                   '1600', '1606', '1618', '1630', '1636', '1657', '1686', '1695', '1714', '1726', '1727', '1728', '1729',
                   '1733', '1734', '1736', '1740', '1741', '1745', '1765', '1771', '1776', '1777', '1785', '1789', '1790',
                   '1792', '1797', '1800', '1803', '1810', '1832', '1834', '1863', '1892', '1899', '1900', '1912', '1913', '1946']
