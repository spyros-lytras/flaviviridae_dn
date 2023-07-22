library(PhylogeneticEM)
library(phytools)

#https://pbastide.github.io/PhylogeneticEM/articles/tutorial.html


setwd('')




############# FN FULL TREE #############


#full tree

mytree <- ape::read.tree('Flaviset_160323_aa_RdRp_10_rt.nwk')
plot(mytree, show.tip.label = TRUE)


#choose most likely model
#https://doi.org/10.1016/j.ympev.2013.02.008
#chronomodels <- c('relaxed', 'correlated', 'clock')

ultratree =chronos(mytree, lambda=1, model = 'relaxed')
#relaxed: log-Lik = -283.1224 PHIIC = 2663.58 --> BEST
#correlated: log-Lik = -297.1143, PHIIC = 2690.59
#clock: log-Lik = -289.2297, PHIIC = 1276.46 

dns <- c('UpA','CpG')

for (d in dns) {
  
  mydata <- read.csv(paste('phyloEM//Flaviset_160323_10_',d, '_RSDUc.csv', sep=''), row.names = 1, header= TRUE)
  
  mydatamat = data.matrix(mydata)
  
  plot(ultratree, show.tip.label = TRUE)
  
  myres <- PhyloEM(phylo = ultratree,
                   Y_data = mydatamat,
                   process = "scOU",                   ## scalar OU model
                   random.root = TRUE,                 ## Root is stationary (true model)
                   stationary.root = TRUE,
                   alpha = tail(find_grid_alpha(ultratree), 5), #get the 5 highest alpha values of the vector
                   K_max = 10,                         ## Maximal number of shifts
                   parallel_alpha = TRUE,              
                   Ncores = 5)
  
  myres
  
  png(paste('phyloEM//Flaviset_160323_10_', d, '_RSDUc_phyloEM_all_set1.png', sep=''), width = 2000, height = 2000)
  treeplot <- plot(myres)
  dev.off() 
  
  png(paste('phyloEM//Flaviset_160323_10_', d, '_RSDUc_crit_all_set1.png', sep=''), width = 500, height = 500)
  critplot <- plot_criterion(myres)
  dev.off() 

}


#####################################################



####  do subtrees  ####


infer_shifts <- function(subtips, rsducdat){
  
  subtree <- drop.tip(mytree, tip = setdiff(mytree$tip.label, subtips))
  
  #http://dx.doi.org/10.1016/j.ympev.2013.02.008
  ultratree <- chronos(subtree, lambda=1, model='relaxed')
  
  
  subdata <- subset(rsducdat, select = ultratree$tip.label)
  subdata <- subdata[, ultratree$tip.label]
  subdata <- data.matrix(subdata)
  
  
  thisrez <- PhyloEM(phylo = ultratree,
                     Y_data = subdata,
                     process = "scOU",                   ## scalar OU model
                     random.root = TRUE,                 ## Root is stationary (true model)
                     stationary.root = TRUE,
                     alpha = tail(find_grid_alpha(ultratree), 5),
                     K_max = 10,                         ## Maximal number of shifts
                     parallel_alpha = TRUE,              
                     Ncores = 5)
  
  return(thisrez)
}




cpgdata <- read.csv(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc.csv', sep=''), row.names = 1, header= TRUE)
upadata <- read.csv(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc.csv', sep=''), row.names = 1, header= TRUE)


####  flaviviruses  ####

flavi = c('NC_017086', 'NC_024806', 'NC_024805', 'NC_003675', 'NC_012932', 'KY320649', 'NC_005064', 'M91671', 'NC_001564', 'NC_034204', 'NC_034017', 'LC567152', 'MN567479', 'JQ268258', 'NC_027817', 'MZ358857', 'JF926699', 'NC_001477', 'NC_001475', 'NC_001474', 'NC_002640', 'NC_005039', 'LC567153', 'NC_030290', 'NC_027819', 'KX669682', 'MH899446', 'NC_012671', 'NC_021069', 'HE574574', 'NC_008604', 'KC505248', 'NC_030400', 'NC_024299', 'NC_035118', 'NC_031327', 'NC_035187', 'LC540441', 'MZ209680', 'KM225263', 'KM225265', 'NC_001437', 'NC_006551', 'NC_001672', 'NC_034007', 'KJ469370', 'NC_004119', 'JX477686', 'NC_001563', 'KC496020', 'NC_026620', 'NC_003635', 'NC_003676', 'NC_034151', 'NC_009028', 'AY632538', 'NC_003690', 'AY632542', 'NC_002031', 'NC_027999', 'DQ462443', 'NC_012534', 'NC_027709', 'NC_008719', 'NC_012735', 'D00246', 'EU082200', 'NC_023424', 'NC_033726', 'DQ859056', 'NC_033698', 'NC_001809', 'DQ235151', 'NC_005062', 'JF416960', 'AF331718', 'NC_003687', 'DQ235149', 'NC_033723', 'DQ859067', 'NC_033697', 'NC_033699', 'NC_033693', 'NC_033725', 'NC_030289', 'NC_040610', 'KF917538', 'MK908103', 'NC_009029', 'NC_033724', 'NC_012532', 'NC_029055', 'NC_015843', 'KM225264', 'EU082199', 'NC_040788', 'MW032264', 'NC_008718', 'NC_026624', 'EU159426', 'NC_024017', 'MT762108', 'AY898809', 'NC_000943', 'OX394156', 'OX394161', 'MK473878', 'MK473877', 'MK473881', 'OX394137', 'OX394149', 'NC_009026', 'NC_012533', 'KC734550', 'JF312912', 'NC_018705', 'NC_007580', 'NC_026623', 'NC_032088', 'KY347801', 'NC_016997', 'KY290249', 'NC_023439', 'NC_033721', 'LC582740')

flavi_cpgrez <- infer_shifts(flavi, cpgdata)
plot_criterion(flavi_cpgrez)
plot(flavi_cpgrez)

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_phyloEM_flavi_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(flavi_cpgrez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_crit_flavi_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(flavi_cpgrez)
dev.off() 


flavi_uparez <- infer_shifts(flavi, upadata)
plot_criterion(flavi_uparez)
plot(flavi_uparez)

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_phyloEM_flavi_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(flavi_uparez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_crit_flavi_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(flavi_uparez)
dev.off() 



####  jingmenviruses  ####

jigmen = c('MG703253', 'NC_024113', 'MK673133', 'NC_034222', 'MH158415', 'LC628180', 'MZ244283', 'MW556730', 'MW314690', 'NC_028398', 'MW208795', 'KP714089', 'MW208799', 'NC_028396', 'NC_028382', 'MW033628', 'KR902721', 'MW023854', 'MZ771214', 'LC505052', 'KM521552', 'NC_028400', 'MW314684', 'MW208798', 'MN551116', 'MW033625', 'MW023851', 'MN764158', 'KX883002', 'MW208796', 'MW208797', 'MN558700', 'MW896893', 'MW896920', 'MW208755')


jigmen_cpgrez <- infer_shifts(jigmen, cpgdata)
plot_criterion(jigmen_cpgrez)
plot(jigmen_cpgrez)

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_phyloEM_jigmen_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(jigmen_cpgrez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_crit_jigmen_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(jigmen_cpgrez)
dev.off() 


jigmen_uparez <- infer_shifts(jigmen, upadata)
plot_criterion(jigmen_uparez)
plot(jigmen_uparez)

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_phyloEM_jigmen_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(jigmen_uparez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_crit_jigmen_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(jigmen_uparez)
dev.off() 



####  long genome flaviviruses  ####

lgf = c('KU754513', 'MH620810', 'MZ852365', 'MW194893', 'MZ209734', 'NC_035071', 'NC_020252', 'MW314681', 'MN714664', 'NC_024077', 'NC_028367', 'KR902737', 'KR902735', 'MW208764', 'MW208762', 'MW208763', 'OK491478', 'MW208765', 'MW208758', 'MW314680', 'NC_028375', 'MW434111', 'MW208757', 'NC_028373', 'MW314679', 'NC_028370', 'OK491477', 'MW208756', 'MW208761', 'LC516844', 'NC_028137', 'MH778148', 'NC_028374', 'MZ209983', 'MZ210000', 'KM405246', 'MW208767', 'MW314682', 'MW561135', 'KR902739', 'KR902730')


lgf_cpgrez <- infer_shifts(lgf, cpgdata)
plot_criterion(lgf_cpgrez)
plot(lgf_cpgrez)

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_phyloEM_lgf_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(lgf_cpgrez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_crit_lgf_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(lgf_cpgrez)
dev.off() 


lgf_uparez <- infer_shifts(lgf, upadata)
plot_criterion(lgf_uparez)
plot(lgf_uparez)

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_phyloEM_lgf_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(lgf_uparez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_crit_lgf_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(lgf_uparez)
dev.off() 




####  pestiviruses  ####

pesti = c('OU592965', 'MH282908', 'NC_038964', 'OM030320', 'NC_035432', 'MK910227', 'AF144618', 'NC_003679', 'NC_025677', 'NC_003678', 'KJ950914', 'OX394172', 'OX394184', 'MG599985', 'NC_001461', 'MZ664273', 'NC_023176', 'OM030319', 'MG599982', 'NC_002657', 'NC_018713', 'MH231127', 'NC_012812', 'NC_024018', 'OM451132', 'MK636874', 'OM451131', 'OM451127', 'KY370099', 'KY370100', 'OM480523', 'MG599984', 'OX394182', 'OX394178')

pesti_cpgrez <- infer_shifts(pesti, cpgdata)
plot_criterion(pesti_cpgrez)
plot(pesti_cpgrez)

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_phyloEM_pesti_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(pesti_cpgrez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_crit_pesti_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(pesti_cpgrez)
dev.off() 


pesti_uparez <- infer_shifts(pesti, upadata)
plot_criterion(pesti_uparez)
plot(pesti_uparez)

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_phyloEM_pesti_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(pesti_uparez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_crit_pesti_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(pesti_uparez)
dev.off() 



####  hepaci-pegiviruses  ####

hepacipegi = c('KY370092', 'MH370348', 'NC_038428', 'OX394165', 'MT210618', 'NC_025679', 'NC_021154', 'NC_038435', 'NC_038434', 'KT439329', 'NC_027998', 'KC796073', 'NC_030291', 'KC796088', 'MK059751', 'KC145265', 'NC_020902', 'MT513216', 'NC_001837', 'NC_001710', 'NC_024377', 'KU351670', 'NC_034442', 'MW897329', 'MW365447', 'MW897327', 'MW091548', 'MT210610', 'NC_031916', 'OM030324', 'MH844500', 'NC_040815', 'NC_038432', 'MG599989', 'MG599988', 'MG599987', 'OX394158', 'MG600000', 'OX394189', 'NC_031947', 'NC_038431', 'NC_038430', 'OX394076', 'OX394143', 'OM451124', 'NC_025673', 'MH824541', 'OX394154', 'OX394181', 'OX394160', 'OX394170', 'OX394167', 'OX394168', 'OX394163', 'OX394180', 'OX394173', 'OX394171', 'MG599997', 'MG599996', 'MG599995', 'MG599994', 'NC_028377', 'MG599991', 'MG599992', 'OX394187', 'MZ545672', 'OM030326', 'MG599999', 'NC_038429', 'NC_026797', 'MN062427', 'NC_021153', 'NC_038427', 'OM030327', 'MH028007', 'JF744991', 'NC_024889', 'NC_004102', 'OX394179', 'MN133813', 'OX394169', 'MK737641', 'MT210605', 'MT210612', 'NC_025672', 'OX394177', 'MZ393518', 'MN635449', 'NC_001655', 'NC_038426', 'NC_031950', 'OX394150', 'MG599993', 'MG334001')


hepacipegi_cpgrez <- infer_shifts(hepacipegi, cpgdata)
plot_criterion(hepacipegi_cpgrez)
plot(hepacipegi_cpgrez)

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_phyloEM_hepacipegi_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(hepacipegi_cpgrez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_CpG_RSDUc_crit_hepacipegi_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(hepacipegi_cpgrez)
dev.off() 


hepacipegi_uparez <- infer_shifts(hepacipegi, upadata)
plot_criterion(hepacipegi_uparez)
plot(hepacipegi_uparez)

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_phyloEM_hepacipegi_set1.png', sep=''), width = 1000, height = 1000)
treeplot <- plot(hepacipegi_uparez)
dev.off() 

png(paste('phyloEM//Flaviset_160323_10_UpA_RSDUc_crit_hepacipegi_set1.png', sep=''), width = 500, height = 500)
critplot <- plot_criterion(hepacipegi_uparez)
dev.off() 







