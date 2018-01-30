rm(list=ls())
setwd("/Users/anju/Desktop/UEGP_stats")
uprural<-read.table("Metaphlan_species_models_UpRural.txt", header = T, sep ="\t")
uprural<-uprural[uprural$pValuesUpRuralAdjusted<0.05 , ]
uprural_shape<-cbind(uprural$names, rep("ring_shape", length(uprural$names)), rep("1", length(uprural$names)), rep("v", length(uprural$names)))
uprural_color<-cbind(uprural$names, rep("ring_color", length(uprural$names)), rep("1", length(uprural$names)), rep("b", length(uprural$names)))


infpce<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_PceInf.txt", header = T, sep ="\t")
infpce<-infpce[infpce$pValuesPceInfAdjusted<0.05 |infpce$pValuesLocationAdjusted<0.05 | infpce$pValuesTimepointAdjusted<0.05 ,]
infpce$bugs<-infpce$names
infpce_shape<-data.frame(clade = infpce$bugs, ring_option=rep("ring_shape", length(infpce$bugs)), ring_level=rep("1", length(infpce$bugs)), parameter=rep("v", length(infpce$bugs)))
infpce_color<-data.frame(clade = infpce$bugs, ring_option=rep("ring_color", length(infpce$bugs)), ring_level=rep("1", length(infpce$bugs)), parameter=rep("#0000FF", length(infpce$bugs)))
#infpce_color<-cbind(infpce$names[1:202], rep("ring_color", length(infpce$names)), rep("1", length(infpce$names)), rep("#FF0000", length(infpce$names)))

pceate<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_PceAte.txt", header = T, sep ="\t")
pceate<-pceate[pceate$pValuesPceAteAdjusted<0.05 |pceate$pValuesLocationAdjusted<0.05 | pceate$pValuesTimepointAdjusted<0.05 ,]
pceate$bugs<-pceate$names
pceate_shape<-data.frame(clade = pceate$bugs, ring_option=rep("ring_shape", length(pceate$bugs)), ring_level=rep("2", length(pceate$bugs)), parameter=rep("v", length(pceate$bugs)))
pceate_color<-data.frame(clade = pceate$bugs, ring_option=rep("ring_color", length(pceate$bugs)), ring_level=rep("2", length(pceate$bugs)), parameter=rep("#FFA500", length(pceate$bugs)))

atefce<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_FceAte.txt", header = T, sep ="\t")
atefce<-atefce[atefce$pValuesFceAteAdjusted<0.05 |atefce$pValuesLocationAdjusted<0.05 | atefce$pValuesTimepointAdjusted<0.05 ,]
atefce$bugs<-atefce$names
atefce_shape<-data.frame(clade = atefce$bugs, ring_option=rep("ring_shape", length(atefce$bugs)), ring_level=rep("3", length(atefce$bugs)), parameter=rep("v", length(atefce$bugs)))
atefce_color<-data.frame(clade = atefce$bugs, ring_option=rep("ring_color", length(atefce$bugs)), ring_level=rep("3", length(atefce$bugs)), parameter=rep("#FF0000", length(atefce$bugs)))

inffce<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_FceInf.txt", header = T, sep ="\t")
inffce<-inffce[inffce$pValuesFceInfAdjusted<0.05 |inffce$pValuesLocationAdjusted<0.05 | inffce$pValuesTimepointAdjusted<0.05 ,]
inffce$bugs<-inffce$names
inffce_shape<-data.frame(clade = inffce$bugs, ring_option=rep("ring_shape", length(inffce$bugs)), ring_level=rep("4", length(inffce$bugs)), parameter=rep("v", length(inffce$bugs)))
inffce_color<-data.frame(clade = inffce$bugs, ring_option=rep("ring_color", length(inffce$bugs)), ring_level=rep("4", length(inffce$bugs)), parameter=rep("#800000", length(inffce$bugs)))

fcedsa<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_FCEDowna.txt", header = T, sep ="\t")
fcedsa<-fcedsa[fcedsa$pValuesFCEDownAdjusted<0.05 |fcedsa$pValuesLocationAdjusted<0.05 | fcedsa$pValuesTimepointAdjusted<0.05 ,]
fcedsa$bugs<-fcedsa$names
fcedsa_shape<-data.frame(clade = fcedsa$bugs, ring_option=rep("ring_shape", length(fcedsa$bugs)), ring_level=rep("5", length(fcedsa$bugs)), parameter=rep("v", length(fcedsa$bugs)))
fcedsa_color<-data.frame(clade = fcedsa$bugs, ring_option=rep("ring_color", length(fcedsa$bugs)), ring_level=rep("5", length(fcedsa$bugs)), parameter=rep("#006400", length(fcedsa$bugs)))

fcedsb<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_FCEDownb.txt", header = T, sep ="\t")
fcedsb<-fcedsb[fcedsb$pValuesFCEDownAdjusted<0.05 |fcedsb$pValuesLocationAdjusted<0.05 | fcedsb$pValuesTimepointAdjusted<0.05 ,]
fcedsb$bugs<-fcedsb$names
fcedsb_shape<-data.frame(clade = fcedsb$bugs, ring_option=rep("ring_shape", length(fcedsb$bugs)), ring_level=rep("6", length(fcedsb$bugs)), parameter=rep("v", length(fcedsb$bugs)))
fcedsb_color<-data.frame(clade = fcedsb$bugs, ring_option=rep("ring_color", length(fcedsb$bugs)), ring_level=rep("6", length(fcedsb$bugs)), parameter=rep("#800080", length(fcedsb$bugs)))

dsadsb<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_dsAdsB.txt", header = T, sep ="\t")
dsadsb<-dsadsb[dsadsb$pValuesdsAdsBAdjusted<0.05 |dsadsb$pValuesLocationAdjusted<0.05 | dsadsb$pValuesTimepointAdjusted<0.05 ,]
dsadsb$bugs<-dsadsb$names
dsadsb_shape<-data.frame(clade = dsadsb$bugs, ring_option=rep("ring_shape", length(dsadsb$bugs)), ring_level=rep("7", length(dsadsb$bugs)), parameter=rep("v", length(dsadsb$bugs)))
dsadsb_color<-data.frame(clade = dsadsb$bugs, ring_option=rep("ring_color", length(dsadsb$bugs)), ring_level=rep("7", length(dsadsb$bugs)), parameter=rep("#696969", length(dsadsb$bugs)))


fceds<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_FCEDown.txt", header = T, sep ="\t")
fceds<-fceds[fceds$pValuesFCEDownAdjusted<0.05 |fceds$pValuesLocationAdjusted<0.05 | fceds$pValuesTimepointAdjusted<0.05 ,]
fceds$bugs<-fceds$names
fceds_shape<-data.frame(clade = fceds$bugs, ring_option=rep("ring_shape", length(fceds$bugs)), ring_level=rep("8", length(fceds$bugs)), parameter=rep("v", length(fceds$bugs)))
fceds_color<-data.frame(clade = fceds$bugs, ring_option=rep("ring_color", length(fceds$bugs)), ring_level=rep("8", length(fceds$bugs)), parameter=rep("#DDA6FF", length(fceds$bugs)))

updown<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_UpDown.txt", header = T, sep ="\t")
updown<-updown[updown$pValuesUpDownAdjusted<0.05 |updown$pValuesLocationAdjusted<0.05 | updown$pValuesTimepointAdjusted<0.05 ,]
updown$bugs<-updown$names
updown_shape<-data.frame(clade = updown$bugs, ring_option=rep("ring_shape", length(updown$bugs)), ring_level=rep("9", length(updown$bugs)), parameter=rep("v", length(updown$bugs)))
updown_color<-data.frame(clade = updown$bugs, ring_option=rep("ring_color", length(updown$bugs)), ring_level=rep("9", length(updown$bugs)), parameter=rep("#6ADB48", length(updown$bugs)))

hospres<-read.table("/Users/anju/Desktop/UEGP_stats/all_old_anal/metaphlan_species_filtered/Metaphlan_species_models_HospRes.txt", header = T, sep ="\t")
hospres<-hospres[hospres$pValuesHospResAdjusted<0.05 |hospres$pValuesLocationAdjusted<0.05 | hospres$pValuesTimepointAdjusted<0.05 ,]
hospres$bugs<-hospres$names
hospres_shape<-data.frame(clade = hospres$bugs, ring_option=rep("ring_shape", length(hospres$bugs)), ring_level=rep("10", length(hospres$bugs)), parameter=rep("v", length(hospres$bugs)))
hospres_color<-data.frame(clade = hospres$bugs, ring_option=rep("ring_color", length(hospres$bugs)), ring_level=rep("10", length(hospres$bugs)), parameter=rep("#203E5F", length(hospres$bugs)))

rings<-rbind(infpce_shape, pceate_shape, atefce_shape, inffce_shape, fcedsa_shape,fcedsb_shape,dsadsb_shape,fceds_shape,updown_shape,hospres_shape, infpce_color, pceate_color, atefce_color, inffce_color, fcedsa_color,fcedsb_color,dsadsb_color,fceds_color,updown_color,hospres_color)

write.table(rings, file= "pvalues_rings_colors.txt", row.names = F, sep ="\t")






