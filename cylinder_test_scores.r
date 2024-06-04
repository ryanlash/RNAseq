cylinder<-read.table('grouped_cylinder_scores.csv', header=T, sep=',')

library(ggplot2)
head(cylinder)
ct<-ggplot(cylinder, aes(x=Treatment,y=Percent_Contralateral.Unilateral)) +
  geom_bar(stat="identity", aes(fill=Treatment))+
  labs(title = "Cylinder Test Results")
ct+labs(x='Treatment Group', y='% Contralateral/Unilateral',title='Unilateral Limb Use')
ct2<-ggplot(cylinder, aes(x=Treatment,y=Percent_Bilateral_Touches)) +
  geom_bar(stat="identity", aes(fill=Treatment))+
  labs(title = "Cylinder Test Results")
ct2+labs(x='Treatment Group', y='% Bilateral/Total', title='Bilateral Limb Use')
library(dplyr)
summarise(cylinder, )