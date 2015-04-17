##make basic metadata plots

load("data//ph_raw_data.RData")
##merge data by image number
library(ggplot2)

all_alp_data<-merge(cell.ftrs,image.data, by="ImageNumber")
##make ggplots cell number per unit 

ggplot(image.data, aes(x = Image_Count_Cells))+
  geom_histogram(aes(y=..density..),     
                 binwidth=1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="darkgreen")

##make ggplots cell number per unit per repaeta

ggplot(image.data,aes(factor(Image_Metadata_ArrayNumber),Image_Count_Cells ))+
  geom_boxplot(aes(fill = factor(Image_Metadata_ArrayNumber)))

ggplot(image.data, aes(Image_Count_Cells, 
      fill = as.factor(Image_Metadata_ArrayNumber))) + geom_density(alpha = 0.2) 

ggplot(image.data, aes(Image_Count_Cells, 
                       fill = as.factor(Image_Metadata_ArrayNumber))) + geom_density() 

ggplot(image.data, aes(x = Image_Count_Cells)) +
  stat_density(aes(ymax = ..density..,  ymin = -..density..),
               fill = "grey50", colour = "grey50",
               geom = "ribbon", position = "identity") +
  facet_grid(. ~ Image_Metadata_ArrayNumber) +
  coord_flip()

