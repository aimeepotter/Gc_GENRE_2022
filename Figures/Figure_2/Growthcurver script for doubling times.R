# First, load the package. 
library(growthcurver)
library(writexl)

# Load the sample growth curve data provided in the Growthcurver package.
# The first column is the time in hours, and there is one column 
# for each well in a 96-well plate.
growthdataCFUsGCBL<-read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/growthdataCFUsGCBL.txt", header = TRUE)
growthdataCFUsMDM<-read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/growthdataCFUsMDM.txt", header = TRUE)
growthdataCFUsRPMI<-read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/growthdataCFUsRPMI.txt", header = TRUE)
growthdataODsGCBL<- read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/growthdataODsGCBL.txt", header = TRUE)
growthdataODsMDM<- read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/growthdataODsMDM.txt", header = TRUE)


a<- growthdataCFUsGCBL
b<- growthdataCFUsMDM
c<- growthdataCFUsRPMI
d<- growthdataODsGCBL
e<- growthdataODsMDM

gc_outCFUsGCBL <- SummarizeGrowthByPlate(a, t_trim=300,bg_correct= "none", plot_fit = TRUE, 
                                 plot_file = "gc_plotsCFUsGCBL.pdf")
gc_outCFUsMDM <- SummarizeGrowthByPlate(b,bg_correct= "none", plot_fit = TRUE, 
                                         plot_file = "gc_plotsCFUsMDM.pdf")
gc_outCFUsRPMI <- SummarizeGrowthByPlate(c, bg_correct= "none", plot_fit = TRUE, 
                                         plot_file = "gc_plotsCFUsRPMI.pdf")


gc_outODsGCBL <- SummarizeGrowthByPlate(d,  t_trim=300, bg_correct= "none", plot_fit = TRUE, 
                                 plot_file = "gc_plotsODsGCBL.pdf")
gc_outODsMDM <- SummarizeGrowthByPlate(e, bg_correct= "none", plot_fit = TRUE, 
                                        plot_file = "gc_plotsODsMDM.pdf")


library(qpdf)
pdf_combine(c("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/gc_plotsCFUsGCBL.pdf", "~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/gc_plotsCFUsMDM.pdf", "~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/gc_plotsCFUsRPMI.pdf"),
                  output = "growthcurver_output_CFUs.pdf")
pdf_combine(c("~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/gc_plotsODsGCBL.pdf", "~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/gc_plotsODsMDM.pdf" ),
            output = "growthcurver_output_ODS.pdf")

write_xlsx(gc_outCFUsGCBL ,"~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/GCBL_CFU_DoublingTimes.xlsx")

write_xlsx(gc_outCFUsMDM ,"~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/MDM_CFU_DoublingTimes.xlsx")
write_xlsx(gc_outCFUsRPMI ,"~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/RPMI_CFU_DoublingTimes.xlsx")
write_xlsx(gc_outODsGCBL ,"~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/GCBL_OD_DoublingTimes.xlsx")

write_xlsx(gc_outODsMDM ,"~/UVA/Metabolic_Modeling/organized/Figures/Fig_2/MDM_OD_DoublingTimes.xlsx")


