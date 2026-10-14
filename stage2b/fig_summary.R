suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
oi<-c(pair="#CC79A7",B="#0072B2",F="#009E73")
th<-theme_classic(base_size=8)+theme(plot.tag=element_text(face="bold",size=11),legend.position="none",plot.title=element_text(size=8,hjust=.5))
a<-read.csv("../stage2a/s2a_threshold_vs_aB_v4.csv")|>filter(mode=="pair")
pa<-ggplot(a,aes(aB,threshold))+geom_vline(xintercept=c(3.59,9.97),linetype="dashed",colour="grey55",linewidth=.3)+
  geom_line(colour=oi["pair"],linewidth=.7,na.rm=TRUE)+geom_point(colour=oi["pair"],size=1.1,na.rm=TRUE)+
  annotate("text",x=10.5,y=0.75,label="critical values\n(3.6 and 10.0)",size=2.1,hjust=0,colour="grey35")+
  annotate("text",x=2,y=0.05,label="pathogen\npersists",size=2.1,colour="#A04000")+
  scale_y_continuous(limits=c(0,0.9))+labs(x="Cross-feeding strength (acetate released by B)",y="Pair threshold dose (per species)",title="Cross-feeding must exceed a critical value")+th
b<-read.csv("../stage2a/s2a_threshold_vs_P0_v4.csv")|>filter(mode=="pair")
pb<-ggplot(b,aes(P0,threshold))+geom_line(colour=oi["pair"],linewidth=.7)+geom_point(colour=oi["pair"],size=1.1)+
  scale_x_log10()+labs(x="Pathogen inoculum",y="Pair threshold dose (per species)",title="More pathogen needs more dose")+th
c<-read.csv("../stage2a/s2a_dose_heterogeneity_v4.csv")|>mutate(P0=factor(P0,labels=c("Low pathogen","10x pathogen")))
pc<-ggplot(c,aes(dose,frac_protected,colour=P0))+geom_line(linewidth=.7)+geom_point(size=1.1)+scale_x_log10()+
  scale_colour_manual(values=c("#0072B2","#D55E00"))+labs(x="Dose per species",y="Fraction of virtual hosts protected",title="Host variation gives graded dose-response")+
  th+theme(legend.position=c(.75,.25),legend.title=element_blank(),legend.key.size=unit(.3,"cm"),legend.text=element_text(size=6))
m<-read.csv("s2b_matched_dose.csv")|>mutate(grp=case_when(Bfrac==0~"F alone",Bfrac==1~"B alone",TRUE~"Mix"))
pd<-ggplot(m,aes(factor(total),factor(Bfrac),fill=tolower(as.character(win_spatial))))+geom_tile(colour="white")+scale_fill_manual(values=c(`false`="#F3D9CC",`true`="#CDE9DF"))+labs(x="Total inoculum (B + F)",y="B fraction (0 = F alone, 1 = B alone)",title="Singles fail, mixtures succeed (spatial)")+th
fig<-(pa|pb)/(pc|pd)+plot_annotation(tag_levels="a")
ggsave("fig_summary_main.png",fig,width=7.2,height=5.4,dpi=300);ggsave("fig_summary_main.pdf",fig,width=7.2,height=5.4)
