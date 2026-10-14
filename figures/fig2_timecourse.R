suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
d<-read.csv("../stage2a/s2a_timecourse_v4.csv")|>mutate(case=factor(case,c("B alone","F alone","B + F pair"),c("B alone","F alone","B + F pair")))
oi<-c(B="#0072B2",F="#009E73",P="#D55E00")
th<-theme_classic(base_size=8)+theme(plot.tag=element_text(face="bold",size=11),plot.title=element_text(size=8,hjust=.5),legend.key.size=unit(.3,"cm"),legend.text=element_text(size=6.5),legend.title=element_blank())
mk<-function(cs,leg=FALSE){
  x<-d|>filter(case==cs,variable%in%c("B","F","P"))
  p<-ggplot(x,aes(time,value,colour=variable,linetype=variable))+geom_line(linewidth=.7)+
    scale_colour_manual(values=oi,labels=c(B="B (Bifidobacterium)",F="F (Faecalibacterium)",P="P (pathogen)"))+
    scale_linetype_manual(values=c(B="22",F="solid",P="solid"),labels=c(B="B (Bifidobacterium)",F="F (Faecalibacterium)",P="P (pathogen)"))+
    scale_x_continuous(breaks=c(0,100,200,300),expand=expansion(mult=c(0,.04)))+coord_cartesian(ylim=c(0,1.15))+
    labs(title=cs,x="Time (h)",y=if(leg)"Biomass (a.u.)" else NULL)+th
  if(leg)p+theme(legend.position=c(.68,.55)) else p+theme(legend.position="none")}
pa<-mk("B alone",TRUE);pb<-mk("F alone");pc<-mk("B + F pair")
y<-d|>filter(case=="B + F pair",variable%in%c("Butyrate","Oxygen"))|>mutate(variable=recode(variable,Butyrate="Butyrate (mM)",Oxygen="Oxygen (a.u.)"))
pd<-ggplot(y,aes(time,value))+geom_line(linewidth=.7,colour="#CC79A7")+facet_wrap(~variable,scales="free_y",ncol=1,strip.position="left")+
  scale_x_continuous(breaks=c(0,100,200,300),expand=expansion(mult=c(0,.04)))+labs(title="Pair: butyrate rises, oxygen falls",x="Time (h)",y=NULL)+th+
  theme(strip.placement="outside",strip.background=element_blank(),strip.text=element_text(size=7))
fig<-(pa|pb)/(pc|pd)+plot_annotation(tag_levels="a")
ggsave("fig2_timecourse.png",fig,width=7.2,height=5,dpi=300);ggsave("fig2_timecourse.pdf",fig,width=7.2,height=5)
