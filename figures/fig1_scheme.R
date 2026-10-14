suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
nodes<-tribble(~x,~y,~lab,~fill,
 1.0,4.5,"Sugar","#EEEEEE", 3.6,4.5,"B\nBifidobacterium","#CFE3F3", 6.0,4.5,"Acetate","#EEEEEE", 8.4,4.5,"F\nFaecalibacterium","#CDEBDD",
 10.9,4.5,"Host\ncarbon","#EEEEEE", 8.4,2.8,"Butyrate","#EEEEEE", 8.4,1.1,"Host colon wall","#E8DCEB",
 5.2,1.1,"Oxygen","#EEEEEE", 2.0,1.1,"P\npathogen","#F6D6C5")
A<-function(x,y,xe,ye,lab=NULL,lx=NULL,ly=NULL,dash=FALSE,bar=FALSE,col="grey25") data.frame(x,y,xe,ye,lab=lab%||%NA,lx=lx%||%NA,ly=ly%||%NA,dash=dash,bar=bar,col=col)
`%||%`<-function(a,b) if(is.null(a)) b else a
e<-bind_rows(
 A(1.8,4.5,2.7,4.5), A(4.5,4.5,5.3,4.5,"releases",4.9,4.85), A(6.7,4.5,7.4,4.5,"uses",7.05,4.85),
 A(10.2,4.5,9.4,4.5), A(8.4,4.1,8.4,3.2,"butyrate\nformed",9.2,3.65), A(8.4,2.4,8.4,1.5,"fuels",8.95,1.95),
 A(7.3,1.1,5.9,1.1,"uses O2",6.6,1.38,bar=TRUE,col="#7B3F98"),
 A(4.5,1.25,2.8,1.25,"respiration",3.65,1.55,col="#D55E00"),
 A(2.7,0.95,4.5,0.95,"inflammation\nraises O2",3.6,0.45,col="#D55E00"),
 A(1.0,4.1,1.7,1.6,"shared sugar",0.55,2.9,col="#D55E00"),
 A(4.9,1.6,3.5,4.1,"O2 harms\nB and F",3.55,2.85,dash=TRUE,bar=TRUE,col="grey40"),
 A(5.8,1.6,7.7,4.1,NULL,dash=TRUE,bar=TRUE,col="grey40"))
p<-ggplot()+
 geom_segment(data=filter(e,!bar),aes(x=x,y=y,xend=xe,yend=ye,colour=col,linetype=dash),arrow=arrow(length=unit(.16,"cm"),type="closed"),linewidth=.5)+
 geom_segment(data=filter(e,bar),aes(x=x,y=y,xend=xe,yend=ye,colour=col,linetype=dash),arrow=arrow(length=unit(.12,"cm"),angle=90),linewidth=.5)+
 geom_text(data=filter(e,!is.na(lab)),aes(lx,ly,label=lab),size=2.2,colour="grey25",lineheight=.9)+
 geom_label(data=nodes,aes(x,y,label=lab,fill=fill),size=2.6,label.size=.3,label.padding=unit(.28,"lines"),lineheight=.9)+
 scale_fill_identity()+scale_colour_identity()+scale_linetype_manual(values=c(`FALSE`="solid",`TRUE`="22"),guide="none")+
 annotate("text",x=0.1,y=5.7,label="Lumen",hjust=0,size=3,fontface="bold",colour="grey40")+
 annotate("text",x=11.5,y=0.35,label="Mucosa",hjust=1,size=3,fontface="bold",colour="grey40")+
 coord_cartesian(xlim=c(0,11.6),ylim=c(0.1,5.9))+theme_void()
ggsave("fig1_scheme.png",p,width=7.2,height=4,dpi=300,bg="white");ggsave("fig1_scheme.pdf",p,width=7.2,height=4)
