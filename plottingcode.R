library(TMB)
library(tidyverse)
library(oneTimeStepPredict)
library(numDeriv)
library(lmtest)
library(reshape2)
library(cowplot)


Nay_p <- function(dat){

    p11=ggplot({{dat}}, aes(x=year,y=YEhat))+
        labs(x = "Year",y = expression(hat(N)[y])) + 
        geom_ribbon(aes(ymin = L_ran,ymax = U_ran),fill = "black",alpha=0.5) +  
        geom_line(aes(y = YE),size=1,color='blue')  +  
        geom_line(aes(y = YEhat),size=1,color='black')+ 
        theme(axis.title.x=element_blank())

    p11
}


resid_plot <- function(dat,column,std_res_col,clr,label){
    p21 = ggplot({{dat}}, aes(x=year,y={{std_res_col}},color=age))+
        labs(x = "Year",y = 'Residual') + 
        geom_smooth(formula = y ~ x,span=0.1,alpha=0.25,color='black',se=F)+
        geom_point(aes(y = {{column}},color=label),size=1,alpha=0.5)+  
        geom_smooth(aes(y = {{column}}),formula = y ~ x,span=0.1,fill=clr,alpha=0.25,color=clr)+
        geom_hline(yintercept = 0,linetype='dashed')+           
        scale_color_manual(values=c(clr)) +
        guides(colour =  guide_legend(position = "inside"))+
        theme(legend.position.inside = c(0.1, 0.95),legend.title=element_blank(),
              legend.text = element_text(size=6),legend.key.size = unit(0.25, "cm"), 
              legend.background=element_blank(),
              axis.title.x=element_blank())
    p21
}

resid_comp_p <- function(dat,column1,column2,column3,column4,column5,labels,colors){
    p41 = ggplot({{dat}}, aes(x={{column1}},y={{column2}},color=age))+
        labs(x = "Pearson Residual",y = 'Predict Residual') +          
        geom_point(aes(y = {{column2}}, color=labels[2]),size=1,alpha=0.5)+ 
        geom_point(aes(y = {{column3}}, color=labels[3]),size=1,alpha=0.5)+  
        geom_point(aes(y = {{column4}}, color=labels[4]),size=1,alpha=0.5)+ 
        geom_point(aes(y = {{column5}}, color=labels[5]),size=1,alpha=0.5)+ 
        geom_abline(intercept = 0,slope=1,linetype='dashed')+           
        scale_color_manual(values=c(colors[4],colors[3],colors[1],colors[2])) +
        guides(colour =  guide_legend(position = "inside"))+
        theme(legend.position.inside = c(0.1, 0.8),legend.title=element_blank(), 
              legend.background=element_blank(),
              legend.text = element_text(size=6),legend.key.size = unit(0.25, "cm"))
    p41
}


dens_p <- function(dat,res_sd_mean,colors){
    p42 = ggplot({{dat}},aes(x = value, color = variable)) +
        geom_density(alpha = 0.4,linewidth=1)+
        geom_boxplot(aes(x = value, color = variable, group=variable,y=box,fill=variable),
                     color = "black", alpha = 0.25,notch="TRUE") +        
        scale_color_manual(values=c(colors[1],colors[2],colors[3],colors[4],colors[5]))      +        
        scale_fill_manual(values=c(colors[1],colors[2],colors[3],colors[4],colors[5]))+ 
        guides(fill="none",color="none") +
        theme(axis.title.y=element_blank()) 
    
    p42 = p42+ geom_text(data={{res_sd_mean}}, aes(x=3.3, y=y, label=label),color='black',size=2)
    p42
}

fix_d <- function(dat,resid_col_names){
            sdat = {{dat}}[,c("year","age",resid_col_names)]    
            sdat1=melt(sdat,id.vars=c("year","age"))
            
            sdat1$box = rep(c(0,0.1,0.2,0.3,0.4),each=200)

            resid.sd= aggregate(data=sdat1,value ~ variable, sd)
            resid.means = aggregate(data=sdat1,value ~ variable + year, mean)
            resid.sd_mean = aggregate(data=resid.means,value ~ variable, sd)
            resid.sd_mean$y = c(0,0.1,0.2,0.3,0.4)+0.02
            resid.sd_mean$label = round(resid.sd_mean$value,2)
            ret = list(sdat1,resid.sd_mean=resid.sd_mean,resid.sd=resid.sd)
}

bubble_plot <- function(dat,column,pvcolumn,dwtest.dat,name){
    datums = select({{dat}},everything())
    datums = mutate(datums,size=abs({{column}}))
    datums = mutate(datums,clr=ifelse({{column}} < 0,'-','+'))

    ggplot(datums,aes(x=year, y=age, size=size, color=clr)) +
        geom_point(alpha=0.5) +       
                                        #  facet_wrap(~survey, ncol = 2) +                     
        scale_size_continuous(range = c(.05, 5), name="Pearson Resid")  +   
        theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
              legend.position = "top")+
        scale_color_manual(values=c("blue","red"))+
        scale_y_continuous(breaks=1:A) + guides(color = "none", size="none") +
        theme(axis.title.x=element_blank()) +
        labs(y = name) + 
        xlim(1,T+5)+ 
        geom_text(data={{dwtest.dat}}, aes(x=T+4, y=age,size=1,label = round({{pvcolumn}}, digits=2)),color='black')

}

qq_p <- function(dat,column,clr,pvalue){
    p42=ggplot({{dat}}, aes(sample = {{column}})) + stat_qq(color=clr) + stat_qq_line(linetype='dashed')+
        theme(axis.title.y=element_blank(),axis.title.x=element_blank())+
        geom_text(aes(x = -2, y = 2, label = round(pvalue, 2)))
    p42
}

resid_plot2 <- function(dat,column,clr,clr2){
    p12 = ggplot(dat, aes(x=year,y={{column}},color=fsurvey))+
        labs(x = "Year",y = 'Residual') + 
        geom_point(aes(y = {{column}},color=fsurvey),size=1,alpha=0.5)+  
        geom_smooth(aes(y = {{column}},color=fsurvey,group=fsurvey,fill=fsurvey),formula = y ~ x,span=0.1,alpha=0.25)+
        geom_hline(yintercept = 0,linetype='dashed')+           
        scale_color_manual(values=c(clr,clr2)) +         
        scale_fill_manual(values=c(clr,clr2)) +
        guides(colour =  guide_legend(position = "inside"))+
        theme(legend.position.inside = c(0.1, 0.95),legend.title=element_blank(),
              legend.text = element_text(size=6),legend.key.size = unit(0.25, "cm"), 
              legend.background=element_blank(),
              axis.title.x=element_blank(),axis.title.y=element_blank())
}


vio_m_plot <- function(dat,fyear,value,fact,f_lab,x_lab,y_lab){
    p1=ggplot(data={{dat}}, aes(x={{fyear}}, y={{value}}, fill={{fact}})) +
    facet_wrap(~method, ncol=1)+     
    geom_hline(yintercept = 0,linetype='dotted')+           
    geom_violin()+   
    theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
    legend.position = "top") + labs(fill=f_lab,x=x_lab,y=y_lab)+
    geom_line(data = mdat,aes(x={{fyear}}, y={{value}},group = {{fact}})) 

    p1
}


##The colors
clr1 ='orange' 
clr2 ='blue'
clr3 ='burlywood3'  
clr4 ='darkolivegreen4'

