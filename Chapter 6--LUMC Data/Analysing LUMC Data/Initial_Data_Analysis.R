source('Read_data.R')

data<-create_data()

################################################################################
#Covariate Transformation
################################################################################
transformations<-data.frame(NO=c(apply(data[,c(4:6)], 2, function(i) cor(i, data$severity_score_pseudo))),
                            ln=c(apply(log(data[,c(4:6)]), 2, function(i) cor(i, data$severity_score_pseudo))),
                            sqrt=c(apply(sqrt(data[,c(4:6)]), 2, function(i) cor(i, data$severity_score_pseudo))))
rownames(transformations)<-c('Cytokine 1', 'Cytokine 2', 'Cytokine 3')
transformations
#Log transformation gives largest correlation, thus we stick with log transformation
#Cytokine 1 has largest correlation with outcome. Thus, we first stick with cytokine 1

#Average number of days patient
total_days_per_patient<-data$day_no[c((match(unique(data$ID), data$ID)-1)[-1], length(data$ID))]-data$day_no[match(unique(data$ID), data$ID)]
mean(total_days_per_patient, na.rm=TRUE)
summary(total_days_per_patient)
#The average is 13, with a median of 10

#Figure 6.1 in thesis (page 39): The Spaghetti plot is created here
p0<-ggplot(data = data, aes(x = day_no, y = Cy_1, group = ID))+geom_line()+stat_smooth(aes(group = 1), method='loess', se=FALSE)+theme_bw()+
  labs(x = "Day", y='Cytokine 1')+theme(legend.text=element_text(size=16), text = element_text(size=12))
p1<-ggplot(data = data, aes(x = day_no, y = severity_score_pseudo, group = ID))+geom_line()+stat_smooth(aes(group = 1), method='loess', se=FALSE)+theme_bw()+
  labs(x = "Day", y='Severity Score')+theme(legend.text=element_text(size=16), text = element_text(size=12))
ggarrange(p1,p0, ncol=2)

################################################################################
#Time-dependency structure: Linear vs Linear Splines vs Quadratic vs Quadratic Splines
################################################################################
#We note that a linear profile with time does not seem to fit the data well. Thus, different approaches are tried

#Linear change over time
lme_linear<-lme(severity_score_pseudo ~day_no,
                random=  ~1| ID, 
                data=data, na.action=na.omit, method='ML')

#Piecewise linear change
#The mars function can tell us where the best breakpoints are located
data_full<-data[!is.na(data$day_no), ]
temp<-mars(data_full[c('day_no')],data_full[c('severity_score_pseudo', 'Cy_1')])
temp$cuts[temp$sel, ]
#We see the best breakpoint occurs at day 28
b1<-function(x, bp) ifelse(x<bp, 0, x-bp)
lme_piecewise<-lme(severity_score_pseudo ~day_no+b1(day_no, 28),
                   random=  ~1+day_no| ID, 
                   data=data, na.action=na.omit, method='ML')

#QUadratic change over time
lme_quadratic<-lme(severity_score_pseudo ~day_no+I(day_no^2),
                   random=  ~1| ID, 
                   data=data, na.action=na.omit, method='ML',control=lmeControl(msMaxIter = 200))

#Quadratic Spline models
lme_cubic<-lme(severity_score_pseudo ~day_no+b1(day_no, 17)+b1(day_no, 38)+I(day_no^2)+I(b1(day_no, 17)^2)+I(b1(day_no, 38)^2)+Cy_1,
               random=  ~1| ID, 
               data=data, na.action=na.omit, method='ML')

#We compare the different time-structures:
day_no_range<-c(min(data$day_no, na.rm=TRUE): 60) 
newdat <- data.frame(day_no = day_no_range, 
                     Cy_1=predict(loess(Cy_1~day_no, data=data), day_no_range))
data=data.frame(Days=day_no_range, Linear=predict(lme_linear, newdat, level=0),
                Piecewise_linear=predict(lme_piecewise, newdat, level=0),
                Quadratic=predict(lme_quadratic, newdat, level=0),
                Piecewise_Quadratic=predict(lme_cubic, newdat, level=0))
df <- data %>%
  gather(key = "Type_of_fit", value = "value", -Days)
# Obtain figure 6.2 in thesis: Comparison of Linear, Linear Spline, Quadratic and Quadratic Spline models
ggplot(df, aes(x = Days, y = value)) + 
  geom_line(aes(color = Type_of_fit)) + 
  scale_color_manual(values = c("red", "blue", 'green', 'black'))+
  xlab("Day since Infection") + ylab("Severity Score")+theme_bw()+
  theme(legend.text=element_text(size=16), text = element_text(size=16))

comparison_time_severity<-data.frame(Linear=unlist(summary(lme_linear)[c('AIC','BIC', 'logLik')]),
                                     Piecewise=unlist(summary(lme_piecewise)[c('AIC','BIC', 'logLik')]),
                                     Quadratic=unlist(summary(lme_quadratic)[c('AIC','BIC', 'logLik')]),
                                     Cubic=unlist(summary(lme_cubic)[c('AIC','BIC', 'logLik')]))
xtable(round(comparison_time_severity, 2))

#Piecewise and Quadratic appear to be the best models. Because of the easier interpretation of the piecewise model this model is preferred!

################################################################################
#Overall conclusions obtained during Data Analysis
################################################################################
#ln is the best transformation for all cytokines

#Since cytokine 1 has the largest correlation with the outcome we shall be using it first

#Time is not linear!!
#For Cy_1 Piecewise linear time-profile with breakpoint at day 28 is best