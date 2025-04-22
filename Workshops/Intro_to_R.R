### R Workshop
### Daniel Castaneda Mogollon
### University of Calgary
### October 4th, 2024

################################  BASICS OF R ################################ 
number1<-3
word1<-"Kananaskis"
vector1<-c(2,3,4,5,6)
vector2<-c(1,6,2,1,9)
vector3<-c("hello","this","is","daniel")

print(number1)
print(word1)
print(vector1)
print(vector2)
print(vector3)

sum_of_vectors1<-vector1+vector2
sum_of_vectors2<-vector1+vector3

mean_sum_of_vectors1<-mean(sum_of_vectors1)
print(mean_sum_of_vectors1)

matrix1<-as.matrix(sum_of_vectors1)
print(matrix1)

################ GENERATING RANDOM DATA AND INTRO TO NORMALITY ################ 
#Setting a seed for our random generated data
set.seed(123)
#Generating four data sets. Let's pretend these are gene expression values
data1 <-rnorm(n =100, mean = 50, sd=10)
data2 <-rnorm(n =100, mean = 54, sd=10)
data3 <-rnorm(n =100, mean = 3, sd = 15)
data4 <-sample(1:100, 100, replace=TRUE)
#First we want to know if the four genes have a normal distribution.
qqnorm(data1)
qqline(data1, col = "red")
qqnorm(data2)
qqline(data2, col = "red")
qqnorm(data3)
qqline(data3, col = "red")
qqnorm(data4)
qqline(data4, col = "red")
#If we see a diagonal 45 degree line, then it is most likely that the data is 
#normally distributed. But how can we be sure of this?
#Shapiro-Wilk normality test is a good start!
shapiro.test(x = data1)
shapiro.test(x = data2)
shapiro.test(x = data3)
shapiro.test(x = data4)
#The null hypothesis of a shapiro wilk test is that the population is normally 
#distributed.
#So what does it mean when you have a p-value < 0.05, is your data normally 
#distributed?

###############################  TESTING OUR DATA ###############################

#Let's test if the differences between genes is significant or not
t.test(x = data1, y = data2, alternative = "two.sided")
boxplot(x = data1)
t.test(x = data2, y = data3, alternative = "greater")
wilcox.test(x = data2, y = data4, alternative = "two.sided")

#Now let's move on to an ANOVA of the three data points that are normally distributed
gene_expression_data<-c(data1, data2, data3)
group<-factor(rep(c("Gene 1", "Gene 2","Gene 3"), each = 100))
anova_results<-aov(gene_expression_data ~ group)

print(anova_results)
summary(anova_results)

tukey_posthoc_comparison<-TukeyHSD(anova_results)
print(tukey_posthoc_comparison)

boxplot(gene_expression_data ~ group)
boxplot(gene_expression_data ~ group,
        main = "T1D gene expression",
        xlab = "Group of genes",
        ylab = "Expression levels",
        col = c("lightblue","pink","salmon"))

###No need to run this segment of the code, this is for visualization only ###
library(ggplot)
library("ggsignif")
df<-data.frame(group = group, gene_expression_data = gene_expression_data)
ggplot(df, aes(x = group, y = gene_expression_data, fill = group))+
  geom_boxplot() +
  geom_signif(comparisons = list(c("Gene 1","Gene 2"),
                                  c("Gene 1","Gene 3"),
                                  c("Gene 2","Gene 3")),
               map_signif_level = TRUE) +
  labs(title = "GGplot2 boxplot with Tukey Test Significance",
       x = "Genes",
       y = "Expression values") + theme_minimal()

######################EXERCISE##################################################
var1 = rnorm(1000, 70, 10)
var2 = rnorm(1000, 60, 25)
var3 = rnorm(1000, 95, 40)

#1.Are the means from var1 and var2 significantly different from one another? Why?
#2.What are the posthoc Tukey values of the three variables? Is one-way ANOVA the right test?
#3.Generate a boxplot of the three variables.

################################################################################
################################################################################
################################################################################

