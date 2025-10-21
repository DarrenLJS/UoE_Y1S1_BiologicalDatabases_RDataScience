##################################################################
#         R for Data Science- R Answers   Class 4          #
##################################################################

#Question 0
#Write a function that prints "Hello World" 
#Design this function so that it can be passed an integer that 
#determines how many times to print out "Hello World"
#Provide an example of how to call this function

#We might also choose to put all functions in a separate file
#then source("myfile.R") to load the functions

print_hello_world <-function(val){
  
  
  if(getOption("warn") >0){
    stopifnot(is.numeric(val),val>0)
  }
  
  for(i in 1:val) print("Hello World")
}

#example call the function
print_hello_world(5)

#Question 1
myvector1 <-c(1,2,3,5,77,4,3,5,7,5,-1)
#Write R code that multiples all elements by two
myvector1 <-myvector1*2

#Write R code delete any elements that are <0
myvector1 <-myvector1[myvector1 >=0]

#Output the final result with values log transformed
myvector1 <-log(myvector1)


#Question 2
myvector2 <-c("apple", "orange", "mango" ,"lemon","cabage", "carrot", "onion", "potato")

#What is the index position of "carrot" ?
#Of course we could also loop through the vector until we find mango
pos <-which(myvector2 == "mango")

#Replace "mango" with "pear"
myvector2[pos]<-"pear"

#Fix the typo- change "cabage" to "cabbage"
#You could also edit the input!!
pos <-which(myvector2 == "cabage")
myvector2[pos]<-"cabbage"


#Question 3a
myvector3a <-c(1,2,3,4,5,6,7,45,2,3)
myvector3b <-c(1,11,3,34,5,61,7,45,2,3)
myvector3c <-c(1,2,3,4,5,6,7,45,2,3)
myvector3d <-c(1,2,3,41,5,61,7,45,21,3)

#Write a program to calculate the median of each index position across the vectors
#Find the row index with the largest mean (hint consider using a data frame and apply

#make  a data frame and convert to a matrix
mytable <- data.frame(myvector3a,myvector3b,myvector3c,myvector3d)
mytable <-as.matrix(mytable)

#medians using apply "1" signifies apply per row
mytable_medians <-apply(mytable,1,median)

#also means
mytable_means <-apply(mytable,1,mean)

#max mean
mytable_means_max <-max(mytable_means)

#max_mean_index
which(mytable_means==mytable_means_max)

#The question was a bit ambiguous -it was really was asking for column means and medians so...
#medians using apply "2" signifies apply per column
mytable_col_medians <-apply(mytable,2,median)

#also means
mytable_col_means <-apply(mytable,2,mean)

#etc etc




#Question 3b
#From the table of results from Question 3a combine these into a table
##Add the following header

#a bit of an ambiguous question- expected answer would be, but truncation of the header or 
#addition of extra columns would also be OK!

myvector4c <-c("A","B","C","D","E","F","G","H","I","J")

rownames(mytable)<-myvector4c

#Write the table to a comma delimited file
write.csv(mytable,file="mytable.csv")

#Question 5
myvector5 <-c(1,2,3,41,5,61,7,45,21,3,7,8,11)
#Write a simple loop that will multiply the first element and then every second element in the vector by 2 and then add 7
#Output the result of this transformation with each number presented as the remainder after dviding each number by 7


#create an empty vector to hold the results
res <-vector(mode="numeric",length=length(myvector5))

#build the loop- lots of possible ways of doing this!
for (count in 0:length(myvector5)-1) {

  if(count%%2==0){
    res[count+1]=(myvector5[count+1]*2)+7

  }else{
    res[count+1]=myvector5[count+1]
  }
  count=count+1
}

#Dump this
res%%7


#Question 6
#You have a list of filenames
myvector6 <-c("file5", "file4", "file7" ,"File2","file13", "file20", "file5", "file20")

#You would like all the names to begin with a capital "F" rather than lowercase
#This can be solved with a loop or using command gsub
myvector6 <-gsub("^f", "F", myvector6)

#You want to append the year "_2021" to the end of each file name
myvector6 <-paste0(myvector6,"_2021")

#Question 7
mygenes <-c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9","gene10","gene10","gene9")
mychrs <-c("chr1","chr3","chr3","chr4","chr5","chr1","chr2","chr3","chr4","chr7","chr7","chr4")

#make a vector of gene names and chromosomes that can be used to label a plot
#We assume that the index position in mygenes is a gene name and the same index in mychrs is the chromosome
myvector7 <-paste(mygenes,mychrs,sep="_")

#Question 7b
#Some of the labels made for Question 7 are not unique- write some R code that fixes this by appending "_b" to the duplicates
#Note this is a weak solution- it only deals with duplicates- we could append a unique number to fix this...

dups<-duplicated(myvector7)

for(i in 1:length(myvector7)){
  if(dups[i])myvector7[i]<-paste0(myvector7[i],"_b")
}

#There is also the command make.unique, which can be applied to a character vector- so no loop required!

#Also consider- just because the name is duplicated is the data the same across the duplicated data row- so should we rename
#or remove duplicates??


#Question 7c
#Make a new table of all genes on chr2 or chr7 and write out a table containing only these gene names to a file
selected_chr2_chr7 <-which(mychrs=="chr7")
selected_chr2_chr7  <-c(selected_chr2_chr7,which(mychrs=="chr2"))
mygenes <-mygenes[selected_chr2_chr7]

#Question 8 results multiplication two vectors- question was not very clear!
vec_v <-c(1,2,3,4)
vec_w <-c(3,4,2,4)

vec_vw <-vec_v *vec_w
#element wise multiplication


vec_vv <-c(1,21,3,4)
vec_ww <-c("3",4,2,4)

vec_result <-vec_vv[4] * vec_ww[4]

#This restults in an error
#The reason is that vectora can only hold one type ww has a mix- so it's type is converted to a vector of strings
> vec_ww
[1] "3" "4" "2" "4"

#Then the multiplication fails because you cannot multiply a number by  a string
vec_result <-vec_vv[4] * as.integer(vec_ww[4])

#this works with an explicit type conversion 4*4
#> vec_result
#[1] 16

#Question 9
#These two commands grab some example data called "mtcars" that  describes the
#performance characteristics of different types of motor car engine.

#Make a data frame for the "mpg","hp" columns of this data (include all rows)
#Plot the values of mpg vs hp in a scatter plot.
#What do you conclude from this plot?

data(mtcars)
head(mtcars)

#head(mtcars)
#                   mpg cyl disp  hp drat    wt  qsec vs am gear carb
#Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1

mydata <-mtcars[,c("mpg","hp")]
plot(mydata)

#Avoid high HP cars is you want to kepp petrol costs lowercase
#but the relationship is not linear so you can choose a car with
#the lowest HP if you want good fuel economy (something like this)



#Question 10  
#Go back over the code you have written and make sure it conforms to the style guide and that follows the defensive coding approach 
#There is an extea sheet with some examples for this...