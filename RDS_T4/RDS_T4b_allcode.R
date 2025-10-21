library(data.table)
########################################################################
#         R for Data Science- Programming Concepts                     #
#           How to think about and make sense of R code....            #
#                       Simon T RDS week 4                             #
#                                                                      #
########################################################################

# Two simple statements
s1 <- "test"
print(paste0("This is a ", s1))

# A Statement- a command in R that the interpreter executes to achieve some action
# s1 is a Variable- a variable can be seen a container to store data values
# Variables have two components:
#     -The variable Name- used to refer to the variable in code
#     -An object - physical space in memory that stores something and has a type
# The Name references the Object containing the data.

# If I write code like
"test2"

# The element with the value "test2" does not have a name, although it has a value
# At some point in execution "test2" o must exist as an object
# It cannot be referenced and unreferenced objects are soon deleted from the R memory
# by something called the Garbage Collector (GC)

# Consider just the paste0 part of the code
paste0("This is a ", s1)

# paste0 is a function, that here is passed two values, "This is a " and s1, 
# and returns a concatenated string.
# s1 is a variable that is created by the s1 <-"test" and once created will exist 
# in the workspace environment, sometimes even after the code has finished execution
# paste0 is a function [function (..., collapse = NULL, recycle0 = FALSE)] 
# that has several parameters set to default values (collapse and recycle0)
# we can access the help for this using
?paste0
# Paste0 also takes other parameters.
# The return value of the statement "paste0("This is a ",s1)" is not assigned to 
# a variable name and so will cease to exist almost as soon as the statement returns
# Inside the call to paste0, the string literal "This is a " is referenced as a 
# function parameter, and exists at least until the function returns

# We can even make statements from other statements (compound statements), 
# although the code starts to get messy!
print(paste0("This", paste0(" IS ", paste0("a ", "TEST"))))
# To make sense of this, consider the part- paste0("a ","TEST")
# - this just returns the character "a TEST"
# The second paste0 statement then reads paste0(" IS ","a TEST")
# -this just returns " IS a TEST"
# The third paste0 command is equal to paste0("This", " IS a TEST")
# - which returns "This IS a TEST"
# then the print statement takes this output and prints it on the screen.
# statements like
print(paste0("This", paste0(" IS ", paste0("a ", "TEST"))))
# Are poor coding style because there are easier ways to achieve the same result 
# that are much more readable eg
print("This IS a TEST")
# In this statement we do not compute anything from code to output
# In this statement, we just enter string literals (fixed value strings)!


# In R Strings are Immutable
# Immutable elements cannot be changed once they are created
# But I can change the value of string called s1 once it is created?

s1 <- "test- "

print(paste0(s1, address(s1)))

s1 <- "Another test- "

print(paste0(s1, address(s1)))

# So here we print out the memory address of the object referenced by the name s1.
# Note how the address changes when we change the value of s1
# So really what has happened is that we create a first object with the value called "test- "
# We then assign this to the name s1
# We then create a new object with the value "Another test- " and then assign this to s1
# The old value of s1 "test- " is no longer referenced and so sent to the GC

# We can also so something like this
s1 <- "test- "
s2 <- "Another test- "
s3 <- s1

print(paste0("s1 ", s1, address(s1)))
print(paste0("s2 ", s2, address(s2)))
print(paste0("s3 ", s3, address(s3)))

s3 <- "Yet another test- "
print(paste0("s1 ", s1, address(s1)))
print(paste0("s2 ", s2, address(s2)))
print(paste0("s3 ", s3, address(s3)))

# So now s3 has changed but really the object referenced by the name s3 has changed
# Note that the memory previously referenced by s3 still exists as s1- it does not go to the GC


# This is also how R manages to dynamically change object type that the name references
s1 <- "test- "
print(paste0("s1 ", s1, address(s1)))

s1 <- 420
print(paste0("s1 ", s1, " ", address(s1)))

# We are actually changing the object that the name references- the object itself 
# has an associated type!!

s1 <- 42
print(paste0("s1 ", s1, " ", address(s1)))

# Vectors are also interesting
vec1 <- c("t", "e", "s", "t")
print(paste0(address(vec1), " ", toString(vec1)))
print(paste0("vec1[1] ", vec1[1], " ", address(vec1[1])))
print(paste0("vec1[2] ", vec1[2], " ", address(vec1[2])))
print(paste0("vec1[3] ", vec1[3], " ", address(vec1[3])))
print(paste0("vec1[4] ", vec1[4], " ", address(vec1[4])))


vec1[2] <- "S"
print(paste0(address(vec1), " ", toString(vec1)))
print(paste0(address(vec1), " ", toString(vec1)))
print(paste0("vec1[1] ", vec1[1], " ", address(vec1[1])))
print(paste0("vec1[2] ", vec1[2], " ", address(vec1[2])))
print(paste0("vec1[3] ", vec1[3], " ", address(vec1[3])))
print(paste0("vec1[4] ", vec1[4], " ", address(vec1[4])))

# Notice all the memory addresses have changed- for the whole vector and all elements
# So actually make a single assignment to an element and every element is copied
# The vector is immutable
# Copying takes time!

# but what about lists
l1 <- as.list(vec1)
class(l1)

print(paste0("mylist: ", toString(unlist(l1)), " ", address(l1)))
print(paste0("l1[1] ", l1[[1]], " ", address(l1[[1]])))
print(paste0("l1[2] ", l1[[2]], " ", address(l1[[2]])))
print(paste0("l1[3] ", l1[[3]], " ", address(l1[[3]])))
print(paste0("l1[4] ", l1[[4]], " ", address(l1[[4]])))


l1[[3]] <- "P"


print(paste0("mylist: ", toString(unlist(l1)), " ", address(l1)))
print(paste0("l1[1] ", l1[[1]], " ", address(l1[[1]])))
print(paste0("l1[2] ", l1[[2]], " ", address(l1[[2]])))
print(paste0("l1[3] ", l1[[3]], " ", address(l1[[3]])))
print(paste0("l1[4] ", l1[[4]], " ", address(l1[[4]])))
# Maybe I should write a function for this??
# Notice that only the address changes for the single element that is updated
# This is because the list is like a cupboard, only elements swapped on shelves change.
# The list is mutable.


# Writing our own version of paste0- just to explore how it might work internally
# To concatenate two strings we exploit that we can convert a string to a vector
# and if we use c(vec1,vec2) the result is a concatenated vector
# We can then convert the result back to a string and return it
# The library function for pasteO is really the best way to concatenate strings
# in Base R (the basic R functionality that comes with the R install)

# My first attempt at a paste0 function
mypaste0 <- function(astring, bstring) {
  # create a character vector from each string
  tt <- strsplit(astring, split = "")
  yy <- strsplit(bstring, split = "")

  # concatenate the vectors-unlist forces the production of a simple vector
  zz <- unlist(c(tt, yy))

  # Note the return uses capture.output- this grabs something normally printed 
  # out and puts it back into a variable
  # Ends up as a vector type- cat almost manages the concatenation of strings 
  # directly but not quite, so we go through vectors
  return(capture.output(cat(zz, sep = "")))
}

# Test mypaste0
print(mypaste0("test of mypaste0 replacing ", "paste0"))

# But wait paste0 allows numerical values to be converted to strings as well- does this?
print(mypaste0("test of ", 32))
# No!

# Modify the function so that the values passed are converted to characters
# Conversion is always applied- even if they are already character type

mypaste0b <- function(aval, bval) {
  # Note the conversions here
  tt <- strsplit(as.character(aval), split = "")
  yy <- strsplit(as.character(bval), split = "")

  zz <- c(tt, yy)

  zz <- unlist(c(tt, yy))

  return(capture.output(cat(zz, sep = "")))
}

# Test mypaste0b
print(mypaste0b("test of mypaste0b replacing ", "paste0"))

# But wait paste0 allows numerical values to be converted to strings as well- does this?
print(mypaste0b("test of ", 32))
# Yes!




# But wait- for paste0 I can concatenate two or more elements- not just two
print(paste0("this ", "is ", " a", " test"))

# To code this I need to use an Ellipsis (...) in the function
# This means one or more parameters are passed eg one or more strings
# https://www.r-bloggers.com/2015/02/r-three-dots-ellipsis/

# more than two commands

mypaste0c <- function(...) {
  # This grabs the arguments as a list
  arguments <- list(...)

  # here strsplit is applied to the whole list at once
  yy <- strsplit(as.character(arguments), split = "")

  # zz is the whole list
  zz <- unlist(yy)

  return(capture.output(cat(zz, sep = "")))
}

# Test mypaste0c
print(mypaste0c("test of mypaste0c replacing ", "paste0"))
print(mypaste0c("test of ", 32))
print(mypaste0c("test of ", "mypaste0C replacing ", "paste0 ", "with number: ", 32))

# To make the strings look more readable I now decided to provide an option to
# capitalise the first letter of each work

mypaste0d <- function(..., first_letter_cap = FALSE) {
  arguments <- list(...)

  yy <- strsplit(as.character(arguments), split = "")


  if (first_letter_cap) {
    for (i in 1:length(yy)) {
      # first letter of each word in the list-
      yy[[i]][1] <- toupper(yy[[i]][1])
    }
  }

  zz <- unlist(yy)
  return(capture.output(cat(zz, sep = "")))
}

#
print(mypaste0d("test ", "mypaste0d ", "with ", "number: ", 32, first_letter_cap = TRUE))

# notice that I can't just write
print(mypaste0d("test ", "mypaste0d ", "with ", "number: ", 32, TRUE))
# The system does not know if TRUE is a part of the message or setting first_letter_cap
# So here we have to be explicit in passing this parameter to the function

# Finally a working set of functions.......
# Not quite.  Time to think through the code and try to consider possible ways
# a user would trigger an error?
# Then design the code to avoid this error-this is 'defensive coding'
# Here we might think about what happens if we call the function with only one parameter?
# Doing this seems pointless but it will be a common mistake when using the function.
# Users might also try to use our paste0 as a convenient way to convert any combination of parameters
# To an output for printing
# Ideally in this case we should pass through the single element as a character 
# that we return from our function

# eg
print(paste0("this is a test"))
print(mypaste0d("this is a test"))
print(mypaste0b("this is a test"))

# Happily our most advanced function (mypaste0d) handles this correctly but not mypaste0b
# Can you see why?  How might you tweak the input to make this function work as required??


# Overall here we have introduced some programming concepts (variables, objects, names, references, GC)
# We have looked at immutable and mutable objects
# We have looked at writing functions and using fancy methods such as using an Ellipsis
# Try to practice thinking through code in the way illustrated- it will help advance your coding skills.
# Avoid "Googling" or using LLMs to solve problems if you can- when you think though simple code you
# are learning skills that allow you to handle more complex cases.  If you offload your thinking
# to computer resources all the time you will struggle with more advanced coding tasks...



########################################################################
#         R for Data Science- R Questions                              #
#                                                                      #
#        Load this worksheet into RStudio and fill in the answers     #
#                                                                      #
########################################################################

# Question 0
# Write a function that prints "Hello World" one or more times
# Design this function so that it can be passed an integer that
# determines how many times to print out "Hello World"
# Provide an example of how to call this function
# Document the function
hello_world <- function(count) {
  rep(print("Hello World"), times=count)
}
rep(1:3, times=2)
hello_world(2)
hello_world(10)

# Question 1
myvector1 <- c(1, 2, 3, 5, 77, 4, 3, 5, 7, 5, -1)
# Write R code that multiples all elements by two
myvector1 * 2
# Write R code delete any elements that are <0
myvector1[myvector1 >= 0]
# Output the final result with values log base 2 transformed
log(myvector1, base = 2)


# Question 2
myvector2 <- c("apple", "orange", "mango", "lemon", "cabage", "carrot", "onion", "potato")

# What is the index position of "carrot" ? 2
# Replace "mango" with "pear"
for (i in 1:length(myvector2)) {
  if (myvector2[i] == "mango") {
    myvector2[i] = "pear"
  }
}
myvector2
# Fix the typo- change "cabage" to "cabbage"
for (i in 1:length(myvector2)) {
  if (myvector2[i] == "cabage") {
    myvector2[i] = "cabbage"
  }
}
myvector2
# Question 3a
myvector3a <- c(1, 2, 3, 4, 5, 6, 7, 45, 2, 3)
myvector3b <- c(1, 11, 3, 34, 5, 61, 7, 45, 2, 3)
myvector3c <- c(1, 2, 3, 4, 5, 6, 7, 45, 2, 3)
myvector3d <- c(1, 2, 3, 41, 5, 61, 7, 45, 21, 3)

# Write a program to calculate the median of each index position across the vectors
median_finder <- function(vec) {
  sort_vec <- sort(vec)
  vec_length <- length(sort_vec)
  if (vec_length %% 2 == 1) {
    index <- vec_length %/% 2
    median <- (sort_vec[index] + sort_vec[index + 1])/2
  } else {
    index <- vec_length / 2
    median <- sort_vec[index]
  }
  median
}
median_finder(myvector3a)
median_finder(myvector3b)
median_finder(myvector3c)
median_finder(myvector3d)

# Find the row index with the largest mean (hint consider using a data frame and apply
test_01 <- rbind(myvector3a, myvector3b, myvector3c, myvector3d)
means <- apply(test_01, 1 , mean)
test_02 <- cbind(test_01, means)
test_02
which(test_02[, "means"]==max(means), arr.ind=TRUE)


# Question 3b
# From the table of results from Question 3a combine these into a table
## Add the following header

myvector4c <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
colnames(test_01) <- myvector4c
test_01

# Write the table to a comma delimited file
write.csv(test_01, file="test_01.csv")

# Question 5
myvector5 <- c(1, 2, 3, 41, 5, 61, 7, 45, 21, 3, 7, 8, 11)
# Write a simple loop that will multiply the first element and then every second
# element (1st, 3rd, 5th etc) in the vector by 2 and then add 7
# Output the result of this transformation with each number presented as the 
# remainder after dividing each number by 7
for (i in seq(1, length(myvector5), by=2)) {
  myvector5[i] <- myvector5[i] * 2 + 7
  test_vector_5 <- myvector5 %% 7 
}
myvector5
test_vector_5


# Question 6
# You have a list of filenames
myvector6 <- c("file5", "file4", "file7", "File2", "file13", "file20", "file5", "file20")

# You would like all the names to begin with a capital "F"
# You want to append the year "_2021" to the end of each file name
for (i in 1:length(myvector6)) {
  substr(myvector6[i], 1, 1) <- "F"
  myvector6[i] <- paste0(myvector6[i], "_2021")
}
myvector6

# Question 7
mygenes <- c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10", "gene10", "gene9")
mychrs <- c("chr1", "chr3", "chr3", "chr4", "chr5", "chr1", "chr2", "chr3", "chr4", "chr7", "chr7", "chr4")

# make a vector of gene names and chromosomes that can be used to label a plot
my_merge <- mygenes
names(my_merge) <- mychrs
my_merge

# Question 7b
# Some of the labels made for Question 7 are not unique- write some R code that 
# fixes this by appending "_b" to the duplicates
dup_genes <- duplicated(mygenes)
dup_chrs <- duplicated(mychrs)
for (i in 1:length(mygenes)) {
  if (dup_genes[i]) {
    mygenes[i] <- paste0(mygenes[i], "_b")
  }
}
for (i in 1:length(mychrs)) {
  if (dup_chrs[i]) {
    mychrs[i] <- paste0(mychrs[i], "_b")
  }
}
mygenes
mychrs
my_merge2 <- mygenes
names(my_merge2) <- mychrs
my_merge2

# Question 7c
# Make a new table of all genes on chr2 or chr7 and write out a table 
# containing only these gene names to a file
my_list <- list(chr=mychrs, genes=mygenes)
my_mat <- array()
for (i in seq_along(my_list$chr)) {
  if (grepl("chr2|chr7", my_list$chr[i])) {
    vec <- c(my_list$chr[i], my_list$genes[i])
    my_mat <- cbind(my_mat, vec)
  }
}
my_mat
rownames(my_mat) <- seq_len(nrow(my_mat))
colnames(my_mat) <- seq_len(ncol(my_mat))
write.table(my_mat[, -1], file="my_mat.csv", row.names=FALSE, col.names=FALSE)


# Question 8 What is the result of multipying these vectors?

vec_v <- c(1, 2, 3, 4)
vec_w <- c(3, 4, 2, 4)
vec_result <- vec_v * vec_w
vec_result

# What do you expect the value of vec_result to be of the following multiplication?
# Why?

vec_vv <- c(1, 21, 3, 4) # numeric vector
vec_ww <- c("3", 4, 2, 4) # as "3" character is included, R auto coerces all elements in the vector to character type.
# So vec_ww is a character vector, not numeric.

vec_vv[4]
vec_ww[4]
vec_result <- vec_vv[4] * vec_ww[4]
# Error in vec_vv[4] * vec_ww[4] : non-numeric argument to binary operator

# Question 9
data(mtcars)
head(mtcars)

# These two commands grab some example data called "mtcars" that  describes the
# performance characteristics of different types of motor car.

# Make a data frame for the "mpg","hp" columns of this data (include all rows)
# Plot the values of mpg vs hp in a scatter plot.
# What do you conclude from this plot?

my_data <- mtcars[, c("mpg", "hp")]
plot(my_data)


# Question 10
# Go back over the code you have written and make sure it conforms to the style 
# guide and that follows the defensive coding approach
