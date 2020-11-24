# https://stackoverflow.com/questions/22531477/using-lists-inside-data-table-columns
# http://brooksandrew.github.io/simpleblog/articles/advanced-data-table/
# https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
# https://r4ds.had.co.nz/many-models.html
# https://www.kaggle.com/timib1203/managing-many-models-with-tidyr-purrr-and-broom

dt = data.table(id = 1:2, comment = vector("list", 2L))

# assign value 1 to just the first column of 'comment'
dt[1L, comment := 1L]

# assign value of 1 and "a" to rows 1 and 2
dt[, comment := list(1, "a")]

# assign value of "a","b" to row 1, and 1 to row 2 for 'comment'
dt[, comment := list(c("a", "b"), 1)]

# assign list(1, "a") to just 1 row of 'comment'
dt[1L, comment := list(list(list(1, "a")))]

dt = data.table(id = 1:2, comment = vector("list", 2L))

# assign value 1 to just the first column of 'comment'
set(dt, i=1L, j="comment", value=1L)

# assign value of 1 and "a" to rows 1 and 2
set(dt, j="comment", value=list(1, "a"))

# assign value of "a","b" to row 1, and 1 to row 2 for 'comment'
set(dt, j="comment", value=list(c("a", "b"), 1))

# assign list(1, "a") to just 1 row of 'comment'
set(dt, i=1L, j="comment", value=list(list(list(1, "a"))))


dt <- data.table(mtcars)[, .(cyl, gear)]
dt[,unique(gear), by=cyl]


dt <- data.table(mtcars)[,.(gear, cyl)]
dt[,gearsL:=list(list(unique(gear))), by=cyl] # original, ugly
dt[,gearsL:=.(list(unique(gear))), by=cyl] # improved, pretty
head(dt)

dt[,gearL1:=lapply(gearsL, function(x) x[2])]
dt[,gearS1:=sapply(gearsL, function(x) x[2])] 

head(dt)

str(head(dt[,gearL1])) 

str(head(dt[,gearS1]))

dt[,gearL1:=lapply(gearsL, `[`, 2)]
dt[,gearS1:=sapply(gearsL, `[`, 2)]

dt[,other_gear:=mapply(function(x, y) setdiff(x, y), x=gearsL, y=gear)]
head(dt)

dt[,other_gear:=mapply(setdiff, gearsL, gear)]

