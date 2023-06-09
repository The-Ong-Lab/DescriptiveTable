#will need to load table functions and its dependent libraries
library(...)

#load data
df <- mtcars

#set classes
df$cyl <- as.factor(df$cyl)
df$vs <- as.factor(df$vs)
df$am <- as.factor(df$am)

#make table
tbl <- combined_tbl('cyl', c('mpg','hp','wt','vs','am','gear'), df)
tbl

#clean up names
##rows
tbl[1] <- c('Total','Miles/Gallon', 'Horsepower', 'Weight (lb/1000)', 
            'Straight Line Cylinder Configuration', 'Manual Transmition',
            'Number of Gears')
##cols
colnames(tbl) <- c('Variables','Count(%)/Mean(SD)',
                   'Four','Six','Eight',
                   'p-value') 

#a kinda plain table with kable
tbl %>% 
  kbl() %>% #pipe to a knittable table (aka, kable)
  kable_classic(c('striped')) #add some predefined styling and stripes
  
#add some more features and formatting
tbl %>% kbl() %>% kable_classic(c('striped')) %>% 
  #add headers above the columns
  add_header_above(c(" "=1, "Total Cohort" = 1, "Number of Cylinders" = 3, " "=1)) %>% 
  #add a split between the full counts and the strata counts
  column_spec(3, border_left = TRUE) %>% 
  #group some rows together
  pack_rows("Group these rows together", 4, 5, label_row_c = "border-bottom: 0px solid;") %>% 
  #add some indentations to the mpg and weight rows
  add_indent(2) %>% add_indent(4)
