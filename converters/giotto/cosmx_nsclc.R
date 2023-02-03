# load the giotto file
load("SMI_Giotto_Object.RData")

# prints the variables found in memory; this tells us that one variable called "gem" has been loaded
ls()

# this tells us that "gem" is a Giotto object
class(gem)

# let's see what's inside the Giotto object 
slotNames(gem)

# we can access a slot with
gem@`parameters`

# if we try to access "custom_expr" we can't because the slot is not there
gem@`custom_expr`

# hacky line to check which slots are available
for(v in slotNames(gem)) {
	tryCatch({
		content = eval(str2expression(sprintf("gem@`%s`", v)))
		print(paste(v, "FOUND"))
		print(eval(str2expression(sprintf("names(gem@`%s`)", v))))
		# this commented line prints everything
		# print(content)
	}, error=function(error_message){
		print(paste(v, "not present"))
	})
	print("************************************")
}

# example, let's get the umap data. We get NULL (there is actually nothing there)
gem@`parameters`$`7_umap`

# example, with something inside
gem@`cell_metadata`

