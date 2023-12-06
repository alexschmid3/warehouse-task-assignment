
figuredirectory = string("outputs/combine")

#Get list of files
folderlist = readdir(figuredirectory)

#Get the file types to combine
filefirstwords = [file[1:3] for file in folderlist]
filetypes = unique!(filefirstwords)
dataDict = Dict()
for t in filetypes
	dataDict[t] = []
end
for file in folderlist
	filetype = file[1:3]
	push!(dataDict[filetype], file)
end

#Combine each file type to one csv
for t in ["run"] #filetypes
	
	#Name of combo file
	combofile = string(figuredirectory, "/", t,"_combined.csv")

	#Combine files and save to combofile
	for file in dataDict[t]
		#Add all rows (including header) from first file
		if file == dataDict[t][1]
			try
				open(string(figuredirectory, "/", file, "/output.csv")) do input
					open(combofile, "a") do output
						for line in eachline(input)
							println(output, line)
						end
					end
				end
				println("Added ", file)
			catch
				println("No output for ", file)
			end
		#Drop header row from other files
		else
			try
				open(string(figuredirectory, "/", file, "/output.csv")) do input
					open(combofile, "a") do output
						for line in Iterators.drop(eachline(input), 1)
							println(output, line)
						end
					end
				end
				println("Added ", file)
			catch 
				println("No output for ", file)
			end
		end
	end

end 
