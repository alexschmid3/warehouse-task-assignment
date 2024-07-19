
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if item[2][1] + item[2][2] + item[2][3] <= 615
		push!(hueList, item) 
	end
end

#drawingname, arclistlist, colorlist, thicknesslist, fractlist, x_size, y_size = string("outputs/viz/order", i,".png"), arclistlist, colorlist, thicknesslist, fractlist, 2000, 1200

function timespacenetwork(drawingname, arclistlist, colorlist, thicknesslist, dashlist, fractlist, x_size, y_size)

	#Find coordinates for each time-space node
	nodelist = []
	x_size_trimmed, y_size_trimmed = x_size*0.9, y_size*0.9
	k1 = x_size_trimmed/(horizon/tstep + 2) 
	k2 = y_size_trimmed/(numlocs + 2)
	for i in 1:numnodes
		ycoord = nodelookup[i][1]
		xcoord = (nodelookup[i][2]/tstep)+1

		#Scaling to image size
		tup = (-x_size_trimmed/2 + xcoord*k1, -y_size_trimmed/2 + ycoord*k2)   
		
		push!(nodelist,tup)
	end

	#Create actual points as a Luxor object
	nodePoints = Point.(nodelist)

	#---------------------------------------------------------------------------------------#

	#Arcs for visualization
	#Duplicate for multiple input arc lists with different colors/thickness/dash if you're trying to show m
	arcinfo = []
    for j in 1:length(arclistlist), a in intersect(1:numarcs,arclistlist[j])
        startPoint = nodePoints[arclookup[a][1]]
        endPoint = nodePoints[arclookup[a][2]]
        
        #Set arc attributes
        arcColor = colorlist[j] # (0,0,255) #RGB tuple 
        arcDash = dashlist[j] #"solid" #"solid", "dashed"			
        arcThickness = thicknesslist[j]
		#if fractlist[j] != []
		#	arcLabel = fractlist[j][a]
		#else
		arcLabel = ""
		#end
        
        #Add to arcinfo list to be used in the drawing 
        push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness, arcLabel))
    end

	#-------------------------------------------------------------------------#

	#Initiailize drawing
	Drawing(x_size, y_size, drawingname)
	origin()
	background("white")

	#Draw arcs
	for i in arcinfo #union(arcinfo,legendarcs)
		
		#Set arc attributes from the arcinfo
		r_val, g_val, b_val = i[3][1]/255, i[3][2]/255, i[3][3]/255
		setcolor(convert(Colors.HSV, Colors.RGB(r_val, g_val, b_val)))  #You can also use setcolor("colorname")
		setdash(i[4])
		setline(i[5])

		#Draw the line from the start node to end node
		Luxor.line(i[1], i[2] , :stroke)
		
		#Figure out the angle of the arrow head
		theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
		dist = distance(i[1], i[2])
		arrowhead = (1-12/dist)*i[2] + (12/dist)*i[1] #12 pixels from the end node
		local p = ngon(arrowhead, max(10, i[5]), 3, theta, vertices=true)
		
		#Draw the arrow head
		poly(p, :fill,  close=true)

	end

	#Draw node points
	setcolor("black")
	Luxor.circle.(nodePoints, 9, :fill)

	#Set font size for labels
	fontsize(40)

	#Add location labels
	for l in 1:numlocs
		coord = nodePoints[nodes[(l,0.0)]]
		label("Location $l       ", :W , coord)
	end
	#for l in 13:14
	#	coord = nodePoints[nodes[(l,0.0)]]
	#	label("Workstation $l       ", :W , coord)
	#end

	#Add time labels
	for t in 0:tstep*2:horizon
		coord = nodePoints[nodes[(1,t)]] + Point(0,-30)
		label("t = $t", :N , coord)
	end

	#Add arc labels
	for i in arcinfo
		if i[6] != ""
			fontsize(20)
			coord = Point((i[1][1] + i[2][1])/2, (i[1][2] + i[2][2])/2)
			sethue("black")
			Luxor.rect(coord+Point(-30,-25), 60, 25, :fill)
			sethue("white")
			label(i[6], :N , coord)
		end
	end

	#Legend labels
	#fontsize(60)
	#Luxor.text("Standard journey", Point(-500,800), halign=:left, valign = :middle)
	#Luxor.text("Extended journey", Point(400,800), halign=:left, valign = :middle)

	finish()
	preview()

end
